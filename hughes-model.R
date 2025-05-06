library(tidyverse)

# dij is distance between i and j
# dxbar is average dispersal distance of organism
# distances is all dij for one point. AKA ith row of distance matrix
# delta_x is a scaling factor, it represents the area of one grid cell
# D is the distance matrix
# H / P is the herbivore / parasitoid density after growth step
# H_star / P_star is the herbivore / parasitoid density after dispersal
# lambda is herbivore intrinsic growth rate
# a is parasitoid search rate
# I is a habitat suitability vector
# phi is parasitoid phenology. phi = 0 means P escape H density dependent mortality

# 10x10 tester
load(file = "data/10x10_landscape_info.RData")
load(file = "data/10x10_dist_matrix.RData")
D <- gp_movement_matrix
omega <- nrow(D)
landscape_dataframe <- gp_landscape
habitat_suitability <- 1-round(landscape_dataframe$value)

# Dispersal kernel ####
# 1. movement kernel
m_kernel <- function(dij, dxbar) { 
  2 / (pi * dxbar^2) * exp(-2 * dij / dxbar)
}

# 2. Cut‑off distance d_max(dxbar)
d_max <- function(dxbar) {
  -0.5 * dxbar * log(pi * dxbar^2 * 1e-10 / 2)
}

# 3. Single‐entry of normalized kernel k_ij of species x
kernel_element <- function(dij, distances, dxbar, delta_x) {
  dm   <- d_max(dxbar)
  if (dij > dm) return(0)
  num  <- m_kernel(dij, dxbar) * delta_x^2
  denom <- sum(m_kernel(distances[distances <= dm], dxbar) * delta_x^2)
  num / denom
}

# 4. Full kernel matrix K (Ω×Ω) for species x
compute_kernel_matrix <- function(distance_matrix, dxbar, delta_x) {
  Ω <- nrow(distance_matrix)
  K <- matrix(0, nrow=Ω, ncol=Ω)
  for (i in seq_len(Ω)) {
    d_row <- distance_matrix[i, ]
    for (j in seq_len(Ω)) {
      K[i,j] <- kernel_element(dij      = d_row[j],
                               distances = d_row,
                               dxbar     = dxbar,
                               delta_x   = delta_x)
    }
  }
  K
}


# Compute dispersal ####
# H* = K_H %*% H
compute_H_star <- function(H, K_H) {
  as.vector(K_H %*% H)
}

# P* = K_P %*% P
compute_P_star <- function(P, K_P) {
  as.vector(K_P %*% P)
}

# Herbivore growth ####
# 1. Density‐dependence: exp(−log(λ)·H*)
density_regulation <- function(H_star, lambda) {
  exp(2 * -log(lambda) * H_star)
}

# 2. Survival after parasitism: exp(−a·P*)
host_parasitism_survival <- function(P_star, a = 1) {
  exp(-a * P_star)
}

# 3. Full H‐update
update_H <- function(H_star, P_star, habitat_suitability, lambda, a = 1) {
  habitat_suitability * lambda * H_star * 
    density_regulation(H_star, lambda) * 
    host_parasitism_survival(P_star, a)
}

# Parasitoid growth ####
# 1. Functional response: 1 − exp(−a·P*)
parasitoid_attack <- function(P_star, a = 1) {
  1 - exp(-a * P_star)
}

# 2. Regulation via host‐density feedback: exp(−ϕ·log(λ)·H*)
parasitoid_regulation <- function(H_star, lambda, phi = 0) {
  exp(-phi * log(lambda) * H_star)
}

# 3. Parasitoid recruitment rate: xi = a * b * K
compute_xi <- function(a, b, K) {
  xi <- a * b * K
  return(xi)
}

# 4. Full P‐update
update_P <- function(H_star, P_star, habitat_suitability, xi, lambda, a = 1, phi = 0) {
  habitat_suitability * xi * H_star *
    parasitoid_attack(P_star, a) *
    parasitoid_regulation(H_star, lambda, phi)
}

# Step simulator ####
step_host_parasitoid <- function(H, P, # densities
                                 dist_mat, # distance matrix
                                 dxbar_H, dxbar_P, # avg dispersal distances
                                 delta_x, # pixel scale
                                 habitat_suitability, # habitat suitability
                                 lambda,       # host growth
                                 a_H  = 1,     # host parasitism coefficient
                                 xi,           # parasitoid recruitment
                                 a_P  = 1,     # parasitoid attack
                                 phi  = 0      # feedback strength
) {
  # build kernels
  K_H <- compute_kernel_matrix(dist_mat, dxbar_H, delta_x)
  K_P <- compute_kernel_matrix(dist_mat, dxbar_P, delta_x)
  
  # dispersal
  Hs <- compute_H_star(H, K_H)
  Ps <- compute_P_star(P, K_P)
  
  # local growth
  H_next <- update_H(Hs, Ps, habitat_suitability, lambda, a_H)
  P_next <- update_P(Hs, Ps, habitat_suitability, xi,     lambda, a_P, phi)
  
  list(H = H_next, P = P_next)
}


H_init <- runif(omega, min = 0, max = 1)   # e.g. random between 0–5
P_init <- runif(omega, min = 0, max = 1)   # e.g. random between 0–2

out <- step_host_parasitoid(H_init, P_init,
                            dist_mat = D,
                            dxbar_H  = 5,  dxbar_P = 10,
                            delta_x  = 1,
                            habitat_suitability        = habitat_suitability,
                            lambda   = 2,
                            a_H      = 1,
                            xi       = 0.5,
                            a_P      = 1,
                            phi      = 0.3)

H_next <- out$H
P_next <- out$P
