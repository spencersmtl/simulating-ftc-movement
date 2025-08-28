library(tidyverse)

# dij is distance between i and j
# dxbar is average dispersal distance of organism x
# distances is all dij for one point. AKA ith row of distance matrix
# delta_x is a scaling factor, it represents the area of one grid cell
# D is the distance matrix
# H / P is the host / parasitoid density after growth step
# H_star / P_star is the host / parasitoid density after dispersal
# a is parasitoid search rate
# I is a habitat suitability vector

# 10x10 tester
load(file = "data/10x10_landscape_info.RData")
load(file = "data/10x10_dist_matrix.RData")
K_H <- gp_movement_matrix
K_P <- gp_movement_matrix/2
omega <- nrow(K)
landscape_dataframe <- gp_landscape
habitat_suitability <- 1-round(landscape_dataframe$value)
lambda <- 2.61 # Host intrinsic growth rate
xi <- 8.06 # Parasitoid intrinsic growth rate
phi <- 0.35 # Parasitoid phenology. phi = 0 means P escape H density dependent mortality


# Dispersal kernel ####
# 1. movement kernel
m_kernel <- function(dij, dxbar) { 
  2 / (pi * dxbar^2) * exp(-2 * dij / dxbar)
}

# 2. Cut‑off distance d_max(dxbar)
d_max <- function(dxbar) {
  -0.5 * dxbar * log(pi * dxbar^2 * 10e-10 / 2)
}

# 3. Single‐entry of normalized kernel k_ij of species x
kernel_element <- function(dij, dxbar, delta_x) {
  dm   <- d_max(dxbar)
  if (dij > dm) return(0)
  numerator  <- m_kernel(dij, dxbar) * delta_x^2
  denominator <- sum(m_kernel(dij[dij <= dm], dxbar) * delta_x^2)
  numerator / denominator
}

# 4. Full kernel matrix K (Ω×Ω) for species x
compute_kernel_matrix <- function(distances, dxbar, delta_x) {
  omega <- nrow(distances)
  K <- matrix(0, nrow=omega, ncol=omega)
  for (i in seq_len(omega)) {
    d_row <- distances[i, ]
    for (j in seq_len(omega)) {
      K[i,j] <- kernel_element(dij      = d_row[j],
                               dxbar     = dxbar,
                               delta_x   = delta_x)
    }
  }
  K
}


# Compute dispersal ####
compute_H_star <- function(Hj, K_H) { # Density of hosts after dispersing
  as.vector(K_H * Hj) # Summed (dispersal probabilities from i to all j's) * (densities at j)
}

compute_P_star <- function(Pj, K_P) { # Density of parasitoids after dispersing
  as.vector(K_P * Pj)
}

# host growth ####
density_regulation <- function(H_star, lambda) { # Density‐dependence: exp(−log(λ)·H*)
  exp(2 * -log(lambda) * H_star)
}

host_parasitism_survival <- function(P_star, a = 1) { # Survival after parasitism: exp(−a·P*)
  exp(-a * P_star)
}

update_H <- function(H_star, P_star, habitat_suitability, lambda, a = 1) { # Full host growth. h(H*,P*)
  habitat_suitability * lambda * H_star * 
    density_regulation(H_star, lambda) * 
    host_parasitism_survival(P_star, a)
}

# Parasitoid growth ####
parasitoid_attack <- function(P_star, a = 1) { # T2 functional response: 1 − exp(−a·P*)
  1 - exp(-a * P_star)
}

parasitoid_regulation <- function(H_star, lambda, phi = 0) { # Host‐density feedback of parasitoids: exp(−ϕ·log(λ)·H*)
  exp(-phi * log(lambda) * H_star)
}

compute_xi <- function(a, b, K) { # Parasitoid recruitment rate: xi = a * b * K
  xi <- a * b * K
  return(xi)
}

update_P <- function( # Full parasitoid growth
    H_star, P_star, habitat_suitability, xi, lambda, a = 1, phi = 0) {
  habitat_suitability * xi * H_star *
    parasitoid_attack(P_star, a) *
    parasitoid_regulation(H_star, lambda, phi)
}

# Step simulator ####

step_HP_kernel <- function(H, P, # densities
                           K_H, K_P, #  dispersal kernels
                           habitat_suitability, # habitat suitability
                           lambda,       # host growth
                           a_H  = 1,     # host parasitism coefficient
                           xi,           # parasitoid recruitment
                           a_P  = 1,     # parasitoid attack
                           phi  = 0) { # feedback strength
  # Computes a single time step, t for all four equations in model

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
