pre_host_size <- list()

# Host growth equation
host_growth <- function(pre_host_size, 
                        suitability, 
                        intrinsic_host_growth,
                        pre_host_size,
                        pre_parasitoid_growth,
                        attack_rate)
{
  host_size <- 
    suitability * 
    intrinsic_host_growth * 
    pre_host_size * 
    exp(-log(intrinsic_host_growth) * pre_host_size) *
    exp(-pre_parasitoid_growth * attack_rate)
}

# Parasitoid growth equation
parasitoid_growth <- function(){
  
}

# Dispersal
dispersal <- function(){
  
}