get_M_LEN <- function(col_num) {
  return((col_num-6)/5)
}
get_F_END <- function(M_LEN) {
  return(2+(M_LEN+1)*4-1)
}
get_E_START <- function(M_LEN) {
  return(2+(M_LEN+1)*4)
}
get_E_END <- function(M_LEN) {
  return(2+(M_LEN+1)*4-1+M_LEN)
}

get_fisher <- function(d) {
  M_LEN <- get_M_LEN(ncol(d))
  F_END <- get_F_END(M_LEN)
  E_START <- get_E_START(M_LEN)
  E_END <- get_E_END(M_LEN)
  fisher_data <- d[, 2:F_END]
  return(fisher_data)
}

get_pca_fisher <- function(d, center=T, scale=F) {
  M_LEN <- get_M_LEN(ncol(d))
  F_END <- get_F_END(M_LEN)
  E_START <- get_E_START(M_LEN)
  E_END <- get_E_END(M_LEN)
  fisher_data <- d[, 2:F_END]
  pca <- prcomp(fisher_data, center=center, scale=scale)
  return(pca)
}
