# File GramQuad/R/utils.R
# internal (utility) functions



alpha_func <- compiler::cmpfun(function(n, m){
  a1 <- m / (n + 1)
  a2 <- (4 * (n + 1)^2) - 1
  a3 <- (m + 1)^2 - (n + 1)^2
  a4 <- a2 / a3
  
  res <- a1 * sqrt(a4)
  return(res)
})

update_w <- compiler::cmpfun(function(w, A_p, A_n, q_n, q_p, alphas, xs, xs_lg, j){
  w <- w + sum(q_n) * A_n
  
  # update  for j = 1, 2, 3, ...,  max_d+1
  A_temp <- (alphas[j+1] * (xs * A_n)) - ((alphas[j+1] / alphas[j]) * A_p)
  A_p <- A_n
  A_n <- A_temp
  
  q_temp <- q_n
  q_n <- (alphas[j+1] * xs_lg * q_n) - (alphas[j+1] / alphas[j] * q_p)
  
  q_p <- q_temp
  
  return(list(w = w, A_p = A_p, A_n = A_n, q_n = q_n, q_p = q_p))
})
