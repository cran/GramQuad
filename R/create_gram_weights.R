# File GramQuad/R/create_gram_weights.R

create_gram_weights <- function(m){
  # param  m: the degree  of  the  largest possible Gram  polynomialfor m+1 points.
  # return : gram  weights for gram  quadrature of m+1 points .
  max_d <- as.integer(sqrt(m))
  xs <- seq(-1, 1, length.out = m + 1)
  
  alphas <- c(1, sapply(0:(max_d), alpha_func, m))
  
  # get GQ  points  and  weights
  lg <- pracma::gaussLegendre(max_d / 2 + 1, -1, 1)
  xs_lg <- lg$x
  w_lg <- lg$w
  
  # denote  A_p as  the  previous row of A being  worked  on. A_n is the newest  row.
  uw <- list(
             w = double(m + 1),
             A_p = double(m + 1),
             A_n = rep((m + 1)^(-0.5), m + 1),
	     q_p = double(as.integer(max_d / 2 + 1)),
	     q_n = w_lg * rep((m + 1)^(-0.5), as.integer(max_d / 2 + 1))
  )
  
  
  # update  w for n = 0, 1, 2, 3, ..., max_d
  for(j in 1:(max_d + 1)){
    uw <- update_w(uw$w, uw$A_p, uw$A_n, uw$q_n, uw$q_p, alphas, xs, xs_lg, j)
  }
  return(uw$w)
}

