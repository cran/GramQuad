\name{create_gram_weights}

\alias{create_gram_weights}

\title{Compute Gram weights}

\description{Computes weights for Gram quadrature of\code{m+1} points.}

\usage{create_gram_weights(m)}

\arguments{
\item{m}{a positive integer value, the number of points minus one for
  which weights will be computed. See 'Details'.}
} 

\value{
  A double-precision vector of the specified length plus one, whose
  elements are the weights for the Gram quadratures of the \code{m+1}
  points in the interval \eqn{[-1,1]}.
}


\details{
  The numerical integration of an analytical function \eqn{f(x)} using a
finite set of equidistant points can be performed by quadrature formulas
like the Newton-Cotes. Unlike Gaussian quadrature formulas however,
higher-order Newton-Cotes formulas are not stable, limiting the usable
order of such formulas. Existing work showed that by the use of
orthogonal polynomials, stable high-order quadrature formulas with
equidistant points can be developed. This algorithm improves upon such
work by making use of (orthogonal) Gram polynomials and deriving an
iterative algorithm, together allowing us to reduce the space-complexity
of the original algorithm significantly.
}


\references{
  Muhammad, Irfan (2021) \emph{Gram quadrature: Numerical integration with Gram
    polynomials}. \href{https://arxiv.org/pdf/2106.14875.pdf}{arXiv:2106.14875 [math.NA]}
}

\seealso{\code{\link[pracma]{gaussLegendre}}}

\examples{
m <- 100
xs <- seq(-1, 1, length.out = m + 1)
gram_weights <- create_gram_weights(m)

# the sum of  stable  weights is  equal  to 2.
cat("Sum of Gram weights:", sum(gram_weights), "\n")

# test integration, integrate f below between [-1,1]
f = function(x){ 9 * x ^ 2 + 45 * 13 * x ^ 3 + 16 * x ^ 4}
gram_quad <- sum(gram_weights * f(xs))
cat("Approx. integration:", gram_quad, "\n")
}
\author{
  Iago Giné-Vázquez, \email{iago.gin-vaz@protonmail.com}
}
\keyword{Gram quadrature}
