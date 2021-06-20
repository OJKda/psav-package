\name{fit_model}
\alias{fit_model}
\title{Pilable Regression Spline model using B-splines and auxiliary variables}
\usage{
fit_model(x, z, y, knots, order, max_iter = 1000, epsilon = 1e-5)
}
\arguments{
   \item{x}{\eqn{n \times 1} predictor vector. \eqn{n} is a sample size.}
   \item{z}{\eqn{n \times K} auxiliary binary variable matrix.
            \eqn{K} is the number of levels for auxiliary variable.}
   \item{y}{\eqn{n \times 1} response vector.}
   \item{knots}{selected knots among the quantiles of \code{x}}
   \item{order}{order of B-splines. 1 fits constant, 2 fits linear,
                3 fits quadratic, 4 fits cubic splines.}
   \item{max_iter}{maximum iteration of Coodinate Descent Algorithm}
   \item{epsilon}{stopping criteria of Coodinate Descent Algorithm.
                  If the difference in the residual sum of squares is less                        than epsilon, the algorithm is stopped.}
}
\value{
The results are output as a list. The list consists of the following:

   \item{beta}{coefficient vectors for predictor terms}
   \item{gamma}{coefficient matrix for auxiliary variable terms}
   \item{residuals}{the residuals, that is response minus fitted values.}
   \item{fitted_values}{the fitted values}
   \item{iter}{number of iteration}
   \item{inter_knots}{interior knots}
   \item{n}{sample size}
   \item{BIC}{Bayesian Information Criterion}
   \item{AIC}{Akaike Information Criterion}
}

\description{
Fit a pilable regression spline model
}

\details{
We conduct a study on a statistical learning for modeling and analysis of data
with another auxiliary variable in addition to the response and predictor
variables. The estimation is expressed as a linear combination of B-splines.
Coefficients are estimated by minimizing the residual sum of squares via
Coordinate Descent Algorithm. For more details,see Oh and Jhong (0000).
}

\author{
Jae-Kwon Oh and Jae-Hwan Jhong
}

