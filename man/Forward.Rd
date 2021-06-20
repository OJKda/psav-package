\name{Forward}
\alias{Forward}
\title{Knot selection via Forward}
\usage{
Forward(x, z, y, inter_knots, left_knots, right_knots,
        pre_fit, order, criter, selection_number)
}
\arguments{
   \item{x}{\eqn{n \times 1} predictor vector. \eqn{n} is a sample size.}
   \item{z}{\eqn{n \times K} auxiliary binary variable matrix.
            \eqn{K} is the number of levels for auxiliary variable.}
   \item{y}{\eqn{n \times 1} response vector.}
   \item{inter_knots}{interior knots}
   \item{left_knots}{left boundary knots}
   \item{right_knots}{right boundary knots}
   \item{pre_fit}{model before choosing a new knot (\code{fit_model} results)}
   \item{order}{order of B-splines. 1 fits constant, 2 fits linear,
                3 fits quadratic, 4 fits cubic splines.}
   \item{criter}{\code{criter} has a numeric of 8 or 9.
                 8 indicates \code{"BIC"} and 9 indicates \code{"AIC"}}
   \item{selection_number}{indicates the numberth of knots selected.}
}
\description{
Choose a knot from the knots other than the knot in \code{pre_fit}
}

