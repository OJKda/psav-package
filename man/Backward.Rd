\name{Backward}
\alias{Backward}
\title{Knot elimination via Backward}
\usage{
Backward(x, z, y, inter_knots, left_knots, right_knots, add_knot,
         existing_knot, order, criter, selection_number)
}
\arguments{
   \item{x}{\eqn{n \times 1} predictor vector. \eqn{n} is a sample size.}
   \item{z}{\eqn{n \times K} auxiliary binary variable matrix.
            \eqn{K} is the number of levels for auxiliary variable.}
   \item{y}{\eqn{n \times 1} response vector.}
   \item{inter_knots}{interior knots}
   \item{left_knots}{left boundary knots}
   \item{right_knots}{right boundary knots}
   \item{add_knot}{knots selected by forward}
   \item{existing_knot}{existing knots excluding newly selected knot
                        (\code{add_knot})}
   \item{order}{order of B-splines. 1 fits constant, 2 fits linear,
                3 fits quadratic, 4 fits cubic splines.}
   \item{criter}{\code{criter} has a numeric of 8 or 9.
                 8 indicates \code{"BIC"} and 9 indicates \code{"AIC"}}
   \item{selection_number}{indicates the numberth of knots selected.}
}
\description{
\code{Backward} is used as a function to remove knots in \code{stepwise}.
Among the selected knots, the oldest knot is removed and a new knot is
selected again.
}
