\name{StepWise}
\alias{StepWise}
\title{Knot Selection via Stepwise Method}
\usage{
StepWise(x, z, y, max_knots = 10, order = 4,
         criterion = "BIC", max_iter = 1000, epsilon = 1e-5)
}
\arguments{
   \item{x}{\eqn{n \times 1} predictor vector. \eqn{n} is a sample size.}
   \item{z}{\eqn{n \times K} auxiliary binary variable matrix.
            \eqn{K} is the number of levels for auxiliary variable.}
   \item{y}{\eqn{n \times 1} response vector.}
   \item{max_knots}{the number of initial knots. should specify as
                    the number of that can be divided into small enough segments.}
   \item{order}{order of B-splines. 1 fits constant, 2 fits linear,
                3 fits quadratic, 4 fits cubic splines.}
   \item{criterion}{scale used as a criterion for knot selection.
                    Criteria can be selected as \code{"BIC"} or \code{"AIC"}.}
   \item{max_iter}{maximum iteration of Coodinate Descent Algorithm}
   \item{epsilon}{stopping criteria of Coodinate Descent Algorithm.
                  If the difference in the residual sum of squares is less                        than epsilon, the algorithm is stopped.}
}
\description{
\code{StepWise} is used for knot selection. At every step, the model fits with the \code{fit_model} function.
}
\details{
Starting from the NULL model, the optimal model is calculated by repeating the selection and removal of knot. Backward is performed when there are more than two knots.
}
\author{
Jae-Kwon Oh and Jae-Hwan Jhong
}
\examples{
# Example for Bone Mineral Density Data

library(psav)
library(loon.data)
data("bone")

# When the auxiliary variable is sex
par(mfrow = c(1,2))
x = bone$age
y = bone$rspnbmd
z = bone$sex
tz = transform(z,
               class1 = ifelse(z == "male", 1, 0),
               class2 = ifelse(z == "female", 1, 0))[,2:3]
x_m = x[z == "male"]
x_f = x[z == "female"]

fit = StepWise(x, tz, y, max_knots = 50, order = 4,
               criterion = "BIC", max_iter = 1000, epsilon = 1e-5)
fit_m = fit$fitted_values[z == "male"]
fit_f = fit$fitted_values[z == "female"]

plot(x, y, bty = "n", type = "n", xlab = "age",
     ylab = "Relative spinal bone mineral density",
     cex.lab = 1.3, cex.axis = 1.2)
points(x, y, col = "grey")
lines(sort(x_m), fit_m[order(x_m)], lty = 1, lwd = 2)
lines(sort(x_f), fit_f[order(x_f)], lty = 3, lwd = 2)
abline(v = fit$knots, col = "grey", lty = 2)
legend("topright", c("male", "female"),
       col = "black",
       lty = c(1, 3), lwd = 2,
       bty = "n", cex = 1.1)


# When the auxiliary variable is ethnic
bone_ethnic <- na.omit(bone)
x <- bone_ethnic$age
y <- bone_ethnic$rspnbmd
z <- bone_ethnic$ethnic

tz = transform(z,
               class1 = ifelse(z == "White", 1, 0),
               class2 = ifelse(z == "Hispanic", 1, 0),
               class3 = ifelse(z == "Asian", 1, 0),
               class4 = ifelse(z == "Black", 1, 0))[, 2:5]
x_w = x[z == "White"]
x_h = x[z == "Hispanic"]
x_a = x[z == "Asian"]
x_b = x[z == "Black"]

fit = StepWise(x, tz, y, max_knots = 50, order = 4,
               criterion = "BIC", max_iter = 1000, epsilon = 1e-5)

fit_w = fit$fitted_values[z == "White"]
fit_h = fit$fitted_values[z == "Hispanic"]
fit_a = fit$fitted_values[z == "Asian"]
fit_b = fit$fitted_values[z == "Black"]

plot(x, y, bty = "n", type = "n", xlab = "age",
     ylab = "Relative spinal bone mineral density",
     cex.lab = 1.3, cex.axis = 1.2)
points(x, y, col = "grey")
lines(sort(x_w), fit_w[order(x_w)], lty = 1, lwd = 2)
lines(sort(x_h), fit_h[order(x_h)], lty = 3, lwd = 2)
lines(sort(x_a), fit_a[order(x_a)], lty = 4, lwd = 2)
lines(sort(x_b), fit_b[order(x_b)], lty = 5, lwd = 2)
abline(v = fit$knots, col = "grey", lty = 2)
legend("topright", c("White", "Hispanic", "Asian", "Black"),
       col = "black",
       bty = "n", cex = 1.1, lwd = 2,
       lty = c(1,3,4,5))
}

