# spline box R-version 1.7

# bsplines(x, knots, order, derivative)
# bspline(x, knots, order, derivative, j)
# bspline_zero(x, knots, order, j)
# bspline_derivative(x, knots, order, derivative, j)
# bspline_jump(knots, order, transpose)
# nsplines(x, knots, derivative)
# nspline_ab(knots)
# nspline_jump(knots)
# knots_quantile(x, dimension, type, order)
# add_boundary_knots(x, interior_knots, order, tiny)
# pbsplines(x, knots, order, derivative)

# matrix of a bspline basis
bsplines = function(x, knots, order = 2, derivative = 0)
{
   dimension = length(knots) - order
   bspline_matrix = matrix(0, length(x), dimension)
   for (j in 1 : dimension)
      bspline_matrix[, j] = bspline(x, knots, order, derivative, j)
   return(bspline_matrix)
}

# a bspline function
bspline = function(x, knots, order = 1, derivative = 0, j = 1)
{
   bspline_vector = rep(0, length(x))
   support = knots[j] <= x & x < knots[j + order]
   if (sum(support) > 0)
   {
      if (derivative == 0)
         bspline_vector[support] = bspline_zero(x[support], knots, order, j)
      else
         bspline_vector[support] = bspline_derivative(x[support], knots, order,
                                                      derivative, j)
   }
   return(bspline_vector)
}

# a bspline function when derivative equals to zero
bspline_zero = function(x, knots, order, j)
{
   if (order > 1)
   {
      if (knots[j + order - 1] > knots[j])
         a = (x - knots[j]) / (knots[j + order - 1] - knots[j])
      else
         a = 0
      if (knots[j + order] > knots[j + 1])
         b = (knots[j + order] - x) / (knots[j + order] - knots[j + 1])
      else
         b = 0
      return(a * bspline_zero(x, knots, order - 1, j) +
                b * bspline_zero(x, knots, order - 1, j + 1))
   }
   else
   {
      bspline = rep(0, length(x))
      bspline[knots[j] <= x & x < knots[j + 1]] = 1
      return(bspline)
   }
}

# derivative of a bspline function
bspline_derivative = function(x, knots, order, derivative, j)
{
   if (derivative > 0)
   {
      if (knots[j + order - 1] > knots[j])
         a = (order - 1) / (knots[j + order - 1] - knots[j])
      else
         a = 0
      if (knots[j + order] > knots[j + 1])
         b = (order - 1) / (knots[j + order] - knots[j + 1])
      else
         b = 0
      return(a * bspline_derivative(x, knots, order - 1, derivative - 1, j) -
                b * bspline_derivative(x, knots, order - 1, derivative - 1, j + 1))
   }
   else
      return(bspline_zero(x, knots, order, j))
}

# compute transposed jump matrix for bspline
bspline_jump = function(knots, order)
{
   dimension = length(knots) - order
   number_interior_knots = dimension - order
   midpoint_between_knots = (knots[order : (length(knots) - order)]
                             + knots[(order + 1) : (length(knots) - order + 1)]) / 2
   # transpose of theoritical jump matrix
   jump_matrix = matrix(nrow = dimension, ncol = number_interior_knots)
   for (j in 1 : dimension)
   {
      derivatives = bspline(midpoint_between_knots, knots, order, order - 1, j)
      for (l in 1 : number_interior_knots)
         jump_matrix[j, l] = derivatives[l + 1] - derivatives[l]
   }
   # re-transpose the transposed jump matirx
   return(jump_matrix)
}

# collection(matrix) of n-spline basis(each column)
nsplines = function(x, knots, derivative = 0)
{
   order = 4
   d = length(knots) - order # = dimension of bsplines
   bspline_matrix = bsplines(x, knots, order, derivative)
   ab = nspline_ab(knots)
   a = ab$a
   b = ab$b
   bspline_matrix[, 3] = bspline_matrix[, 3] +
      a[1] * bspline_matrix[, 1] + b[1] * bspline_matrix[, 2]
   bspline_matrix[, 4] = bspline_matrix[, 4] +
      a[2] * bspline_matrix[, 1] + b[2] * bspline_matrix[, 2]
   bspline_matrix[, d - 2] = bspline_matrix[, d - 2] +
      a[3] * bspline_matrix[, d] + b[3] * bspline_matrix[, d - 1]
   bspline_matrix[, d - 3] = bspline_matrix[, d - 3] +
      a[4] * bspline_matrix[, d] + b[4] * bspline_matrix[, d - 1]
   nspline_matrix = bspline_matrix[, 3 : (d - 2)]
   return(nspline_matrix)
}

# linear functions at knots a and b for n-spline
nspline_ab = function(knots)
{
   order = 4
   dimension = length(knots) - order
   a = b = rep(0, 4)
   x_left = knots[4] + 0.5 * (knots[5] - knots[4])
   second_deriv_left = third_deriv_left = rep(0, 4)
   for (j in 1 : 4)
   {
      second_deriv_left[j] = bspline_derivative(x_left, knots, order, 2, j)
      third_deriv_left[j] = bspline_derivative(x_left, knots, order, 3, j)
   }
   denominator_left = second_deriv_left[2] * third_deriv_left[1] -
      third_deriv_left[2] * second_deriv_left[1]
   b[1] = -(second_deriv_left[3] * third_deriv_left[1] -
               third_deriv_left[3] * second_deriv_left[1]) / denominator_left
   b[2] = -(second_deriv_left[4] * third_deriv_left[1] -
               third_deriv_left[4] * second_deriv_left[1]) / denominator_left
   a[1] = -(second_deriv_left[3] + b[1] * second_deriv_left[2]) / second_deriv_left[1]
   a[2] = -(second_deriv_left[4] + b[2] * second_deriv_left[2]) / second_deriv_left[1]
   # right coefficients of natural splines
   x_right = knots[dimension] + 0.5 * (knots[dimension + 1] - knots[dimension])
   second_deriv_right = third_deriv_right = rep(0, 4)
   for (j in 1 : 4)
   {
      k = dimension - j + 1
      second_deriv_right[j] = bspline_derivative(x_right, knots, order, 2, k)
      third_deriv_right[j] = bspline_derivative(x_right, knots, order, 3, k)
   }
   denominator_right = second_deriv_right[2] * third_deriv_right[1] -
      third_deriv_right[2] * second_deriv_right[1]
   b[3] = -(second_deriv_right[3] * third_deriv_right[1] -
               third_deriv_right[3] * second_deriv_right[1]) / denominator_right
   b[4] = -(second_deriv_right[4] * third_deriv_right[1] -
               third_deriv_right[4] * second_deriv_right[1]) / denominator_right
   a[3] = -(second_deriv_right[3] + b[3] * second_deriv_right[2]) / second_deriv_right[1]
   a[4] = -(second_deriv_right[4] + b[4] * second_deriv_right[2]) / second_deriv_right[1]
   return(list(a = a, b = b))
}

# compute jump matrix for natural spline
nspline_jump = function(knots)
{
   order = 4
   dimension_bs = length(knots) - order   # J (= bs_num_interior - order)
   dimension_ns = dimension_bs - order   # J - 4 (= ns_num_interior = bs_num_interior)
   ns_num_interior_knots = dimension_ns
   jump_bs = bspline_jump(knots, order)   # bs jump matrix, (dim_bs - order) X dim_bs
   jump_ns = matrix(nrow = ns_num_interior_knots, ncol = dimension_ns)   # ns jump matrix, (dim_bs - order) X (dim_bs - order)
   ab = nspline_ab(knots)
   a = ab$a
   b = ab$b
   # at left boundary knots
   jump_ns_temp = rep(0, dimension_ns)
   jump_ns_temp = jump_bs[3, ] + a[1] * jump_bs[1, ] + b[1] * jump_bs[2, ]
   jump_bs[3, ] = jump_ns_temp
   jump_ns_temp = jump_bs[4, ] + a[2] * jump_bs[1, ] + b[2] * jump_bs[2, ]
   jump_bs[4, ] = jump_ns_temp
   # at right boundary constraints
   jump_ns_temp = jump_bs[dimension_bs - 2, ] + 
      a[3] * jump_bs[dimension_bs, ] + b[3] * jump_bs[dimension_bs - 1, ]
   jump_bs[dimension_bs - 2, ] = jump_ns_temp
   jump_ns_temp = jump_bs[dimension_bs - 3, ] + 
      a[4] * jump_bs[dimension_bs, ] + b[4] * jump_bs[dimension_bs - 1, ]
   jump_bs[dimension_bs - 3, ] = jump_ns_temp
   for (j in 1:dimension_ns)
      jump_ns[j, ] = jump_bs[j + 2, ]
   return(jump_ns)
}

# setting knots by quantile of data
knots_quantile = function(x, dimension, order = 2, type = "bs")
{
   if (type == "bs")
   {
      dimension = max(dimension, order)
      number_interior_knots = dimension - order
      if (number_interior_knots > 0)
         probs = (1 : number_interior_knots) / (number_interior_knots + 1)
      else
         probs = NULL
   }
   else
   {
      dimension = max(dimension, 2)
      probs = seq(0, 1, length = dimension)
      order = 4
   }
   interior_knots = quantile(x, probs, type = 1)
   knots = add_boundary_knots(x, interior_knots, order)
   return(knots)
}

# adding boundary knots to interior knots
add_boundary_knots = function(x, interior_knots, order = 1, tiny = 1e-5)
{
   knots = c(rep(min(x) - tiny, order), interior_knots, rep(max(x) + tiny, order))
   return(knots)
}

# periodic bsplines
pbsplines = function(x, knots, order = 1, 
                     lower = 0, upper = 1, derivative = 0)
{
   dimension = length(knots)
   if (dimension < order)
      stop("dimension < order")
   x = x - lower
   knots = knots - lower
   period = upper - lower
   knots_extension = c(knots, period + knots[1 : order])
   pbspline_matrix = matrix(0, length(x), dimension)
   for (j in 1 : (dimension - order))
   {
      knots_j = knots_extension[j : (j + order)]
      pbspline_matrix[, j] = bspline(x, knots_j, order, derivative, 1) 
   }
   x_extention = x
   left = x <= knots[order]
   x_extention[left] = x[left] + period
   for (j in (dimension - order + 1) : dimension)
   {
      knots_j = knots_extension[j : (j + order)]
      pbspline_matrix[, j] = bspline(x_extention, knots_j, order, derivative, 1) 
   }
   return(pbspline_matrix)
}

# work notes (maintainer: JKS)
#  version 1.9 (2020.05.01)
#     'nspline_jump' is revised so that 
#      row: dimension, col: the number of constraints
#  version 1.8 (2020.03.15)
#     removed 'transpose' option on 'bspline_jump()'
#     revised all arguments 'dim' -> 'dimension'
#     revised 'nspline_jump' to be equal to Rcpp version splineBox2.1
#  version 1.7 (2020.03.09)
#     added pbsplines
#     bspline => bspline_zero
#     comments are revised
#     order of functions are changed
#  version 1.6 (2019.08.22)
#     added 'nspline_jump()' which is converted from 'ns_jump()' in KYB's 'psdce_natural.cpp'
#  version 1.5 (2019.08.16)
#     removed 'bspline_derivative_point()' and 'bspline_point()',
#        because they are same to 'bspline_derivative()' and 'bspline_zero()'
#  version 1.4 (2019.07.21)
#     added 'transpose' option on 'bspline_jump()'
#  version 1.3 (2019.07.20)
#     added 'bspline_jump()', the converted version of 'jump_bbspline()' in KYB's 'psdce_natural.cpp'
#  version 1.2 (2019.07.20)
#     added 'bspline_derivative_point()' and 'bspline_point()' made by JYK
#  version 1.1 (2019.03.26)
#     added comments for functions
#  version 1.0 (2019.03.25)
#     received original code from JYK