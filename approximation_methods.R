#################################################################
# Approximate algorithms and methods for different depths in R^d
#################################################################

#' Universal Random Search for Depth Approximation
#'
#' @param data Data matrix (n x d)
#' @param point Point to evaluate (vector of length d)
#' @param depth.func Depth function with signature: function(data, point, direction, ...)
#' @param n_it Number of iterations (default: 1000)
#' @param ... Additional parameters passed to depth.func
#' @return Minimum depth value found

random.search.universal <- function(data, point, depth.func, n_it = 1000, ...) {
  # Input validation
  stopifnot(
    is.matrix(data),
    ncol(data) == length(point),
    is.function(depth.func)
  )
  
  # Direction generator function
  rnd.sphere <- function(d) {
    v <- rnorm(d)
    v / sqrt(sum(v^2))
  }
  
  # Wrapper function to handle additional parameters
  depth_func_wrapper <- function(data, point, direction, ...) {
    # Check if the depth function accepts the parameter 'p'
    if ("p" %in% as.character(formals(depth.func))) {
      do.call(depth.func, list(data = data, point = point, direction = direction, p = 0.1, ...))
    } else {
      do.call(depth.func, list(data = data, point = point, direction = direction, ...))
    }
  }
  
  # Initialize minimum depth
  min_depth <- 1  # Initialize with maximum possible depth
  
  for (i in 1:n_it) {
    # Generate random direction
    direction <- rnd.sphere(ncol(data))
    
    # Calculate depth with error handling
    current_depth <- tryCatch({
      depth_func_wrapper(data = data,
                         point = point,
                         direction = direction,
                         ...)
    }, error = function(e) {
      warning("Error in depth calculation: ", e$message)
      NA
    })
    
    # Update minimum depth (ignoring NAs)
    if (!is.na(current_depth) && current_depth < min_depth) {
      min_depth <- current_depth
    }
  }
  
  return(min_depth)
}


#' Universal Refined Random Search for Depth Approximation in R^d
#'
#' @param data Data matrix (n x d)
#' @param point Point to evaluate (vector of length d)
#' @param depth.func Depth function with signature: function(data, point, direction, ...)
#' @param n_it Number of iterations per refinement step (default: 100)
#' @param n_ref Number of refinement steps (default: 10)
#' @param alpha Shrinking factor for the neighborhood (default: 0.5)
#' @param ... Additional parameters passed to depth.func
#' @return Minimum depth value found

refined.random.search.universal <- function(data, point, depth.func, n_it = 100, n_ref = 10, alpha = 0.5, ...) {
  # Input validation
  stopifnot(
    is.matrix(data),
    ncol(data) == length(point),
    is.function(depth.func)
  )
  
  # Direction generator function
  rnd.sphere <- function(d) {
    v <- rnorm(d)
    v / sqrt(sum(v^2))
  }
  
  # Spherical cap generator function
  rnd.spherical.cap <- function(p, epsilon) {
    d <- length(p)
    x <- cos(epsilon * runif(1))
    v <- c(x, sqrt(1 - x^2) * rnd.sphere(d - 1))
    householder(v, p)
  }
  
  # Householder transformation
  householder <- function(x, p) {
    if (p[1] == 1) {
      return(x)
    } else {
      lambda <- (sum(p * x) - x[1]) / (1 - p[1])
      x[1] <- x[1] + lambda
      return(x - lambda * p)
    }
  }
  
  # Wrapper function to handle additional parameters
  depth_func_wrapper <- function(data, point, direction, ...) {
    # Check if the depth function accepts the parameter 'p'
    if ("p" %in% as.character(formals(depth.func))) {
      do.call(depth.func, list(data = data, point = point, direction = direction, p = 0.1, ...))
    } else {
      do.call(depth.func, list(data = data, point = point, direction = direction, ...))
    }
  }
  
  # Initialize parameters
  p_min <- rep(0, length(point))
  p_min[1] <- 1
  min_depth <- depth_func_wrapper(data = data, point = point, direction = p_min, ...)
  epsilon <- pi / 2
  
  # Refined search
  for (i in 1:n_ref) {
    for (j in 1:n_it) {
      # Generate random direction within spherical cap
      direction <- rnd.spherical.cap(p_min, epsilon)
      
      # Calculate depth with error handling
      current_depth <- tryCatch({
        depth_func_wrapper(data = data,
                           point = point,
                           direction = direction,
                           ...)
      }, error = function(e) {
        warning("Error in depth calculation: ", e$message)
        NA
      })
      
      # Update minimum depth and direction (ignoring NAs)
      if (!is.na(current_depth) && current_depth < min_depth) {
        min_depth <- current_depth
        p_min <- direction
      }
    }
    # Shrink the neighborhood
    epsilon <- epsilon * alpha
  }
  
  return(min_depth)
}



#' Universal Spherical Nelder-Mead for Depth Approximation
#'
#' @param data Data matrix (n x d)
#' @param point Point to evaluate (vector of length d)
#' @param depth.func Depth function with signature: function(data, point, direction, ...)
#' @param n_it Maximum number of iterations (default: 1000)
#' @param alpha Reflection coefficient (default: 1)
#' @param gamma Expansion coefficient (default: 2)
#' @param rho Contraction coefficient (default: 0.5)
#' @param sigma Shrink coefficient (default: 0.5)
#' @param start.Mn Logical: Start from mean direction (default: TRUE)
#' @param ... Additional parameters passed to depth.func
#' @return Minimum depth value found

spherical.nelder.mead.universal <- function(data, point, depth.func, n_it = 1000,
                                            alpha = 1, gamma = 2, rho = 0.5, sigma = 0.5,
                                            start.Mn = TRUE, ...) {
  # Input validation
  stopifnot(
    is.matrix(data),
    ncol(data) == length(point),
    is.function(depth.func)
  )
  
  # Great circle function
  great.circle <- function(xx, yy, t, bound = TRUE) {
    sp <- sum(xx * yy)
    if (sp >= 0) {
      suma <- sum((xx - yy)^2)
      alpha <- 2 * asin(sqrt(suma) / 2)
      sina <- sqrt(suma * (1 + sp) / 2)
    } else {
      suma <- sum((xx + yy)^2)
      alpha <- pi - 2 * asin(sqrt(suma) / 2)
      sina <- sqrt(suma * (1 - sp) / 2)
    }
    gx <- (1 - t) * alpha
    gy <- t * alpha
    if (bound & (abs(gy) > pi / 2)) {
      if (gy > 0) {
        gy <- pi / 2
      } else {
        gy <- -pi / 2
      }
      gx <- alpha - gy
    }
    cx <- sin(gx) / sina
    cy <- sin(gy) / sina
    zz <- cx * xx + cy * yy
    return(zz)
  }
  
  # Householder transformation
  householder <- function(xx, pp) {
    if (pp[1] == 1) {
      return(xx)
    } else {
      lambda <- (sum(pp * xx) - xx[1]) / (1 - pp[1])
      xx[1] <- xx[1] + lambda
      return(xx - lambda * pp)
    }
  }
  
  # Random sphere function
  rnd.sphere <- function(d) {
    v <- rnorm(d)
    v / sqrt(sum(v^2))
  }
  
  # Spherical cap generator function
  rnd.spherical.cap <- function(p, epsilon) {
    d <- length(p)
    x <- cos(epsilon * runif(1))
    v <- c(x, sqrt(1 - x^2) * rnd.sphere(d - 1))
    householder(v, p)
  }
  
  # Wrapper function to handle additional parameters
  depth_func_wrapper <- function(data, point, direction, ...) {
    # Check if the depth function accepts the parameter 'p'
    if ("p" %in% as.character(formals(depth.func))) {
      do.call(depth.func, list(data = data, point = point, direction = direction, p = 0.1))
    } else {
      do.call(depth.func, list(data = data, point = point, direction = direction))
    }
  }
  
  # Initialize parameters
  if (start.Mn) {
    uu <- point - colMeans(data)
    uu <- uu / sqrt(sum(uu^2))
  } else {
    uu <- rnd.sphere(length(point))
  }
  
  epsilon <- (pi / 2) / 2
  depth <- rep(0, length(point) + 1)
  pp <- matrix(nrow = length(point) + 1, ncol = length(point))
  
  # Initial simplex
  for (i in 1:(length(point) + 1)) {
    pp[i, ] <- rnd.spherical.cap(uu, epsilon)
    depth[i] <- depth_func_wrapper(data = data, point = point, direction = pp[i, ], ...)
  }
  
  # Order simplex
  pp <- pp[order(depth), ]
  depth <- sort(depth)
  
  # Nelder-Mead optimization
  it <- 0
  while (it < n_it & sum((pp[1, ] - pp[length(point) + 1, ])^2) > 1e-20) {
    it <- it + 1
    if (length(point) > 1) {
      xx0 <- colMeans(pp[-(length(point) + 1), ])
    } else {
      xx0 <- pp[1, ]
    }
    xx0 <- xx0 / sqrt(sum(xx0^2))
    xxr <- great.circle(xx = xx0, yy = pp[length(point) + 1, ], t = -alpha, bound = TRUE)
    depthr <- depth_func_wrapper(data = data, point = point, direction = xxr, ...)
    
    if (depth[1] <= depthr & depthr < depth[length(point)]) {
      pp[length(point) + 1, ] <- xxr
      depth[length(point) + 1] <- depthr
    } else if (depthr < depth[1]) {
      xxe <- great.circle(xx = xx0, yy = xxr, t = gamma, bound = TRUE)
      depthe <- depth_func_wrapper(data = data, point = point, direction = xxe, ...)
      if (depthe < depthr) {
        pp[length(point) + 1, ] <- xxe
        depth[length(point) + 1] <- depthe
      } else {
        pp[length(point) + 1, ] <- xxr
        depth[length(point) + 1] <- depthr
      }
    } else {
      if (depthr < depth[length(point) + 1]) {
        xxh <- xxr
      } else {
        xxh <- pp[length(point) + 1, ]
      }
      xxc <- great.circle(xx = xx0, yy = xxh, t = rho, bound = TRUE)
      depthc <- depth_func_wrapper(data = data, point = point, direction = xxc, ...)
      if (depthc < depth[length(point) + 1]) {
        pp[length(point) + 1, ] <- xxc
        depth[length(point) + 1] <- depthc
      } else {
        for (i in 2:(length(point) + 1)) {
          pp[i, ] <- great.circle(xx = pp[1, ], yy = pp[i, ], t = sigma, bound = TRUE)
          depth[i] <- depth_func_wrapper(data = data, point = point, direction = pp[i, ], ...)
        }
        pp[-(length(point) + 1), ] <- pp[order(depth[-(length(point) + 1)]), ]
        depth[-(length(point) + 1)] <- sort(depth[-(length(point) + 1)])
      }
    }
    pp <- pp[order(depth), ]
    depth <- sort(depth)
  }
  
  return(depth[1])
}


# Example of usage
# Define a custom depth function
# custom_depth <- function(data, point, direction, p = 0.1) {
#   # Compute the depth using the provided direction
#   # Your depth computation logic here
# }

# Generate sample data
# set.seed(123)
# data <- matrix(rnorm(2000), ncol = 10)  # 200 x 10
# point <- rnorm(10)  # Length 10

# Run the Random Search Universal algorithm
# result <- random.search.universal(
#   data = data,
#   point = point,
#   depth.func = custom_depth,
#   n_it = 1000,
#   p = 0.1
# )
# 
# print(result)