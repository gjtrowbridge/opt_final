library(MASS)

part1c = function(x, A, b) {
  u = rep(0, n)
  for (i in 1:20) {
    direction = get_descent_direction(u, x, A, b)
    step_size = 0.001
    u = u - direction * step_size
  }
  return(u)
}

get_loss = function(u, x, A, b) {
  part1 = sum(abs(A %*% x - b))
  part2 = 0.5 * t(u - x) %*% (u - x)
  return(part1 + part2)
}

get_descent_direction = function(u, x, A, b) {
  n = length(x)
  m = dim(A)[1]
  
  
  regular_deriv = u - x
  
  pos_check = A %*% u - b
  pos_check = pos_check >= 0
  new_deriv = rep(0, n)
  
  for (j in 1:n) {
    for (i in 1:m) {
      next_term = A[i, j]
      if (!pos_check[j]) {
        next_term = next_term * -1
      }
      new_deriv[j] = new_deriv[j] + next_term
    }
  }
  
  return(new_deriv + regular_deriv)
}
