library(MASS)

part1c = function(x, A, b) {
  u = rep(0, n)
  for (i in 1:20) {
    direction = get_descent_direction(u, x, A, b)
    step_size = armijos_step_size(u, x, A, b, 0.5, direction, tau = 0.5, c = 0.5)
    u = u - direction * 0.001
    print(get_loss(u,x,A,b))
    print(u)
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

armijos_step_size = function(u, x, A, b, max_step_size, gradient, tau = 0.5, c = 0.5) {
  m = gradient %*% -gradient
  t = -c * m
  step_size = max_step_size
  while (get_loss(u, x, A, b) - get_loss(u - step_size * gradient, x, A, b) < step_size * t) {
    step_size = step_size * tau
  }
  return(step_size)
}

m = 5
n = 3
x = rep(1, n)
A = mvrnorm(m, rep(0, n), diag(rep(1, n)))
A = matrix(rep(3, n*m), nrow=m)
b = rep(2, m)

part1c(x, A, b)