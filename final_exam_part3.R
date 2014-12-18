library(MASS)
library(psych)
library(svd)

get_Gt_x = function(X, step_size, gradient, lambda) {
  return(1/step_size * (X - get_proximal(X-step_size*gradient, step_size, lambda)))
}

get_step_size = function(y, X, A, J, gradient, step_size, lambda, beta=0.5) {
  return(0.001)
}

# get_step_size = function(y, X, A, J, gradient, step_size, lambda, beta=0.5) {
#   Gt_x = get_Gt_x(X, step_size, gradient, lambda)
#   while (get_convex_loss(y, X - step_size * Gt_x, A, J) > 
#          get_convex_loss(y, X, A, J) -
#          step_size * t(gradient) %*% Gt_x +
#          step_size / 2 * norm(Gt_x, "F")) {
#     step_size = step_size * beta
#   }
#   return(step_size)
# }

get_convex_loss = function(y, X, A, J) {
  result = 0
  for (j in 1:J) {
    extra = (y[j] - tr(t(X) %*% A[[j]]))^2
    result = result + extra
  }
  return(result)
}

get_proximal_loss = function(y, X, A, J, lambda) {
  
}

get_proximal = function(X, step_size, lambda) {
  #SVD, then threshold the sigmas using the step size
  a = svd(X)
  d = a[['d']]
  U = a[['u']]
  V = a[['v']]
  
  for (i in 1:length(d)) {
    if (d[i] >= step_size) {
      d[i] = d[i] - step_size
    } else if (d[i] <= -step_size) {
      d[i] = d[i] + step_size
    } else {
      d[i] = 0
    }
  }
  
  D = diag(d)
  return(U %*% D %*% t(V))
}

get_gradient = function(y, X, A, J, m, n) {
  #-2 * fxn without sqr * a_ik at ik, j
  result = matrix(rep(0, m*n),nrow=m)
  for (j in 1:J) {
    result = result + -2 * (y[j] - tr(t(X) %*% A[[j]])) * A[[j]]
  }
  return(result)
}

#Generate data
J = 10
A = list()
n = 5
m = 4

true_X = matrix(rep(5, n * m), nrow=m)
true_X = mvrnorm(m, rep(0, n), diag(rep(1, n)))

for (j in 1:J) {
  A[[j]] = mvrnorm(m, rep(0,n), diag(rep(1,n)))
}
errors = rnorm(J, 0, 0.1)
#errors = rep(0, J)
y = rep(0, J)

for (j in 1:J) {
  y[j] = tr(t(true_X) %*% A[[j]]) + errors[j]
}

#Starting estimate
X_est = matrix(rep(0, m*n),nrow=m)
lambda = 1
starting_step_size = 0.2
tol = 0.001
old_loss = get_convex_loss(y, X_est, A, J)
new_loss = old_loss - tol

while (old_loss - new_loss >= tol) {
  #gradient descent
  gradient = get_gradient(y, X_est, A, J, m, n)
  step_size = get_step_size(y, X_est, A, J, gradient, starting_step_size, lambda)
  
  Gt_x = get_Gt_x(X_est, step_size, gradient, lambda)
  X_est = X_est - step_size * Gt_x
  
  old_loss = new_loss
  new_loss = get_convex_loss(y, X_est, A, J)
  print(new_loss)
}
