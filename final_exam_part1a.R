part_1a = function(x) {
  u = rep(0, length(x))
  for (i in 1:length(x)) {
    if (x[i] > 1) {
      u[i] = min(1, x[i] - 1)
    } else if (x[i] < -1) {
      u[i] = max(-1, x[i] + 1)
    } else {
      u[i] = 0
    }
  }
  return(u)
}