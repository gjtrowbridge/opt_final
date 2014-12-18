regular_loss = function(u, x) {
  return(0.5 * (u - x) %*% (u - x))
}

part1b = function(x) {
  get_max_info = function(u) {
    umax = max(u)
    umaxes_flagged = u == umax
    umaxes_count = sum(umaxes_flagged)
    if (umaxes_count == length(u)) {
      next_umax = -Inf
    } else {
      next_umax = max(u[!umaxes_flagged])
    }
    return(list(
      max = umax,
      maxes_flagged = umaxes_flagged,
      maxes_count = umaxes_count,
      next_max = next_umax
    ))
  }
  u = x
  #Derived from gradient of the proximal operator
  reduction_points = 1
  
  loss_from_max = max(u)
  
  while (reduction_points > 0) {
    max_info = get_max_info(u)
    reduce_max_by = max_info$max - max_info$next_max
    reduction_points_required = (reduce_max_by) * max_info$maxes_count
    
    #Reduce only by amount remaining
    if (reduction_points_required > reduction_points) {
      reduce_max_by = reduction_points / max_info$maxes_count
      reduction_points_required = reduction_points
    }
    
    reduction_points = reduction_points - reduction_points_required
    u[max_info$maxes_flagged] = u[max_info$maxes_flagged] - reduce_max_by
  }
  
  return(u)
}