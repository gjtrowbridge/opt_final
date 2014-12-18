#Generate the data

chain_network_omega = function(p) {
  omega = diag(rep(2, p))
  for (i in 1:(p-1)) {
    omega[i, i+1] = 1
    omega[i+1, i] = 1
  }
  
  return(omega)
}

nn_network_omega = function(p, rho=0, m=3) {
  #gen p points on unit square
  #calculate all p(p-1)/2 distances
  #find m nearest neighbors of each pt and link them  
  xloc<-runif(p,-1,1)
  yloc<-runif(p,-1,1)
  dists<-matrix(NA, nrow=p,ncol=p)
  
  for (i in 1:p) {
    for (j in 1:p) {
      dists[i,j]<-sqrt((xloc[i]-xloc[j])^2+(yloc[i]-yloc[j])^2)
    }
  }
  
  omega<-matrix(0,nrow=p,ncol=p)
  for (i in 1:p) {
    hold<-sort.int(dists[,i], index.return=TRUE)$ix
    for (j in 1:m) {
      omega[hold[j],i]<-1
      omega[i,hold[j]]<-1
    }
    omega[i,i]<-m #this makes it invertible. 
  }
  
  #add individual links as in chain
  if (rho!=0) {
    omega<-gen_add_links(omega,rho)
  }
  
  return(round(omega,10))
}

sf_network_omega = function(p, rho=0, m=3) {
  #use Barabasi-Albert algorithm (1999 paper)
  #it's included in igraph package... 
  #Details are given there if someone wants to write it from scratch
  require(igraph)
  g<-barabasi.game(p, directed=FALSE)
  omega<-as.matrix(get.adjacency(g))
  diag(omega)<-p/2 #to make it invertible
  
  # can add links here
  if (rho != 0) {
    omega<-gen_add_links(omega,rho)
  }
  
  
  return(round(omega,10))
}


#Data generation sample
overlaps = 5
for (omega_function in c(chain_network_omega, nn_network_omega, sf_network_omega)) {
  for (p in c(150, 25)) {
    
    for (i in 1:20) {
      omega1 = omega_function(p)
      omega2 = omega1
      
      #Make the networks similar
      for (j in 1:overlaps) {
        to_swap = sample(1:p, 2, replace=FALSE)
        omega2[to_swap[1], to_swap[2]] = 1 - omega2[to_swap[1], to_swap[2]]
      }
    }
  }
}


