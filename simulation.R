set.seed (41792)

#Create and populate 4x10000 matrix with normal random variables
Z = matrix( rnorm(40000,mean=0,sd=1), 4, 10000)

#Create contrast matrix
A = matrix(c(1, 3, 0, 
             0, 1, 3, 
             1, 3, 0.5, 
             1.5, 1, 3), 
          nrow=3,
          ncol=4,
          byrow=TRUE)

#Create design matrix
X = matrix(c(0,0,0,0,
             1,0,0,0,
             1,1,0,0,
             1,2,0,0,
             1,3,0,0,
             1,3,1,0,
             1,3,1,1,
             1,3,1,2,
             1,3,1,3),
           nrow=9,
           ncol=4,
           byrow=TRUE)

ID = diag(4)
V1 = c(1,1,1,1)

simn = dim(Z)[2]

#Make this dynamic
rhos = matrix(c(0.6,0.4,0.7,0.5,0.5,0.3),
              nrow=3,
              ncol=2)
tau.mat = matrix(c(0.6,0.12,0,0,0.6,0.12,-0.005,-0.01),
               nrow=2,
               ncol=4,
               byrow=TRUE)
samp = c(200,300,400,500,600)

#Create matrices to hold simulated data
sims = matrix(nrow=dim(rhos)[1]*dim(tau.mat)[1]*length(samp)*simn,ncol=22)
calc = matrix(nrow=dim(rhos)[1]*dim(tau.mat)[1]*length(samp),ncol=22)


for (p in 1:dim(rhos)[1]){
  
  rho1 = rhos[p,1]
  rho2 = rhos[p,2]
  rho3 = rho1
  rho4 = rho2
  
  Rtilda = ((1-rho3)*ID) + (rho3*V1 %*% t(V1))

  #Create correlation matrix
  R = rbind(cbind(1,rho1*t(V1),rho2*t(V1)),
            cbind(rho1*V1,Rtilda,rho4*V1 %*% t(V1)),
            cbind(rho2*V1,rho4*V1 %*% t(V1),Rtilda))
  
  XR = solve(t(X) %*% solve(R) %*% X)
  XR[XR<=10**-9] = 0
  
  XR.T = chol(XR)
  
  
  for (t in 1:dim(tau.mat)[1]){
    tau = tau.mat[t,]  
    
    for (n in 1:length(samp)){
      
      #Variance
      v.a = (2/samp[n]) * A %*% t(solve(XR.T)) %*% solve(XR.T) 
      v.tau = (2/samp[n]) * t(solve(XR.T)) %*% solve(XR.T) 
    }
 
  }
}
mymat

