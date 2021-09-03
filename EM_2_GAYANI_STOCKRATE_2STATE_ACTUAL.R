STT_ACTUAL = read.csv(file.choose(),header = F)
STT_ACTUAL = STT_ACTUAL$V1
EM2_EXPECTATION_STOCKRATE=function(yt,TRUE_COEF)  # Function Name (Data,True_Parameters)
{
  theta = matrix(0, nrow = 2,ncol=6,byrow = T)
  #TRUE_COEF = c(0.6020,0.5148,0.044,0.2111,0.0250,0.0158,0.02,0.05)
  theta[1,]= c(TRUE_COEF)
  
  S_t <- matrix(0, nrow = length(yt), ncol = 2, byrow = T)
  S_tplus1 <- matrix(0, nrow = length(yt), ncol = 2, byrow = T)
  S_tT <- matrix(0, nrow = length(yt), ncol = 2, byrow = T)
  F_Like <- matrix(0, nrow = length(yt), ncol = 2, byrow = T)
  S_t[1,] = c(0.5,0.5)
  k=1
  for (i in 2:length(yt))
  {
    S_tplus1[i,1] = (1-theta[k,5])* S_t[i-1,1] + theta[k,6] * S_t[i-1,2]
    S_tplus1[i,2] = theta[k,5]*S_t[i-1,1] + (1-theta[k,6]) * S_t[i-1,2] 
    F_Like[i,1] = 1/sqrt(2*pi*theta[k,3]*yt[i-1]^2)*exp(-(yt[i]- yt[i-1]- yt[i-1]*theta[k,1])^2/(2*theta[k,3]*yt[i-1]^2)) #dnorm(yt[1], r0 + theta[k,1]*(theta[k,3]-r0)*dt , theta[k,5]*sqrt(dt) )
    F_Like[i,2]= 1/sqrt(2*pi*theta[k,4]*yt[i-1]^2)*exp(-(yt[i]- yt[i-1]- yt[i-1]*theta[k,2])^2/(2*theta[k,4]*yt[i-1]^2)) 
    S_t[i,1] = F_Like[i,1] * S_tplus1[(i),1] / sum(F_Like[i,1]* S_tplus1[(i),1], F_Like[i,2] * S_tplus1[(i),2])
    S_t[i,2] = F_Like[i,2] * S_tplus1[(i),2] / sum(F_Like[i,1]* S_tplus1[(i),1], F_Like[i,2] * S_tplus1[(i),2])
  }
  
  #applying kim's smoothing , backward filtering
  
  S_tT[length(yt),1] = S_t[length(yt),1]
  S_tT[length(yt),2] = S_t[length(yt),2]
  
  for (i in (length(yt)-1):1)
  {
    S_tT[i,1] <- S_t[i,1] * (1-theta[k,5])* S_tT[i+1,1]/ S_tplus1[i+1,1]+ S_t[i,1]* theta[k,5] * S_tT[i+1,2]/ S_tplus1[i+1,2]
    S_tT[i,2] <- S_t[i,2] * theta[k,6]* S_tT[i+1,1]/ S_tplus1[i+1,1] + S_t[i,2]* (1-theta[k,6]) * S_tT[i+1,2]/ S_tplus1[i+1,2]
  }
  S_tT=round(S_tT,10)
  
  
  #################################################################################################
  trans.mat <- matrix(0,nrow = length(yt)-1, ncol = 4, byrow = T)
  pij.matrix <- matrix(0,nrow = 1,ncol = 4)
  
  for (i in 1:(length(yt)-1))
  {
    trans.mat[i,1] <- S_t[i,1]* (S_tT[i+1,1] / S_tplus1[i+1,1]) * (1-theta[k,5])
    trans.mat[i,2] <- S_t[i,1]* (S_tT[i+1,2] / S_tplus1[i+1,2]) * (theta[k,5])
    trans.mat[i,3] <- S_t[i,2]* (S_tT[i+1,1] / S_tplus1[i+1,1]) * (theta[k,6])
    trans.mat[i,4] <- S_t[i,2]* (S_tT[i+1,2] / S_tplus1[i+1,2]) * (1-theta[k,6])
  }
  
  trans.sum <- colSums(trans.mat) #trans.sum[1]+trans.sum[3]
  #cc_next_state_1 <- sum(S_tT[,1]) - S_tT[1,1]
  #cc_next_state_2 <- sum(S_tT_rev[,2]) - S_tT_rev[1,2]
  
  pij.matrix[1,1] <- trans.sum[1]/(trans.sum[1]+trans.sum[3])
  pij.matrix[1,2] <- trans.sum[2]/(trans.sum[2]+trans.sum[4])
  pij.matrix[1,3] <- trans.sum[3]/(trans.sum[1]+trans.sum[3])
  pij.matrix[1,4] <- trans.sum[4]/(trans.sum[2]+trans.sum[4])
  #invariant.probs <- S_tT_rev[1,]
  ##############################################################################################################
  
  F_Like_1=as.vector(F_Like[,1]);F_Like_2=as.vector(F_Like[,1]);S_tplus1_1=as.vector(S_tplus1[,1]);S_tplus1_2=as.vector(S_tplus1[,2]);S_t_1=as.vector(S_t[,1]);S_t_2=as.vector(S_t[,2]);S_tT_1=as.vector(S_tT[,1]); S_tT_2=as.vector(S_tT[,2]);
  probs=data.frame(F_Like_1,F_Like_2,S_tplus1_1,S_tplus1_2,S_t_1,S_t_2,S_tT_1,S_tT_2)
  loglikelihood=sum(log(F_Like[2:length(yt),1])*S_tT[2:length(yt),1]+log(F_Like[2:length(yt),2])*S_tT[2:length(yt),2])
  probabilities=list(probs,loglikelihood)
  probabilities
}


yt= STT_ACTUAL
results=c(0.001,0.002,0.00002,0.00002,0.02,0.05)
N=1
epsilon=10
while ( N< 9) 
{
  
  E = EM2_EXPECTATION_STOCKRATE(yt,results)
  log_1 = E[[2]]
  
  ### B1 and B2 values (ALPHA VALUE)
  B1 = matrix(0, nrow = length(yt), ncol = 2)
  B2 = matrix(0, nrow = length(yt), ncol = 2)
  
  for (i in 2:length(yt))
  {
    B1[i,1]= yt[i] - yt[i-1]- sum((E[[1]]$S_tT_1[2:length(yt)])*(yt[2:length(yt)]-yt[1:length(yt)-1]))/sum(E[[1]]$S_tT_1[2:length(yt)])
    B1[i,2]= yt[i] - yt[i-1]- sum((E[[1]]$S_tT_2[2:length(yt)])*(yt[2:length(yt)]-yt[1:length(yt)-1]))/sum(E[[1]]$S_tT_2[2:length(yt)])
    B2[i,1]= sum(E[[1]]$S_tT_1[2:length(yt)]*yt[1:length(yt)-1])/sum(E[[1]]$S_tT_1[2:length(yt)]) - yt[i-1]
    B2[i,2]= sum(E[[1]]$S_tT_2[2:length(yt)]*yt[1:length(yt)-1])/sum(E[[1]]$S_tT_2[2:length(yt)]) - yt[i-1]
  }
  
  #K_1 = sum(E[[1]]$S_tT_1[2:length(yt)]*yt[1:length(yt)-1]*B1[2:length(yt),1])/sum(E[[1]]$S_tT_1[2:length(yt)]*yt[1:length(yt)-1]*B2[2:length(yt),1])
  #K_2 = sum(E[[1]]$S_tT_2[2:length(yt)]*yt[1:length(yt)-1]*B1[2:length(yt),2])/sum(E[[1]]$S_tT_2[2:length(yt)]*yt[1:length(yt)-1]*B2[2:length(yt),2])
  K_1= sum((E[[1]]$S_tT_1[2:length(yt)])*((yt[2:length(yt)]-yt[1:length(yt)-1])/yt[1:length(yt)-1]))/sum(E[[1]]$S_tT_1[2:length(yt)])
  K_2= sum((E[[1]]$S_tT_2[2:length(yt)])*((yt[2:length(yt)]-yt[1:length(yt)-1])/yt[1:length(yt)-1]))/sum(E[[1]]$S_tT_2[2:length(yt)])
  T_1= 0 #sum(E[[1]]$S_tT_1[2:length(yt)]*(yt[2:length(yt)]-(1-K_1)*yt[1:length(yt)-1]))/K_1/sum(E[[1]]$S_tT_1[2:length(yt)])
  T_2= 0 #sum(E[[1]]$S_tT_2[2:length(yt)]*(yt[2:length(yt)]-(1-K_2)*yt[1:length(yt)-1]))/K_2/sum(E[[1]]$S_tT_2[2:length(yt)])
  SI_1=sum(E[[1]]$S_tT_1[2:length(yt)]*((yt[2:length(yt)]- yt[1:length(yt)-1] - K_1*yt[1:length(yt)-1])^2))/sum(E[[1]]$S_tT_1[2:length(yt)]*yt[1:length(yt)-1]^2)
  SI_2=sum(E[[1]]$S_tT_2[2:length(yt)]*((yt[2:length(yt)]- yt[1:length(yt)-1] - K_2*yt[1:length(yt)-1])^2))/sum(E[[1]]$S_tT_2[2:length(yt)]*yt[1:length(yt)-1]^2)
  
  pi_11=sum(E[[1]]$S_tT_1[2:length(yt)]*(1-results[5])*E[[1]]$S_t_1[1:length(yt)-1]/E[[1]]$S_tplus1_1[2:length(yt)])/sum(E[[1]]$S_tT_1[1:length(yt)-1])
  pi_12=1-pi_11
  pi_22=sum(E[[1]]$S_tT_2[2:length(yt)]*(1-results[6])*E[[1]]$S_t_2[1:length(yt)-1]/E[[1]]$S_tplus1_2[2:length(yt)])/sum(E[[1]]$S_tT_2[1:length(yt)-1])
  pi_21=1-pi_22
  
  results=c(K_1,K_2,SI_1,SI_2,pi_12,pi_21)
  E=EM2_EXPECTATION_STOCKRATE(yt,results)
  log_2=E[[2]]
  epsilon=abs(log_1-log_2)
  N=N+1
}
