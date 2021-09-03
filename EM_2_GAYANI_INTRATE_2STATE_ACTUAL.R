EM2_EXPECTATION_INTERESTRATE=function(yt,TRUE_COEF)  # Function Name (Data,True_Parameters)
{
  theta = matrix(0, nrow = 2 , ncol = 8,byrow = T)
  #TRUE_COEF = c(0.4,0.2,  0.05238 , 0.00165628,  0.02313579, 0.00073605,  0.05,  0.05)
  theta[1,]= c(TRUE_COEF)
  
  S_t <- matrix(0, nrow = length(yt), ncol = 2, byrow = T)
  S_tplus1 <- matrix(0, nrow = length(yt), ncol = 2, byrow = T)
  S_tT <- matrix(0, nrow = length(yt), ncol = 2, byrow = T)
  F_Like <- matrix(0, nrow = length(yt), ncol = 2, byrow = T)
  S_t[1,] = c(0.5,0.5)
  k=1
  for (i in 2:length(yt))
  {
    S_tplus1[i,1] = (1-theta[k,7])* S_t[i-1,1] + theta[k,8] * S_t[i-1,2]
    S_tplus1[i,2] = theta[k,7]*S_t[i-1,1] + (1-theta[k,8]) * S_t[i-1,2] 
    F_Like[i,1]=1/sqrt(2*pi*theta[k,5])*exp(-(yt[i]-(1-theta[k,1])*yt[i-1]-theta[k,1]*theta[k,3])^2/(2*theta[k,5])) #dnorm(yt[1], r0 + theta[k,1]*(theta[k,3]-r0)*dt , theta[k,5]*sqrt(dt) )
    F_Like[i,2]=1/sqrt(2*pi*theta[k,6])*exp(-(yt[i]-(1-theta[k,2])*yt[i-1]-theta[k,2]*theta[k,4])^2/(2*theta[k,6])) 
    S_t[i,1] = F_Like[i,1] * S_tplus1[(i),1] / sum(F_Like[i,1]* S_tplus1[(i),1], F_Like[i,2] * S_tplus1[(i),2])
    S_t[i,2] = F_Like[i,2] * S_tplus1[(i),2] / sum(F_Like[i,1]* S_tplus1[(i),1], F_Like[i,2] * S_tplus1[(i),2])
  }
  
  #applying kim's smoothing , backward filtering
  
  S_tT[length(yt),1] = S_t[length(yt),1]
  S_tT[length(yt),2] = S_t[length(yt),2]
  
  for (i in (length(yt)-1):1)
  {
    S_tT[i,1] <- S_t[i,1] * (1-theta[k,7])* S_tT[i+1,1]/ S_tplus1[i+1,1]+ S_t[i,1]* theta[k,7] * S_tT[i+1,2]/ S_tplus1[i+1,2]
    S_tT[i,2] <- S_t[i,2] * theta[k,8]* S_tT[i+1,1]/ S_tplus1[i+1,1] + S_t[i,2]* (1-theta[k,8]) * S_tT[i+1,2]/ S_tplus1[i+1,2]
  }
  S_tT=round(S_tT,12)
  
  
  #################################################################################################
  trans.mat <- matrix(0,nrow = length(yt)-1, ncol = 4, byrow = T)
  pij.matrix <- matrix(0,nrow = 1,ncol = 4)
  
  for (i in 1:(length(yt)-1))
  {
    trans.mat[i,1] <- S_t[i,1]* (S_tT[i+1,1] / S_tplus1[i+1,1]) * (1-theta[k,7])
    trans.mat[i,2] <- S_t[i,1]* (S_tT[i+1,2] / S_tplus1[i+1,2]) * (theta[k,7])
    trans.mat[i,3] <- S_t[i,2]* (S_tT[i+1,1] / S_tplus1[i+1,1]) * (theta[k,8])
    trans.mat[i,4] <- S_t[i,2]* (S_tT[i+1,2] / S_tplus1[i+1,2]) * (1-theta[k,8])
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

INTEREST_RATE_DATA = read.csv(file.choose(), header = T)
INTEREST_RATE_DATA= INTEREST_RATE_DATA$Rate
INTEREST_RATE_DATA= INTEREST_RATE_DATA/100
INTEREST_RATE_DATA[1]
ts.plot(yt)
yt=INTEREST_RATE_DATA
results=c(0.5,0.2,  0.05238 , 0.00165628,  0.02313579, 0.023605,  0.05,  0.05) #results = c(0.1179, 0.1499, 0.0491 ,0.0345,0.0214^2, 0.0213^2,0.02,0.05) ##c(-0.1,-0.1,0.0220,0.02,0.03^2,0.04^2,0.02,0.05)
N=1

epsilon=10
while ( epsilon > 0.5  && N < 20)
{
  
  E = EM2_EXPECTATION_INTERESTRATE(yt,results)
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
  
  K_1 = sum(E[[1]]$S_tT_1[2:length(yt)]*yt[1:length(yt)-1]*B1[2:length(yt),1])/sum(E[[1]]$S_tT_1[2:length(yt)]*yt[1:length(yt)-1]*B2[2:length(yt),1])
  K_2 = sum(E[[1]]$S_tT_2[2:length(yt)]*yt[1:length(yt)-1]*B1[2:length(yt),2])/sum(E[[1]]$S_tT_2[2:length(yt)]*yt[1:length(yt)-1]*B2[2:length(yt),2])
  T_1=sum(E[[1]]$S_tT_1[2:length(yt)]*(yt[2:length(yt)]-(1-K_1)*yt[1:length(yt)-1]))/K_1/sum(E[[1]]$S_tT_1[2:length(yt)])
  T_2=sum(E[[1]]$S_tT_2[2:length(yt)]*(yt[2:length(yt)]-(1-K_2)*yt[1:length(yt)-1]))/K_2/sum(E[[1]]$S_tT_2[2:length(yt)])
  SI_1=sum(E[[1]]$S_tT_1[2:length(yt)]*(yt[2:length(yt)]-K_1*T_1-(1-K_1)*yt[1:length(yt)-1])^2)/sum(E[[1]]$S_tT_1[2:length(yt)])
  SI_2=sum(E[[1]]$S_tT_2[2:length(yt)]*(yt[2:length(yt)]-K_2*T_2-(1-K_2)*yt[1:length(yt)-1])^2)/sum(E[[1]]$S_tT_2[2:length(yt)])
  
  pi_11=sum(E[[1]]$S_tT_1[2:length(yt)]*(1-results[7])*E[[1]]$S_t_1[1:length(yt)-1]/E[[1]]$S_tplus1_1[2:length(yt)])/sum(E[[1]]$S_tT_1[1:length(yt)-1])
  pi_12=1-pi_11
  pi_22=sum(E[[1]]$S_tT_2[2:length(yt)]*(1-results[8])*E[[1]]$S_t_2[1:length(yt)-1]/E[[1]]$S_tplus1_2[2:length(yt)])/sum(E[[1]]$S_tT_2[1:length(yt)-1])
  pi_21=1-pi_22
  
  results=c(K_1,K_2,T_1,T_2,SI_1,SI_2,pi_12,pi_21)
  E=EM2_EXPECTATION_INTERESTRATE(yt,results)
  log_2=E[[2]]
  epsilon=abs(log_1-log_2)
  N=N+1
}

##################################################################################################################
AIC = -2*log_2 + 2*8  
AIC
BIC = -2*log_2 + 8*log(length(yt))
BIC








