rtt_actual = read.csv(file.choose(),header = T)
rtt_yt=rtt_actual$Rate/100
EM2_EXPECTATION_INTERESTRATE_4STATE=function(yt,TRUE_COEF)  # Function Name (Data,True_Parameters)
{
  NS = 4
  theta = matrix(0, nrow = 10 , ncol = 12+16, byrow = T)
  # a1,a2,a3,b1,b2,b3,s1,s2,s3,p11,p12,p13, p21,p22, p23,p31,p32,p33
  #TRUE_COEF = c(0.2,0.4,0.6,0.8, 0.05, 0.10, 0.15, 0.2, 0.03^2, 0.03^2,0.03^2,0.03^2, 0.9, 0.07, 0.02,0.01,  0.07,0.9,0.02, 0.01, 0.02, 0.07,0.9,0.01, 0.01,0.02,0.07,0.9) 
  
  theta[1,]= TRUE_COEF
  S_t <- matrix(0, nrow = length(yt), ncol = NS, byrow = T)
  S_tplus1 <- matrix(0, nrow = length(yt), ncol = NS, byrow = T)
  S_tT <- matrix(0, nrow = length(yt), ncol = NS, byrow = T)
  F_Like <- matrix(0, nrow = length(yt), ncol = NS, byrow = T)
  S_t[1,] = c(rep(1/NS, NS))
  k=1
  for (i in 2:length(yt))
  {
    S_tplus1[i,1] = theta[k,13]* S_t[i-1,1] + theta[k,17]* S_t[i-1,2] + theta[k,21]* S_t[i-1,3] + theta[k,25]* S_t[i-1,4]
    S_tplus1[i,2] = theta[k,14]* S_t[i-1,1] + theta[k,18]* S_t[i-1,2] + theta[k,22]* S_t[i-1,3] + theta[k,26]* S_t[i-1,4] 
    S_tplus1[i,3] = theta[k,15]* S_t[i-1,1] + theta[k,19]* S_t[i-1,2] + theta[k,23]* S_t[i-1,3] + theta[k,27]* S_t[i-1,4]
    S_tplus1[i,4] = theta[k,16]* S_t[i-1,1] + theta[k,20]* S_t[i-1,2] + theta[k,24]* S_t[i-1,3] + theta[k,28]* S_t[i-1,4]
    F_Like[i,1]=1/sqrt(2*pi*theta[k,9])*exp(-(yt[i]-(1-theta[k,1])*yt[i-1]-theta[k,1]*theta[k,5])^2/(2*theta[k,9])) #dnorm(yt[1], r0 + theta[k,1]*(theta[k,3]-r0)*dt , theta[k,5]*sqrt(dt) )
    F_Like[i,2]=1/sqrt(2*pi*theta[k,10])*exp(-(yt[i]-(1-theta[k,2])*yt[i-1]-theta[k,2]*theta[k,6])^2/(2*theta[k,10]))
    F_Like[i,3]=1/sqrt(2*pi*theta[k,11])*exp(-(yt[i]-(1-theta[k,3])*yt[i-1]-theta[k,3]*theta[k,7])^2/(2*theta[k,11]))
    F_Like[i,4]=1/sqrt(2*pi*theta[k,12])*exp(-(yt[i]-(1-theta[k,4])*yt[i-1]-theta[k,4]*theta[k,8])^2/(2*theta[k,12]))
    S_t[i,1] = F_Like[i,1] * S_tplus1[(i),1] / sum(F_Like[i,1]* S_tplus1[(i),1], F_Like[i,2] * S_tplus1[(i),2] , F_Like[i,3]* S_tplus1[(i),3], F_Like[i,4]* S_tplus1[(i),4])
    S_t[i,2] = F_Like[i,2] * S_tplus1[(i),2] / sum(F_Like[i,1]* S_tplus1[(i),1], F_Like[i,2] * S_tplus1[(i),2] , F_Like[i,3]* S_tplus1[(i),3], F_Like[i,4]* S_tplus1[(i),4])
    S_t[i,3] = F_Like[i,3] * S_tplus1[(i),3] / sum(F_Like[i,1]* S_tplus1[(i),1], F_Like[i,2] * S_tplus1[(i),2] , F_Like[i,3]* S_tplus1[(i),3], F_Like[i,4]* S_tplus1[(i),4])
    S_t[i,4] = F_Like[i,4] * S_tplus1[(i),4] / sum(F_Like[i,1]* S_tplus1[(i),1], F_Like[i,2] * S_tplus1[(i),2] , F_Like[i,3]* S_tplus1[(i),3], F_Like[i,4]* S_tplus1[(i),4])
  }
  
  #applying kim's smoothing , backward filtering
  
  S_tT[length(yt),1] = S_t[length(yt),1]
  S_tT[length(yt),2] = S_t[length(yt),2]
  S_tT[length(yt),3] = S_t[length(yt),3]
  S_tT[length(yt),4] = S_t[length(yt),4]
  for (i in (length(yt)-1):1)
  {
    S_tT[i,1] <- S_t[i,1] * theta[k,13]* S_tT[i+1,1]/ S_tplus1[i+1,1] + S_t[i,1]* theta[k,14] * S_tT[i+1,2]/ S_tplus1[i+1,2] + S_t[i,1]* theta[k,15] * S_tT[i+1,3]/ S_tplus1[i+1,3] + S_t[i,1]* theta[k,16] * S_tT[i+1,4]/ S_tplus1[i+1,4]
    S_tT[i,2] <- S_t[i,2] * theta[k,17]* S_tT[i+1,1]/ S_tplus1[i+1,1] + S_t[i,2]* theta[k,18] * S_tT[i+1,2]/ S_tplus1[i+1,2] + S_t[i,2]* theta[k,19] * S_tT[i+1,3]/ S_tplus1[i+1,3] + S_t[i,2]* theta[k,20] * S_tT[i+1,4]/ S_tplus1[i+1,4]
    S_tT[i,3] <- S_t[i,3] * theta[k,21]* S_tT[i+1,1]/ S_tplus1[i+1,1] + S_t[i,3]* theta[k,22] * S_tT[i+1,2]/ S_tplus1[i+1,2] + S_t[i,3]* theta[k,23] * S_tT[i+1,3]/ S_tplus1[i+1,3] + S_t[i,3]* theta[k,24] * S_tT[i+1,4]/ S_tplus1[i+1,4]
    S_tT[i,4] <- S_t[i,4] * theta[k,25]* S_tT[i+1,1]/ S_tplus1[i+1,1] + S_t[i,4]* theta[k,26] * S_tT[i+1,2]/ S_tplus1[i+1,2] + S_t[i,4]* theta[k,27] * S_tT[i+1,3]/ S_tplus1[i+1,3] + S_t[i,4]* theta[k,28] * S_tT[i+1,4]/ S_tplus1[i+1,4]
  }
  S_tT=round(S_tT,12)
  
  
  #################################################################################################
  trans.mat <- matrix(0,nrow = length(yt)-1, ncol = 16, byrow = T)
  pij.matrix <- matrix(0,nrow = 1,ncol = 16)
  
  for (i in 1:(length(yt)-1))
  {
    trans.mat[i,1] <- S_t[i,1]* (S_tT[i+1,1] / S_tplus1[i+1,1]) * theta[k,13]
    trans.mat[i,2] <- S_t[i,1]* (S_tT[i+1,2] / S_tplus1[i+1,2]) * theta[k,14]
    trans.mat[i,3] <- S_t[i,1]* (S_tT[i+1,3] / S_tplus1[i+1,3]) * theta[k,15]
    trans.mat[i,4] <- S_t[i,1]* (S_tT[i+1,4] / S_tplus1[i+1,4]) * theta[k,16]
    trans.mat[i,5] <- S_t[i,2]* (S_tT[i+1,2] / S_tplus1[i+1,1]) * theta[k,17]
    trans.mat[i,6] <- S_t[i,2]* (S_tT[i+1,2] / S_tplus1[i+1,3]) * theta[k,18]
    trans.mat[i,7] <- S_t[i,2]* (S_tT[i+1,3] / S_tplus1[i+1,3]) * theta[k,19]
    trans.mat[i,8] <- S_t[i,2]* (S_tT[i+1,4] / S_tplus1[i+1,4]) * theta[k,20]
    trans.mat[i,9] <- S_t[i,3]* (S_tT[i+1,1] / S_tplus1[i+1,1]) * theta[k,21]
    trans.mat[i,10] <- S_t[i,3]* (S_tT[i+1,2] / S_tplus1[i+1,2]) * theta[k,22]
    trans.mat[i,11] <- S_t[i,3]* (S_tT[i+1,3] / S_tplus1[i+1,3]) * theta[k,23]
    trans.mat[i,12] <- S_t[i,3]* (S_tT[i+1,4] / S_tplus1[i+1,4]) * theta[k,24]
    trans.mat[i,13] <- S_t[i,4]* (S_tT[i+1,1] / S_tplus1[i+1,1]) * theta[k,25]
    trans.mat[i,14] <- S_t[i,4]* (S_tT[i+1,2] / S_tplus1[i+1,2]) * theta[k,26]
    trans.mat[i,15] <- S_t[i,4]* (S_tT[i+1,3] / S_tplus1[i+1,3]) * theta[k,27]
    trans.mat[i,16] <- S_t[i,4]* (S_tT[i+1,4] / S_tplus1[i+1,4]) * theta[k,28]
  }
  
  trans.sum <- colSums(trans.mat) #trans.sum[1]+trans.sum[3]
  cc_next_state_1 <- sum(S_tT[,1]) - S_tT[1,1]
  #cc_next_state_2 <- sum(S_tT_rev[,2]) - S_tT_rev[1,2]
  
  pij.matrix[1,1] <- trans.sum[1]/(trans.sum[1]+trans.sum[5]+trans.sum[9]+trans.sum[13])
  pij.matrix[1,2] <- trans.sum[2]/(trans.sum[2]+trans.sum[6]+trans.sum[10]+trans.sum[14])
  pij.matrix[1,3] <- trans.sum[3]/(trans.sum[3]+trans.sum[7]+trans.sum[11]+trans.sum[15])
  pij.matrix[1,4] <- trans.sum[4]/(trans.sum[4]+trans.sum[8]+trans.sum[12]+trans.sum[16])
  pij.matrix[1,5] <- trans.sum[5]/(trans.sum[1]+trans.sum[5]+trans.sum[9]+trans.sum[13])
  pij.matrix[1,6] <- trans.sum[6]/(trans.sum[2]+trans.sum[6]+trans.sum[10]+trans.sum[14])
  pij.matrix[1,7] <- trans.sum[7]/(trans.sum[3]+trans.sum[7]+trans.sum[11]+trans.sum[15])
  pij.matrix[1,8] <- trans.sum[8]/(trans.sum[4]+trans.sum[8]+trans.sum[12]+trans.sum[16])
  pij.matrix[1,9] <- trans.sum[9]/(trans.sum[1]+trans.sum[5]+trans.sum[9]+trans.sum[13])
  pij.matrix[1,10] <- trans.sum[10]/(trans.sum[2]+trans.sum[6]+trans.sum[10]+trans.sum[14])
  pij.matrix[1,11] <- trans.sum[11]/(trans.sum[3]+trans.sum[7]+trans.sum[11]+trans.sum[15])
  pij.matrix[1,12] <- trans.sum[12]/(trans.sum[4]+trans.sum[8]+trans.sum[12]+trans.sum[16])
  pij.matrix[1,13] <- trans.sum[13]/(trans.sum[1]+trans.sum[5]+trans.sum[9]+trans.sum[13])
  pij.matrix[1,14] <- trans.sum[14]/(trans.sum[2]+trans.sum[6]+trans.sum[10]+trans.sum[14])
  pij.matrix[1,15] <- trans.sum[15]/(trans.sum[3]+trans.sum[7]+trans.sum[11]+trans.sum[15])
  pij.matrix[1,16] <- trans.sum[16]/(trans.sum[4]+trans.sum[8]+trans.sum[12]+trans.sum[16])
  
  #invariant.probs <- S_tT_rev[1,]
  ##############################################################################################################
  
  F_Like_1=as.vector(F_Like[,1]);F_Like_2=as.vector(F_Like[,1]);F_Like_3=as.vector(F_Like[,3]);F_Like_4=as.vector(F_Like[,4]);S_tplus1_1=as.vector(S_tplus1[,1]);S_tplus1_2=as.vector(S_tplus1[,2]);S_tplus1_3=as.vector(S_tplus1[,3]);S_tplus1_4=as.vector(S_tplus1[,4]);S_t_1=as.vector(S_t[,1]);S_t_2=as.vector(S_t[,2]);S_t_3=as.vector(S_t[,3]);S_t_4=as.vector(S_t[,4]);S_tT_1=as.vector(S_tT[,1]); S_tT_2=as.vector(S_tT[,2]);S_tT_3=as.vector(S_tT[,3]);S_tT_4=as.vector(S_tT[,4]);
  probs=data.frame(F_Like_1,F_Like_2,F_Like_3,F_Like_4,S_tplus1_1,S_tplus1_2,S_tplus1_3,S_tplus1_4,S_t_1,S_t_2,S_t_3,S_t_4,S_tT_1,S_tT_2,S_tT_3,S_tT_4)
  loglikelihood=sum(log(F_Like[2:length(yt),1])*S_tT[2:length(yt),1],log(F_Like[2:length(yt),2])*S_tT[2:length(yt),2],log(F_Like[2:length(yt),3])*S_tT[2:length(yt),3],log(F_Like[2:length(yt),4])*S_tT[2:length(yt),4])
  probabilities=list(probs,loglikelihood)
  probabilities
}

L=1
PARAS = matrix(NA, nrow = L, ncol=28)
LIKELIHOOD = numeric(L)
EPSILON = numeric(L)

for (j in 1:L) {
  #results=c(0.6020,0.6020,0.5148,0.0440, 0.0440, 0.2111, 0.2111, 0.0250, 0.0158, 0.9200, 0.0500, 0.0300, 0.0500, 0.9200, 0.0300,0.0300, 0.0500, 0.9200)
  results=c(0.2,0.4,0.2,0.6, 0.05, 0.10, 0.1, 0.1, 0.03^2, 0.03^2,0.03^2,0.03^2, 0.9, 0.07, 0.02,0.01,  0.07,0.9,0.02, 0.01, 0.02, 0.07,0.9,0.01, 0.01,0.02,0.07,0.9) 
  N=1
  epsilon=10
  yt = rtt_yt # RESULTS_INTERESTRATE_4state[,j]
  
  while ( epsilon > 5 && N < 20) 
  {
    E = EM2_EXPECTATION_INTERESTRATE_4STATE(yt,results)
    log_1 = E[[2]]
    
    ### B1 and B2 values (ALPHA VALUE)
    B1 = matrix(0, nrow = length(yt), ncol = 4)
    B2 = matrix(0, nrow = length(yt), ncol = 4)
    
    for (i in 2:length(yt))
    {
      B1[i,1]= yt[i] - yt[i-1]- sum((E[[1]]$S_tT_1[2:length(yt)])*(yt[2:length(yt)]-yt[1:length(yt)-1]))/sum(E[[1]]$S_tT_1[2:length(yt)])
      B1[i,2]= yt[i] - yt[i-1]- sum((E[[1]]$S_tT_2[2:length(yt)])*(yt[2:length(yt)]-yt[1:length(yt)-1]))/sum(E[[1]]$S_tT_2[2:length(yt)])
      B1[i,3]= yt[i] - yt[i-1]- sum((E[[1]]$S_tT_3[2:length(yt)])*(yt[2:length(yt)]-yt[1:length(yt)-1]))/sum(E[[1]]$S_tT_3[2:length(yt)])
      B1[i,4]= yt[i] - yt[i-1]- sum((E[[1]]$S_tT_4[2:length(yt)])*(yt[2:length(yt)]-yt[1:length(yt)-1]))/sum(E[[1]]$S_tT_4[2:length(yt)])
      B2[i,1]= sum(E[[1]]$S_tT_1[2:length(yt)]*yt[1:length(yt)-1])/sum(E[[1]]$S_tT_1[2:length(yt)]) - yt[i-1]
      B2[i,2]= sum(E[[1]]$S_tT_2[2:length(yt)]*yt[1:length(yt)-1])/sum(E[[1]]$S_tT_2[2:length(yt)]) - yt[i-1]
      B2[i,3]= sum(E[[1]]$S_tT_3[2:length(yt)]*yt[1:length(yt)-1])/sum(E[[1]]$S_tT_3[2:length(yt)]) - yt[i-1]
      B2[i,4]= sum(E[[1]]$S_tT_4[2:length(yt)]*yt[1:length(yt)-1])/sum(E[[1]]$S_tT_4[2:length(yt)]) - yt[i-1]
    }
    
    K_1 = sum(E[[1]]$S_tT_1[2:length(yt)]*yt[1:length(yt)-1]*B1[2:length(yt),1])/sum(E[[1]]$S_tT_1[2:length(yt)]*yt[1:length(yt)-1]*B2[2:length(yt),1])
    K_2 = sum(E[[1]]$S_tT_2[2:length(yt)]*yt[1:length(yt)-1]*B1[2:length(yt),2])/sum(E[[1]]$S_tT_2[2:length(yt)]*yt[1:length(yt)-1]*B2[2:length(yt),2])
    K_3 = sum(E[[1]]$S_tT_3[2:length(yt)]*yt[1:length(yt)-1]*B1[2:length(yt),3])/sum(E[[1]]$S_tT_3[2:length(yt)]*yt[1:length(yt)-1]*B2[2:length(yt),3])
    K_4 = sum(E[[1]]$S_tT_4[2:length(yt)]*yt[1:length(yt)-1]*B1[2:length(yt),4])/sum(E[[1]]$S_tT_4[2:length(yt)]*yt[1:length(yt)-1]*B2[2:length(yt),4])
    T_1=sum(E[[1]]$S_tT_1[2:length(yt)]*(yt[2:length(yt)]-(1-K_1)*yt[1:length(yt)-1]))/K_1/sum(E[[1]]$S_tT_1[2:length(yt)])
    T_2=sum(E[[1]]$S_tT_2[2:length(yt)]*(yt[2:length(yt)]-(1-K_2)*yt[1:length(yt)-1]))/K_2/sum(E[[1]]$S_tT_2[2:length(yt)])
    T_3=sum(E[[1]]$S_tT_3[2:length(yt)]*(yt[2:length(yt)]-(1-K_3)*yt[1:length(yt)-1]))/K_3/sum(E[[1]]$S_tT_3[2:length(yt)])
    T_4=sum(E[[1]]$S_tT_4[2:length(yt)]*(yt[2:length(yt)]-(1-K_4)*yt[1:length(yt)-1]))/K_4/sum(E[[1]]$S_tT_4[2:length(yt)])
    SI_1=sum(E[[1]]$S_tT_1[2:length(yt)]*(yt[2:length(yt)]-K_1*T_1-(1-K_1)*yt[1:length(yt)-1])^2)/sum(E[[1]]$S_tT_1[2:length(yt)])
    SI_2=sum(E[[1]]$S_tT_2[2:length(yt)]*(yt[2:length(yt)]-K_2*T_2-(1-K_2)*yt[1:length(yt)-1])^2)/sum(E[[1]]$S_tT_2[2:length(yt)])
    SI_3=sum(E[[1]]$S_tT_3[2:length(yt)]*(yt[2:length(yt)]-K_3*T_3-(1-K_3)*yt[1:length(yt)-1])^2)/sum(E[[1]]$S_tT_3[2:length(yt)])
    SI_4=sum(E[[1]]$S_tT_4[2:length(yt)]*(yt[2:length(yt)]-K_4*T_4-(1-K_4)*yt[1:length(yt)-1])^2)/sum(E[[1]]$S_tT_4[2:length(yt)])
    
    pi_11=sum(E[[1]]$S_tT_1[2:length(yt)]*results[13]*E[[1]]$S_t_1[1:length(yt)-1]/E[[1]]$S_tplus1_1[2:length(yt)])/sum(E[[1]]$S_tT_1[1:length(yt)-1])
    pi_12=sum(E[[1]]$S_tT_2[2:length(yt)]*results[14]*E[[1]]$S_t_1[1:length(yt)-1]/E[[1]]$S_tplus1_2[2:length(yt)])/sum(E[[1]]$S_tT_1[1:length(yt)-1])
    pi_13=sum(E[[1]]$S_tT_3[2:length(yt)]*results[15]*E[[1]]$S_t_1[1:length(yt)-1]/E[[1]]$S_tplus1_3[2:length(yt)])/sum(E[[1]]$S_tT_1[1:length(yt)-1])
    pi_14= 1 - pi_11 - pi_12 - pi_13
    pi_21=sum(E[[1]]$S_tT_1[2:length(yt)]*results[17]*E[[1]]$S_t_2[1:length(yt)-1]/E[[1]]$S_tplus1_1[2:length(yt)])/sum(E[[1]]$S_tT_2[1:length(yt)-1])
    pi_22=sum(E[[1]]$S_tT_2[2:length(yt)]*results[18]*E[[1]]$S_t_2[1:length(yt)-1]/E[[1]]$S_tplus1_2[2:length(yt)])/sum(E[[1]]$S_tT_2[1:length(yt)-1])
    pi_23=sum(E[[1]]$S_tT_3[2:length(yt)]*results[19]*E[[1]]$S_t_2[1:length(yt)-1]/E[[1]]$S_tplus1_3[2:length(yt)])/sum(E[[1]]$S_tT_2[1:length(yt)-1])
    pi_24=1- pi_21 - pi_22 - pi_23
    pi_31=sum(E[[1]]$S_tT_1[2:length(yt)]*results[21]*E[[1]]$S_t_3[1:length(yt)-1]/E[[1]]$S_tplus1_1[2:length(yt)])/sum(E[[1]]$S_tT_3[1:length(yt)-1])
    pi_32=sum(E[[1]]$S_tT_2[2:length(yt)]*results[22]*E[[1]]$S_t_3[1:length(yt)-1]/E[[1]]$S_tplus1_2[2:length(yt)])/sum(E[[1]]$S_tT_3[1:length(yt)-1])
    pi_33=sum(E[[1]]$S_tT_3[2:length(yt)]*results[22]*E[[1]]$S_t_3[1:length(yt)-1]/E[[1]]$S_tplus1_3[2:length(yt)])/sum(E[[1]]$S_tT_3[1:length(yt)-1])
    pi_34=1- pi_31 - pi_32 - pi_33
    pi_41=sum(E[[1]]$S_tT_1[2:length(yt)]*results[25]*E[[1]]$S_t_4[1:length(yt)-1]/E[[1]]$S_tplus1_1[2:length(yt)])/sum(E[[1]]$S_tT_4[1:length(yt)-1])
    pi_42=sum(E[[1]]$S_tT_2[2:length(yt)]*results[26]*E[[1]]$S_t_4[1:length(yt)-1]/E[[1]]$S_tplus1_2[2:length(yt)])/sum(E[[1]]$S_tT_4[1:length(yt)-1])
    pi_43=sum(E[[1]]$S_tT_3[2:length(yt)]*results[27]*E[[1]]$S_t_4[1:length(yt)-1]/E[[1]]$S_tplus1_3[2:length(yt)])/sum(E[[1]]$S_tT_4[1:length(yt)-1])
    pi_44=1- pi_41 - pi_42 - pi_43
    
    results=c(K_1,K_2,K_3,K_4,T_1,T_2,T_3,T_4,SI_1,SI_2,SI_3,SI_4,pi_11,pi_12,pi_13,pi_14,pi_21,pi_22,pi_23,pi_24,pi_31,pi_32,pi_33,pi_34,pi_41,pi_42,pi_43,pi_44)
    E=EM2_EXPECTATION_INTERESTRATE_4STATE(yt,results)
    log_2=E[[2]]
    epsilon=abs(log_1-log_2)
    N=N+1
  }
  PARAS[j,] <- results
  LIKELIHOOD[j] <- log_2
  EPSILON[j] <- epsilon
}


PARAS[1,]


BIAS = colMeans(PARAS,na.rm = TRUE) - c(0.2,0.4,0.6,0.8, 0.05, 0.10, 0.15, 0.2, 0.03^2, 0.03^2,0.03^2,0.03^2, 0.9, 0.07, 0.02,0.01,  0.07,0.9,0.02, 0.01, 0.02, 0.07,0.9,0.01, 0.01,0.02,0.07,0.9) 
SD = c(sd(PARAS[1:20,1]) , sd(PARAS[1:20,2]) , sd(PARAS[1:20,3]) , sd(PARAS[1:20,4]) 
       ,sd(PARAS[1:20,5]) , sd(PARAS[1:20,6]) , sd(PARAS[1:20,7]) , sd(PARAS[1:20,8]),
       sd(PARAS[1:20,9]) , sd(PARAS[1:20,10]) , sd(PARAS[1:20,11]) , sd(PARAS[1:20,12]) 
       ,sd(PARAS[1:20,13]) , sd(PARAS[1:20,14]) , sd(PARAS[1:20,15]) , sd(PARAS[1:20,16]),
       sd(PARAS[1:20,17]) , sd(PARAS[1:20,18]), sd(PARAS[1:20,19]) , sd(PARAS[1:20,20]) , 
       sd(PARAS[1:20,21]) , sd(PARAS[1:20,22]) ,sd(PARAS[1:20,23]) , sd(PARAS[1:20,24]) , 
       sd(PARAS[1:20,25]) , sd(PARAS[1:20,26]), sd(PARAS[1:20,27]) , sd(PARAS[1:20,28]))
MSE = BIAS^2+ SD^2
