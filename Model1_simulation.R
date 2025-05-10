library(splines)
library(MASS)

source('fpca_ARL.R')
source('mfpca_ARL.R')
source('top_delta_ARL.R')
source('control limits train.R')


Gen_X<-function(m){
  l=seq(1,20,length.out = 20)
  n=length(l)
  Sigma=matrix(data=NA, nrow=n, ncol=n)
  for(i in 1:n){
    for (j in 1:n){
      Sigma[i,j]<-0.8^abs(i-j)
    }
  }
  vk <-matrix(bs(1:50,df = 6,degree = 2),nrow = 50,ncol = 6)
  X <-array(0,dim = c(m,50,20))
  for (i in 1:m){
    score_k <- mvrnorm(6,mu = rep(0,20),Sigma = Sigma)
    ksi_k<-score_k #whether indicator
    X[i,,] = vk%*%ksi_k + matrix(rnorm(1000,0,0.04),nrow = 50,ncol = 20)
  }
  list(data = X,basis = vk) #add score table maybe,add vk for check
}

set.seed(513) 
X_param<-Gen_X(200)
mf_est<-mfpca_est(X_param$data,d0 = 6)
p<-dim(mf_est$scores)[3] #number of sensors
dpca_list <- vector(mode = "list", length = p)
for (i in 1:p){
  #store parameter list of each sensor individually into a list
  dpca_list[[i]] <- dpca_est(X_param$data[,,i],d0 = 6) 
}

#Estimate the between profile correlation for MFPCA
Covb_mfpca <- array(0,dim = c(p,p))
#aggregate all the d = 1,2,...d0 for the covariance matrix 
for (i in 1:mf_est$d){
  Covb_mfpca = Covb_mfpca + mf_est$covb[i,,]
}
#scaling step for MFPCA
E1 <- solve(diag(sqrt(diag(Covb_mfpca)))) 
#between-correlation matrix after scaling
Corr <- E1%*%Covb_mfpca%*%E1 
umin <- 1.5

#scaling process for FPCA for having unit standard deviation
E <-rep(NA,length = p)
for(j in 1:p){
  Covb = sum(dpca_list[[j]]$covb)
  E[j] <- 1/sqrt(Covb)
}

#estimating in-control L
X_train<-Gen_X(20000)$data
L_8<-L_train(X_train,obnum = 8)
L_6<-L_train(X_train,obnum = 6)
L_4<-L_train(X_train,obnum = 4)

#main function 
fillans<-function(L,obnum){
  #group levels of mean-shift
  oc_run1 <-array(NA,dim = c(100,5))
  oc_run2 <-array(NA,dim = c(100,5))
  oc_run3 <-array(NA,dim = c(100,5))
  oc_run4 <-array(NA,dim = c(100,5))
  oc_run5 <-array(NA,dim = c(100,5))
  
  for(i in 1:100){
    for(mu_degree in c(0.7,1.0,1.4,2.1,2.8)){ #the magnitude of mean-shifts 
      X_oc <- Gen_X(300)$data
      #demonstrate partial shift on Model 2
      X_oc[,11:40,3] <- X_oc[,11:40,3]+mu_degree*array(1,dim = c(300,30))
      X_oc[,11:40,4] <- X_oc[,11:40,4]+mu_degree*array(1,dim = c(300,30))
      X_oc[,11:40,5] <- X_oc[,11:40,5]+mu_degree*array(1,dim = c(300,30))
      X_oc[,11:40,8] <- X_oc[,11:40,8]+mu_degree*array(1,dim = c(300,30))
      X_oc[,11:40,13] <- X_oc[,11:40,13]+mu_degree*array(1,dim = c(300,30))
      
      if (mu_degree == 0.7){ 
        oc_run1[i,1] <- Get_ARL(X_oc,Corr,L[1],obnum)$RL
        oc_run1[i,2] <- Get_ARLi(X_oc,L[2],obnum)$RL
        oc_run1[i,3] <- Get_ARLi_delta(X_oc,L[3], obnum,del = 0.01, q = 4)$RL
        oc_run1[i,4] <- Get_ARLi_delta(X_oc,L[4], obnum,del = 0.05, q = 4)$RL
        oc_run1[i,5] <- Get_ARLi_delta(X_oc,L[5], obnum,del = 0.1, q = 4)$RL}
      
      else if (mu_degree == 1.0){
        oc_run2[i,1] <- Get_ARL(X_oc,Corr,L[1],obnum)$RL
        oc_run2[i,2] <- Get_ARLi(X_oc,L[2],obnum)$RL
        oc_run2[i,3] <- Get_ARLi_delta(X_oc,L[3], obnum,del = 0.01, q = 4)$RL
        oc_run2[i,4] <- Get_ARLi_delta(X_oc,L[4], obnum,del = 0.05, q = 4)$RL
        oc_run2[i,5] <- Get_ARLi_delta(X_oc,L[5], obnum,del = 0.1, q = 4)$RL}
      
      else if (mu_degree == 1.4){
        oc_run3[i,1] <- Get_ARL(X_oc,Corr,L[1],obnum)$RL
        oc_run3[i,2] <- Get_ARLi(X_oc,L[2],obnum)$RL
        oc_run3[i,3] <- Get_ARLi_delta(X_oc,L[3], obnum,del = 0.01,q = 4)$RL
        oc_run3[i,4] <- Get_ARLi_delta(X_oc,L[4], obnum,del = 0.05, q = 4)$RL
        oc_run3[i,5] <- Get_ARLi_delta(X_oc,L[5], obnum,del = 0.1, q = 4)$RL}
      
      else if (mu_degree == 2.1){
        oc_run4[i,1] <- Get_ARL(X_oc,Corr,L[1],obnum)$RL
        oc_run4[i,2] <- Get_ARLi(X_oc,L[2],obnum)$RL
        oc_run4[i,3] <- Get_ARLi_delta(X_oc,L[3], obnum,del = 0.01, q = 4)$RL
        oc_run4[i,4] <- Get_ARLi_delta(X_oc,L[4], obnum,del = 0.05, q = 4)$RL
        oc_run4[i,5] <- Get_ARLi_delta(X_oc,L[5], obnum,del = 0.1, q = 4)$RL}
      
      else {
        oc_run5[i,1] <- Get_ARL(X_oc,Corr,L[1],obnum)$RL
        oc_run5[i,2] <- Get_ARLi(X_oc,L[2],obnum)$RL
        oc_run5[i,3] <- Get_ARLi_delta(X_oc,L[3], obnum,del = 0.01, q= 4)$RL
        oc_run5[i,4] <- Get_ARLi_delta(X_oc,L[4], obnum,del = 0.05, q= 4)$RL
        oc_run5[i,5] <- Get_ARLi_delta(X_oc,L[5], obnum,del = 0.1, q= 4)$RL}
    }
  }
  ansarray<-array(0,dim = c(10,5))
  ansarray[1,]<-apply(oc_run1,2,mean);ansarray[2,]<-apply(oc_run1,2,sd)
  ansarray[3,]<-apply(oc_run2,2,mean);ansarray[4,]<-apply(oc_run2,2,sd)
  ansarray[5,]<-apply(oc_run3,2,mean);ansarray[6,]<-apply(oc_run3,2,sd)
  ansarray[7,]<-apply(oc_run4,2,mean);ansarray[8,]<-apply(oc_run4,2,sd)
  ansarray[9,]<-apply(oc_run5,2,mean);ansarray[10,]<-apply(oc_run5,2,sd)
  ansarray
}

#main function execution 
output<-array(NA,dim = c(30,5))
output[1:10,]<-fillans(L_8,obnum = 8) 
output[11:20,]<-fillans(L_6,obnum = 6) 
output[21:30,]<-fillans(L_4,obnum = 4) 

out_se<-output
for (i in 1:15){
  out_se[2*i,]<-output[2*i,]/sqrt(100) #even line are standard errors
}
#final output in out_se, take about 20 minutes for out_se to get ready
rown<-c('0.7','','1.0','','1.4','','2.1','','2.8','')
rownames(out_se)<-c(rown,rown,rown)
colnames(out_se)<-c('proposed','fpca','fpca001','fpca005','fpca01')

#For line plot plot 
proposed<-out_se[seq(1,29,by = 2),1]
FPCA<-out_se[seq(1,29,by = 2),2]
FPCA001<-out_se[seq(1,29,by = 2),3]
FPCA005<-out_se[seq(1,29,by = 2),4]
FPCA01<-out_se[seq(1,29,by = 2),5]




