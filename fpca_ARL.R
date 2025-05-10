#Get FPC scores for single sample in single sensor 
#X2 has 50 time points, but sample# = 1, and p = 1
dpca_score<-function(X2,efns,mua){ 
  demeaned <- X2-mua
  n <- dim(efns)[1] 
  d<- dim(efns)[2]
  scores <-rep(0,d)
  for(l in 1:d){
    scores[l] = sum(demeaned*efns[,l])/n
  }
  scores #return a d*1 matrix
}
#Estimate fpca parameters 
#input x have samples and time points but only 1 sensor
#d0 is the number of decomposition dimension.
dpca_est<-function(x,d0){
  m = dim(x)[1]
  n = dim(x)[2]
  cova<-array(0,dim=c(n,n))
  mean.hat<-array(0,dim=c(n,1))
  mean.hat <- colMeans(x)
  #construct covariance fucntion 
  for (i in 1:m){
    cova = cova + (x[i,]-mean.hat)%*%t(x[i,]-mean.hat)
  }
  cova = (1/m)*cova
  #Eigen decompostion 
  eio<-eigen(cova)
  eval<-eio$values/n
  evec<-eio$vectors*sqrt(n)
  #truncate for d = d0
  eval0<-eval[1:d0]  
  efns0<-evec[,1:d0]
  for(i in 1:d0){
    if (sum(efns0[1:10,i])<0){
      efns0[,i]<-efns0[,i]*(-1)
    }
  }
  demeaned <- x - t(matrix(rep(mean.hat, m),nrow=length(mean.hat)))
  #Get FPC score
  scores <- matrix(NA, nrow=m, ncol=d0)
  for(i in 1:m){
    for (l in 1:d0){
      scores[i,l] = sum(demeaned[i,]*efns0[,l])/n #**
    }
  }
  #get variance structure 
  covb<-array(0,dim = d0)
  for (l in 1:d0){
    for(i in 1:m){
      covb[l] = covb[l] + sum((x[i,]-mean.hat)*efns0[,l])*sum((x[i,]-mean.hat)*efns0[,l])/(n^2)#**
    }
    covb[l] = covb[l]/m
  }
  #output lsit includes: 
  #1.dimension, 2.template profile, 3.eigenvalue, 4.eigenfucntions, 5.fpc-scores and 6.variacne
  list(d=d0,mua=mean.hat,evalue = eval0,efns = efns0,scores = scores,covb = covb)
}


#FPCA + N(0,1) approach CUSUM procedure

#set up initial table
ARLi_initial<-function(testset,obnum){
cusum_score<-rep(0,length=p)
for (i in 1:obnum){
  ob_init<-dpca_score(testset[1,,i],dpca_list[[i]]$efns,dpca_list[[i]]$mua) #d*1
  cusum_score[i] <- sum(ob_init)*E[i] #stanardlized FPC score 
}

for(i in (obnum+1):p){
  ran_init <- rnorm(1,0,1) ## i.i.d N(0,1) 
  cusum_score[i] <- ran_init  
}
#Two sided CUSUM W_i = max(W_pos,W_neg)
table <- data.frame(sensor_id = seq(1,p,length.out = p),cusum_score)
for(i in 1:length(cusum_score)){
  table[i,"W_positive"]<-max(umin*cusum_score[i]-(umin^2)/2,0)
  table[i,"W_negative"]<-max(-umin*cusum_score[i]-(umin^2)/2,0)
  table[i,"W_i"]<-max(table[i,"W_positive"],table[i,"W_negative"])
}
#reorder local statistics
table["rank"] = rank(-table$W_i,ties.method = "random")
ob_sensor <- subset(table,rank<= obnum)$sensor_id
list(table = table,next_ob = ob_sensor)
}

#Get Running Length of 1 repetition 
#L is the control limits, obnum is the observation number 
Get_ARLi<-function(testset,L,obnum){
  procedure <-vector(mode ="list",length = 300)
  ob_hist<-array(0,dim = c(300,obnum))
  L_hist <- rep(0,dim(testset)[1]) #depend on train size
  init_list <- ARLi_initial(testset,obnum)
  ob_hist[1,] <- seq(1,obnum,length.out = obnum)
  ob_hist[2,] <- init_list$next_ob
  procedure[[1]]<-init_list$table
  s <- as.matrix(init_list$table[,"W_i"])
  #Global index 
  L_hist[1] <- sqrt(t(s)%*%solve(diag(length(s)))%*%s)  
  if(L_hist[1]>=L){
    return (list(RL = 1,ob_hist = ob_hist,procedure = procedure,L_hist = L_hist))
  }
  current <- init_list$table
  for(sample in 2:dim(testset)[1]){
    next_step <- cusum_df(testset,current,sample,obnum) #return the newest table
    next_table<- next_step$table
    s <- as.matrix(next_table[,"W_i"])
    #Global Index 
    L_hist[sample] <- sqrt(t(s)%*%solve(diag(length(s)))%*%s)
    procedure[[sample]]<-next_table
    #Relocate observation 
    if(sample<dim(testset)[1]){
      next_ob<-next_step$observed
      ob_hist[sample+1,] <- next_ob
    }
    
    #Check if exceed the control limits 
    if(L_hist[sample]>=L){
      return(list(RL = sample,ob_hist = ob_hist,procedure = procedure,L_hist = L_hist))
    }  
    current <- next_table
  }
  return(list(RL = 300,ob_hist = ob_hist,procedure = procedure,L_hist = L_hist))
}

#1.get histgram for control limit L to apply 1/200 rule
df_hist<-function(Entire,Init,obnum){
  L_hist <- rep(0,dim(Entire)[1]) #depend on train size
  s <- as.matrix(Init$table[,"W_i"])
  L_hist[1] <- sqrt(t(s)%*%solve(diag(length(s)))%*%s)
  current <- Init$table
  for(sample in 2:dim(Entire)[1]){
    next_step <- cusum_df(Entire,current,sample,obnum) 
    next_table<-next_step$table #return the latest table
    s <- as.matrix(next_table[,"W_i"])
    #compute in-control G(i)
    L_hist[sample] <- sqrt(t(s)%*%solve(diag(length(s)))%*%s)
    current <- next_table
  }
  #return hist of control limits 
  L_hist
}

#CUSUM procedure and sensor redistribution for FPCA
#input is the previous sample step's CUSUM table
cusum_df <- function(Entire,input,example,obnum){
  table <- input
  for (i in 1:p){
    if(table[i,]$rank<=obnum){
      #Compensation for the observed case
      scale <- E[i]
      ob_comp <- dpca_score(Entire[example,,i],dpca_list[[i]]$efns,dpca_list[[i]]$mua) #d*1
      table[i,"cusum_score"] <- sum(ob_comp)*scale #N(0,1) variable
    }
    else{
      #compensation for the unobserved case, which comes from individual N(0,1)
      ran_comp <- rnorm(1,0,1) #N(0,1) normal variable 
      table[i,"cusum_score"] <-ran_comp
    }
    #Local CUSUM
    table[i,"W_positive"]<-max(table[i,"W_positive"]+umin*table[i,"cusum_score"]-(umin^2)/2,0)
    table[i,"W_negative"]<-max(table[i,"W_negative"]-umin*table[i,"cusum_score"]-(umin^2)/2,0)
    table[i,"W_i"]<-max(table[i,"W_positive"],table[i,"W_negative"])
  }
  table$rank <- rank(-1*table$W_i,ties.method = "random")
  #allocate observable sensors for next sample step based on rank of W_i
  ob_sensor <- subset(table,rank<=obnum)$sensor_id
  list(table = table,observed = ob_sensor)
}

