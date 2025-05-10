#X1 will be a single profile sample readings, so dim(X1) = 1*n*p 
#eval,efns,mua are from estimated mfpca model 
m_score<-function(X1,efns,mua){
  demeaned = X1-mua[1,,] #n*p
  n <- dim(efns)[1]
  d <- dim(efns)[2]
  p <- dim(mua)[3]
  scores <-array(0,dim = c(d,p))
  #see E.q 10
  for(l in 1:d){
    scores[l,] = colSums(demeaned*efns[,l])/n 
  }
  scores #d*p matrix
}

#Input X is of m*n*p matrix (number of samples, number of time points, number of sensors)
##d0 is the number of decomposition dimension. in E.q 5
mfpca_est<-function(X,d0){
  m = dim(X)[1]
  n = dim(X)[2]
  p = dim(X)[3]
  mua <- array(0,dim = c(1,n,p))
  # Equation 6 in the paper 
  for(i in 1:p){
    for(j in 1:n){
      mua[1,j,i] = mean(X[,j,i]) 
    }
  }
  #Convarinace function in E.q 7 in the paper 
  cova = array(0,dim = c(n,n))
  for(u in 1:m){
    sum_p = array(0,dim = c(n,n))
    for (i in 1:p){
      sum_p = sum_p + (X[u,,i]-mua[1,,i])%*%t(X[u,,i]-mua[1,,i]) 
    }
    cova = cova + sum_p
  }
  cova = cova/m
  #Eigen Decompostion is in E.q 8  
  eio<-eigen(cova)
  eval<-eio$values/n
  evec<-eio$vectors*sqrt(n)
  #d0 = 6
  d = d0
  eval0 <- eval[1:d]
  efns0 <- evec[,1:d]
  #scores that is calculated by E.q 10
  scores <-array(0,dim = c(m,d,p))
  for(i in 1:m){
    scores[i,,]<-m_score(X[i,,],efns0,mua)
  }
  #Estimation covariance matrix of MFPCA score Ksi, introduced in E.q 9
  covb<-array(0,dim=c(d,p,p))
  for(l in 1:d){
    for (j in 1:p){
      for (h in 1:p){
        for(i in 1:m){
          a = (X[i,,j]-mua[1,,j])*efns0[,l]/n
          b = (X[i,,h]-mua[1,,h])*efns0[,l]/n
          sum_n = sum(a)*sum(b)
          covb[l,j,h] = covb[l,j,h]+sum_n
        }
      }
    } 
  }
  covb = covb/m
  
  #output lsit includes: 
  #1.dimension, 2.template profile, 3.eigenvalue, 4.eigenfucntions, 5.mfpc-scores and 6.covariacne matrix 
  list(d=d0,mua=mua,evalue = eval0,efns = efns0,scores = scores,covb = covb)
}


#mainly for proposed method
ARL_initial<-function(testset,obnum){ 
  r<-obnum
  #set up initial table and setting up initial compensation
  cusum_score<-rep(0,length=p) 
  data_samp = testset[1,,]
  ob_init <- m_score(data_samp,mf_est$efns,mf_est$mua) #d*p matrix score table
  for(i in 1:r){
    scale <- E1[i,i] #in E.q 7
    cusum_score[i] <- sum(ob_init[,i])*scale #N(0,1) variable
  }
  ran_init <- mvrnorm (mu=rep (0, times = p),Sigma = Corr)
  for(i in (r+1):p){
    cusum_score[i] <- ran_init[i] #N(0,1) variable
  }
  #introduced in E.q 1&2
  table <- data.frame(sensor_id = seq(1,p,length.out = p),cusum_score)
  for(i in 1:length(cusum_score)){
    table[i,"W_positive"]<-max(umin*cusum_score[i]-(umin^2)/2,0)
    table[i,"W_negative"]<-max(-umin*cusum_score[i]-(umin^2)/2,0)
    table[i,"W_i"]<-max(table[i,"W_positive"],table[i,"W_negative"])
  }
  #order local statistics in E.q 3
  table["rank"] = rank(-table[,"W_i"],ties.method = "random")
  ob_sensor <- subset(table,rank<= r)$sensor_id
  list(table = table,observed = ob_sensor)
}

Get_ARL <- function(testset,Corr,L,obnum){
  procedure <-vector(mode ="list",length = 300)
  ob_hist<-array(0,dim = c(300,obnum))
  L_hist <- rep(0,dim(testset)[1]) #depend on train size
  Init <- ARL_initial(testset,obnum)
  ob_hist[1,] <- seq(1,obnum,length.out = obnum)
  ob_hist[2,]<-Init$observed
  procedure[[1]]<-Init$table
  s <- as.matrix(Init$table[,"W_i"]) #initiate an CUSUM chart
  L_hist[1] <- sqrt(t(s)%*%solve(Corr)%*%s) #compute G(i) using E.q 9
  if(L_hist[1]>=L){ #check if exceed control limits h in E.q 10
    return (list(RL = 1,observed = ob_hist,procedure = procedure,L_trace = L_hist))
  }
  current <- Init$table
  for(sample in 2:dim(testset)[1]){
    next_step <- cusum_stage(testset,current,sample,obnum)
    next_table<- next_step$table
    
    s <- as.matrix(next_table[,"W_i"])
    L_hist[sample] <- sqrt(t(s)%*%solve(Corr)%*%s) #compute G(i) using E.q 9
    procedure[[sample]]<-next_table
    if(sample<dim(testset)[1]){
      next_ob<-next_step$observed 
      ob_hist[sample+1,] <- next_ob
    }
    if(L_hist[sample]>=L){ #Check if G(i) exceed control limits L in E.q 10
      return(list(RL = sample,ob_hist = ob_hist,procedure = procedure,L_hist = L_hist))
    }
    current <- next_table
  }
  return(list(RL = 300,ob_hist = ob_hist,procedure = procedure,L_hist = L_hist))
}

#get histgram for mfpca method L 
mf_hist<-function(Entire,Init,Corr,obnum){
  L_hist <- rep(0,dim(Entire)[1]) #depend on train size
  s <- as.matrix(Init$table[,"W_i"])
  #in E.q 9
  L_hist[1] <- sqrt(t(s)%*%solve(Corr)%*%s)
  current <- Init$table
  for(sample in 2:dim(Entire)[1]){
    next_step <- cusum_stage(Entire,current,sample,obnum) #return the latest table
    s <- as.matrix(next_step$table[,"W_i"])
    #in E.q 9 to estimate in-control G(i) histgram
    L_hist[sample] <- sqrt(t(s)%*%solve(Corr)%*%s)
    current <- next_step$table
  }
  L_hist
}

#mainly for proposed method 
cusum_stage <- function(Entire,input,example,obnum){
  table <- input
  ob_comp <- m_score(Entire[example,,],mf_est$efns,mf_est$mua) #return d*p
  #random compensation 
  ran_comp = mvrnorm (mu=rep (0, times = p),Sigma = Corr)
  for (i in 1:p){
    if(table[i,]$rank<=obnum){
      #observed sensor and aggregate all the d = 1,2,...d0 for the mfpc scores
      scale = E1[i,i] #scaling in E.q 7
      table[i,"cusum_score"] = sum(ob_comp[,i])*scale #marginal N(0,1) normal variable 
    }
    else{
      #unobserved sensor sample from multivariate normal
      table[i,"cusum_score"] = ran_comp[i] #marginal N(0,1) normal variable 
    }
    #CUSUM procedure in E.q 1&2
    table[i,"W_positive"]<-max(table[i,"W_positive"]+umin*table[i,"cusum_score"]-(umin^2)/2,0)
    table[i,"W_negative"]<-max(table[i,"W_negative"]-umin*table[i,"cusum_score"]-(umin^2)/2,0)
    table[i,"W_i"]<-max(table[i,"W_positive"],table[i,"W_negative"])
  }
  #order local statistics in E.q 3
  table$rank <- rank(-1*table$W_i,ties.method = "random")
  ob_sensor <- subset(table,rank<=obnum)$sensor_id
  list(table = table,observed = ob_sensor)
}