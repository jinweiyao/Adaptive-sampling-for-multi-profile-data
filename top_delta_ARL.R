#FPCA + Del = 0.01,0.06,0.1 approach under CUSUM framework

ARLi_init_delta<-function(testset,obnum,delta,q){
  p<-dim(testset)[3]
  cusum_score<-rep(0,length=p)
  #update two sided local CUSUM
  table <- data.frame(sensor_id = seq(1,p,length.out = p),cusum_score)
  for (i in 1:obnum){
    ob_init<-dpca_score(testset[1,,i],dpca_list[[i]]$efns,dpca_list[[i]]$mua) 
    table[i,"cusum_score"] <- sum(ob_init)*E[i]
    table[i,"W_positive"]<-max(umin*table[i,"cusum_score"]-(umin^2)/2,0)
    table[i,"W_negative"]<-max(-umin*table[i,"cusum_score"]-(umin^2)/2,0)
    table[i,"W_i"]<-max(table[i,"W_positive"],table[i,"W_negative"])
  }
  for(i in (obnum+1):p){
    #unobserved compensation become delta = 0.01,0.05,0.1
    table[i,"cusum_score"] <-NA
    table[i,"W_positive"]<- delta
    table[i,"W_negative"]<- delta
    table[i,"W_i"]<-max(table[i,"W_positive"],table[i,"W_negative"])
  }
  table["rank"] = rank(-table$W_i,ties.method = "random")
  ob_sensor <- subset(table,rank<= obnum)$sensor_id
  top_q<-subset(table,rank<= q)$sensor_id
  list(table = table ,observed = ob_sensor,top_q = top_q)
}


Get_ARLi_delta<-function(testset,L,obnum,delta,q){
  procedure <-vector(mode ="list",length = 300)
  ob_hist<-array(0,dim = c(300,obnum))
  L_hist <- rep(0,dim(testset)[1])
  init_list <- ARLi_init_delta(testset,obnum,delta,q)
  ob_hist[1,] <- seq(1,obnum,length.out = obnum)
  ob_hist[2,] <- init_list$observed
  procedure[[1]]<-init_list$table
  s <- sum(init_list$table[,"W_i"][init_list$top_q]) #top_q global statistics
  L_hist[1] <- s
  if(L_hist[1]>=L){
    return (list(RL = 1,ob_hist = ob_hist,procedure = procedure,L_hist = L_hist))
  }
  current <- init_list
  for(sample in 2:dim(testset)[1]){
    next_step <- cusum_df_delta(testset,current,sample,obnum,delta,q) 
    next_table<- next_step$table  
    s <- sum(next_table[,"W_i"][next_step$next_top]) #top_q global statistics
    L_hist[sample] <- s
    procedure[[sample]]<-next_table
    if(sample<dim(testset)[1]){
      ob_hist[sample+1,]<-next_step$observed
    }
    if(L_hist[sample]>=L){
      return(list(RL = sample,ob_hist = ob_hist,procedure = procedure,L_hist = L_hist))
    }  
    current <- next_step
  }
  return(list(RL = 300,ob_hist = ob_hist,procedure = procedure,L_hist = L_hist))
}

df_hist_delta<-function(Entire,Init,obnum,delta,q){
  procedure<-vector(mode = "list",length = dim(Entire)[1])
  procedure[[1]]<-Init$table
  L_hist <- rep(0,dim(Entire)[1])
  ob_hist<-array(0,dim = c(dim(Entire)[1],obnum))
  s<-sum(Init$table[,"W_i"][Init$top_q]) #top_q global statistics
  L_hist[1] <- s
  ob_hist[1,] <- seq(1,obnum,length.out = obnum)
  ob_hist[2,] <- Init$observed
  current <- Init
  for(sample in 2:dim(Entire)[1]){
    next_step <- cusum_df_delta(Entire,current,sample,obnum,delta,q) 
    next_table<-next_step$table 
    procedure[[sample]]<-next_table
    s <- sum(next_table[,"W_i"][next_step$next_top]) #top_q global statistics
    L_hist[sample] <- s
    if(sample<dim(Entire)[1]){
      ob_hist[sample+1,]<-next_step$observed
    }
    current <-next_step
  }
  return(list(L_hist = L_hist,ob_hist = ob_hist, procedure = procedure))
}

#input is a list,comprive of next_top, next_ob and table 
cusum_df_delta <- function(Entire,input,example,obnum,delta,q){
  table<-input$table
  p<-dim(Entire)[3]
  for (i in 1:p){
    if(table[i,]$rank<=obnum){ #determine ob or not
      ob_comp <- dpca_score(Entire[example,,i],dpca_list[[i]]$efns,dpca_list[[i]]$mua) #d*1
      table[i,"cusum_score"] <- sum(ob_comp)*E[i] #N(0,1) variable
      table[i,"W_positive"]<-max(table[i,"W_positive"]+umin*table[i,"cusum_score"]-(umin^2)/2,0)
      table[i,"W_negative"]<-max(table[i,"W_negative"]-umin*table[i,"cusum_score"]-(umin^2)/2,0)
      table[i,"W_i"]<-max(table[i,"W_positive"],table[i,"W_negative"])
    }
    else{
      #unobserved compensation become delta = 0.01,0.05,0.1 from FPCA
      table[i,"cusum_score"] <-NA
      table[i,"W_positive"]<-table[i,"W_positive"] + delta
      table[i,"W_negative"]<-table[i,"W_negative"] + delta
      table[i,"W_i"]<-max(table[i,"W_positive"],table[i,"W_negative"])
    }
  }
  table$rank = rank(-1*table$W_i,ties.method = "random")
  ob_sensor <- subset(table,rank<=obnum)$sensor_id
  next_top <- subset(table,rank<=q)$sensor_id
  list(table = table,observed = ob_sensor,next_top = next_top)
}


