L_train<-function(Train,obnum){
  L<-rep(NA,length = 5)
  #Estimating control limit h for the five methods 
  mf <-ARL_initial(Train,obnum)
  mf_glb <- mf_hist(Train,mf,Corr,obnum)
  L[1] <- quantile(mf_glb,199/200)  
  
  fp <-ARLi_initial(Train,obnum)
  fp_glb <- df_hist(Train,fp,obnum)
  L[2] <- quantile(fp_glb,199/200) 
  
  d001 <-ARLi_init_delta(Train,obnum,del = 0.01, q = 4)
  d001_glb <- df_hist_delta(Train,d001, obnum,del = 0.01, q = 4)
  L[3] <- quantile(d001_glb$L_hist,199/200)
  
  d05 <-ARLi_init_delta(Train, obnum,del = 0.05, q = 4)
  d05_glb <- df_hist_delta(Train,d05, obnum,del = 0.05, q = 4)
  L[4] <- quantile(d05_glb$L_hist,199/200)
  
  d10 <-ARLi_init_delta(Train, obnum,del = 0.1, q = 4)
  d10_glb <- df_hist_delta(Train,d10, obnum,del = 0.1, q = 4)
  L[5] <- quantile(d10_glb$L_hist,199/200) 
  L
}
