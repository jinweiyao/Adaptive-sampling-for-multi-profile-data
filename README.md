# Adaptive_sampling_for_multi-profile

To replicate the simulation result, one just need to run 'Int_simulation.R" file. We mainly used MASS and Splines packages, also use source function to access the other R file. 

Interpret the result: 

The output of the simulation will be based on partial shift pattern. The output from the 'Model1_simulation.R' correspond to fig.2(a) and 'Model2_simulation.R' correspond to fig.2(b) respectively. Putting them together in 'Int_simulation.R' to complete the fig.2 from the manuscript. The output will have 6 subplots and 3  each column correspond to available resource r = 8,6,4. Each subplot has five lines correspond to the proposed, fpca, fpca with Δ = 0.01,0.05,0.1 methods. These line plots have mu-shift of 0.7,1.0,1.4,2.1,2.8 on the x-axis and corresponded ARL result on the y-axis. A table of ARL_1 is also generated in variable out_se, with ARL_1 information in every odd line and the standard error in every even row (correspond to table A1 of the appendix). 

To modify the simulation setting, one should go to "Model1_simulation.R" and "Model2_simulation.R". The shift type can be changed from partial shift into full shift by changing line 75-85 In the demo code, we used only for 100 repletion for computing ARL for demonstrating purpose, the estimated running time is about an hour depends on the devise you use. (in the actual simulation, we did 500 times). If one want to increase the number of trial, it's easy to edit line 64-70, line 134-138 by changing both i and length of from oc_run1 trough oc_run5. 



Detailed of files layout and functions contained

Demo Folders have 6 other files : 


1.top_delta_ARL.R: include methods required for ARL_1 for fpca Δ = 0.01,0.1,0.5
For fpca with delta: ARLi_init_delta:helper function for initiate the CUSUM chart 
                     Get_ARLi_delta:Get running length
                     df_hist_delta:Estimate control limit L
                     cusum_df_delta:helper function for redistribution
2. fpca_ARL.R: include methods required for ARL_1 for fpca chart
For fpca: dpca_score:calculate FPC scores from a profile 
	  dpca_est:estimate FPCA model and parameters
	  ARLi_initial:helper function for initiate CUSUM chart 
          Get_ARLi:Get running length
          df_hist:Estimate control limit L
          cusum_df:helper function for redistribution
3.mfpca_ARL.R: include methods required for ARL_1 for proposed method 
For proposed: mfpca_score:calculate MFPC scores from a multi-profile samples
	      mfpca_est:estimate MFPCA model and parameters
	      ARL_initial: initiate CUSUM chart
              Get_ARL: Get running length
              mf_hist: Estimate control limit L
              cusum_stage: helper function for redistribution
4.control limits train.R 
	L_train: tune the control limits for all methods

5.Model1_simulation.R
	Run the simulation for Model I Partial shift 

6.Model2_simulation.R
	Run the simulation for Model II Partial shift



Link for data used in case study: https://www.kaggle.com/datasets/podsyp/production-quality?select=data_X.csv


