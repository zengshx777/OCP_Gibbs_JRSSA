rm(list=ls())
args=commandArgs(trailingOnly = TRUE)
if(length(args)==0){
  print("No arguments supplied")
}else{
  for (i in 1:length(args)){
    eval(parse(text=args[i]))
  }
}
source("Data_Reading.R")
###Outcome Index
###1 for confidence
###2 for Anxiety
###4 for Desperation
res.id=1
###Subgroup Index
###1 for rural female
###2 for rural male
###3 for urban female
###4 for urban male
sub_group=1
source("Data_Loading.R")
source("Data_Preprocess.R")
###Sensitity Parameter
###delta.z=0 under assumptions.
###
mcmc.time=50000
t_approx=0
#delta.z=0
source("Gibbs_Sampler.R")
#results_sensitivity=NULL
results_sensitivity=rbind(results_sensitivity,
                          c(delta.z,mean(att),
                            quantile(att,0.025),
                            quantile(att,0.975),
                            mean(prte),
                            quantile(prte,0.025),
                            quantile(prte,0.975)))