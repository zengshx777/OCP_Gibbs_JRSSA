#Main Scrip to Perform Sampling
library(R2jags)

#Data Reading Script
source("Data_Reading.R")
#Loading the definition for model in JAGS
source("JAGS_MODEL.R")

#Define Logit function
logit<-function(x){1/(1+exp(-x))}
#for (res.id in c(1,2,9)){
res.id=4
#Reformat Data to be suitable for JAGS Modeling
  source("DataLoading.R")  
  for (i in 1:4){
    #Response Level:Number of Categories observed
    res.level=length(unique(data.jags[[i]][,response_index[res.id]]))
    
    #Formulate Data for jags input
    data.jags.input = list(
      y = data.jags[[i]][, response_index[res.id]] + res.level - 5,
      X = as.matrix(data.jags[[i]][,covariate_index]),
      t = data.jags[[i]][, treatment_index],
      nsample = nrow(data.jags[[i]]),
      N = length(unique(data.jags[[i]]$province)),
      RX =0.01*diag(length(covariate_index)),
      Z = data.jags[[i]]$ifppr,
      rd.ind = data.jags[[i]]$rd.ind,
      ncov = length(covariate_index),
      RT =0.01*diag(res.level - 1)
    )
    #Run Three Chains With different Initial Values
    if(res.level==5){
      stack.sim <-
        jags(
          data.jags.input,
          par = parameters,
          model = rr.model,
          n.iter = 150000,
          n.chains = 3,
          inits = list(
            list(
              "d" = -99 * (2 * data.jags.input$t - 1)-50,
              "c" = 1:4,
              "y1_obs" = data.jags.input$y,
              "y0_obs" = data.jags.input$y
            ),
            list(
              "d" = -100 * (2 * data.jags.input$t - 1)-50,
              "c" = 1:4 - 0.5,
              "y1_obs" = data.jags.input$y,
              "y0_obs" = data.jags.input$y
            ),
            list(
              "d" = -101 * (2 * data.jags.input$t - 1)-50,
              "c" = 0:3,
              "y1_obs" = data.jags.input$y,
              "y0_obs" = data.jags.input$y
            )
          )
        )
    }else{
      stack.sim <-
        jags(
          data.jags.input,
          par = parameters,
          model = rr.model.den,
          n.iter = 150000,
          n.chains = 3,
          inits = list(
            list(
              "d" = -99 * (2 * data.jags.input$t - 1)-50,
              "c" = 1:3,
              "y1_obs" = data.jags.input$y,
              "y0_obs" = data.jags.input$y
            ),
            list(
              "d" = -100 * (2 * data.jags.input$t - 1)-50,
              "c" = 1:3 - 0.5,
              "y1_obs" = data.jags.input$y,
              "y0_obs" = data.jags.input$y
            ),
            list(
              "d" = -101 * (2 * data.jags.input$t - 1)-50,
              "c" = 0:2,
              "y1_obs" = data.jags.input$y-1,
              "y0_obs" = data.jags.input$y-1
            )
          )
        )    
    }
    #Save MCMC Samples and JAGS Object
    stack.mcmc<-stack.sim$BUGSoutput$sims.matrix
    assign(paste(response_index[res.id],i,"MCMCSample",sep="_"),stack.mcmc)
    assign(paste(response_index[res.id],i,"JAGSObject",sep="_"),stack.sim)
    print(paste("==",res.id,i,"=="))
  }  
  save.image(file=paste(response_index[res.id],"result_monotonic_model.RData",sep="_"))
  
#}
