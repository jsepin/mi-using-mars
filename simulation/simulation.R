#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Simulation ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gc()
rm(list = ls())
if(sessionInfo()$platform=="x86_64-w64-mingw32/x64"){
  setwd("C:/Users/SepinJ/OneDrive - Universität Luzern/mi-using-mars/simulation")
}else{
  setwd("/home/jerome/OneDrive/phd/mi-using-mars/simulation")
}

require(tidyverse)
require(mice)
require(mvtnorm)
require(zeallot)
require(mixgb)
require(parallel)
source("mice.impute.mars.R")
source("functions_for_sim.R")

cores <- parallel::detectCores()
meta <- list("Y_mean_true" = 10, "cores" = cores)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Parameters ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

analysis_form <- "y ~ 1" 

methods_comp <- expand.grid(
  method = c("cc", "full","mars","mixgb","cart","rf","lasso.select.norm", "pmm")
  ,Nsim = 500
  ,stringsAsFactors = FALSE
) %>%
  arrange(method)

methods_comp$IDmethod <- 1:nrow(methods_comp)

exp_factor <- expand.grid(
   structure = c("constant", "linear","additive", "non-additive")
  ,missing =   c("constant", "linear","additive", "non-additive")
  ,error = "gaussian" #c("gaussian","binomial") 
  ,p = c(10,100)
  ,rho = 0
  ,n= 100
  ,m = 50
  ,stringsAsFactors = FALSE
) 
exp_factor2 <- expand.grid(
  structure = "non-additive"
  ,missing =  "non-additive"
  ,error = "gaussian" #c("gaussian","binomial") 
  ,p = 100
  ,rho = 0 
  ,n= c(200,400,800)
  ,m = 50
  ,stringsAsFactors = FALSE
)

exp_factor <- rbind(exp_factor, exp_factor2) %>% distinct()
exp_factor$ID = 1:nrow(exp_factor)

saveRDS(methods_comp, "results/methods_comp.RDS")
saveRDS(exp_factor, "results/exp_factor.RDS")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Run simulation ---- 
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time <- Sys.time()
set.seed(0815)
res <- lapply(1:nrow(exp_factor),function(e){ 
  # e = 1
  print(paste0("expfact: ",e, " of ",nrow(exp_factor)))
  # assign experimental factors
  eval(parse(
    text = paste0("c(",paste0(colnames(exp_factor),collapse = ", "),") %<-% exp_factor[e,]")
  ))
  # run methods on datasets
  res <- lapply(1:nrow(methods_comp),function(i){
    # i = 3
    # get methods
    eval(parse(
      text = paste0("c(",paste0(colnames(methods_comp),collapse = ", "),") %<-% methods_comp[",i,",]")
    ))
    print(paste(i," method: ",method))
    
    # execute simulation Nsim times for certain method
    # res_sim <- lapply(1:Nsim, function(k){
    res_sim <- mclapply(mc.cores = cores, 1:Nsim, function(k){
      # k = 14
      # print(k)
      set.seed(k)
      # simulate data
      data <- dgp_and_miss(n,structure,missing,p,rho,transform=TRUE,error = error)
      
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      # apply method (mice, complete case, mixgb,...)
      #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      time <- Sys.time()
      if(method == "cc"){
        mice_res <- NA
        # print("using cc method")
        mod <- glm(data = data$amp, formula(analysis_form),family = error)
        if(!is.na(coef(mod)["(Intercept)"])){
          mice_res <- summary(mod)$coef["(Intercept)",1:2]
          mice_res <- data.frame(t(mice_res),mod$df.residual)
          names(mice_res) <- c("estimate","std.error","df")
        }else{
          mice_res <- data.frame("estimate"=NA,"std.error"=NA,"df"=NA)
        }
        
      }else if(method == "full"){
        mice_res <- NA
        # print("using full method")
        mod <- glm(data = data$data, formula(analysis_form),family = error)
        if(!is.na(coef(mod)["(Intercept)"])){
          mice_res <- summary(mod)$coef["(Intercept)",1:2]
          mice_res <- data.frame(t(mice_res),mod$df.residual)
          names(mice_res) <- c("estimate","std.error","df")
        }else{
          mice_res <- data.frame("estimate"=NA,"std.error"=NA,"df"=NA)
        }
        
      }else if(method == "mixgb"){
        mice_res <- NA
        xgb.params <- list(
          nthread =1, max_depth = 3, gamma = 0, eta = 0.3, min_child_weight = 1
          , subsample = 0.7 # 0.7 default; important otherwise undercoverage!
          , colsample_bytree = 1, colsample_bylevel = 1, colsample_bynode = 1,
          tree_method = "auto"
        )
        mixgb <- mixgb(data = data$amp, m = m, maxit = 1, xgb.params = xgb.params, nrounds = 100, bootstrap = FALSE )
        mixgb.fit <- lapply(mixgb, function(dataset) glm(formula = analysis_form, data = dataset,family = error))
        imp_info <- mixgb.fit %>% pool()
        mice_res <- imp_info$pooled %>% 
          dplyr::select(-c(term,m))%>%
          mutate(std.error = sqrt(ubar + b + b/m), .after = "estimate")
        
        }else{
        mice_res <- NA
        mice_m <- tryCatch(
          {
            mice(data = data$amp,method = method,m = m, maxit = 1,printFlag = FALSE,eps=0) # sampling = sampling,oob = oob
          },
          error = function(cond) {"error"}
        )
        if(!any(mice_m=="error")){ 
          if(sum(is.na(mice_m %>% mice::complete("long", include = FALSE)))!=0){
            mice_res <- data.frame("estimate"=NA,"std.error"=NA, "df" = NA)
          }else{
            imp_info <- mice_m %>% with(data = ., glm(formula = formula(analysis_form),family = error )) %>% pool()
            
            mice_res <- imp_info$pooled %>% 
              dplyr::select(-c(term,m))%>%
              mutate(std.error = sqrt(ubar + b + b/m), .after = "estimate")
          }
        }else{
          mice_res <- data.frame("estimate"=NA,"std.error"=NA, "df" = NA)
        }
      }
      
      time <- as.numeric(difftime(Sys.time(),time, units = 'secs'))
      names(mice_res) <- paste0("mice_", names(mice_res))
      
      # combine
      res <- mice_res %>%
        mutate(
           y_mean = data$y_mean
          ,n_true = data$n_true
          ,n_cc = data$n_cc
          ,perc_missing_actual = (n_true-n_cc)/n_true
          ,time = time
          ,k = k
        )
      
      return(res)
    })
    res_sim <- do.call("bind_rows", res_sim)
    res_sim$IDmethod <- IDmethod
    
    return(res_sim)
  })
  res <- do.call("bind_rows", res) %>%
    left_join(., methods_comp, by = "IDmethod") %>%
    mutate(ID = ID)
  return(res)
})
res <- do.call("bind_rows", res) %>%#do.call("rbind", res) %>%
  left_join(., exp_factor, by = "ID")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Save results ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

(time <- Sys.time()-time)
meta$time <- time
saveRDS(res, file = "results/sim_fin.RDS")
saveRDS(meta, file = "results/meta.RDS")

