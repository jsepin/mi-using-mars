#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Visualization ----
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
require(ggpubr)
require(scales) 

res <- readRDS("results/sim_fin.RDS")
meta <- readRDS("results/meta.RDS")
methods_comp <- readRDS("results/methods_comp.RDS")
exp_factor <- readRDS("results/exp_factor.RDS")

res <- res %>%
    mutate(true_estimate = meta$Y_mean_true ) %>%
    mutate(structure = factor(structure, levels = c("constant", "linear", "additive", "non-additive"))) %>%
    mutate(missing = factor(missing, levels = c("constant", "linear", "additive", "non-additive"))) %>%
    mutate(facet = paste0("missing: ",missing,"\np: ",p,"\nn: ",n, "\nm: ",m))

res_full <- res %>%
  filter(method =="full") %>%
  dplyr::select(ID, k,error,mice_estimate, mice_std.error) %>%
  rename(full_estimate = mice_estimate
         ,full_std.error = mice_std.error
         )
res <- res %>% 
  left_join(.,res_full,by = c("ID","k","error")) %>%
  mutate(bias = mice_estimate - full_estimate)# full_estimate ----

# method order
method_order <- c("full", "cc", "mars")
last_method  <- c("pmm","lasso.select.norm")
method_order <- c(method_order, unique(res$method)[!(unique(res$method) %in% c(method_order,last_method) )],last_method)
res$method <- factor(res$method, method_order)

# facet order
facet_order <- lapply(c(" constant", " linear", " additive", " non-additive"),function(i){
    unique(res$facet)[str_detect(unique(res$facet), i)]
}) %>% unlist()
res$facet <- factor(res$facet, facet_order)


# are there NA values? 
n_na <- res[rowSums(is.na(res))>0,] %>% dplyr::filter(!(method %in% c("cc", "full"))) %>% nrow()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## function to save plot as png and tiff ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# saver function
save_stuff <- function(p,width, height){
  name <- deparse(substitute(p))
  name <- str_replace(name, "\\$","_")
  # saving png
  png(paste0("../little_sim/figures/",name,".png")
      ,width=width, height=height, units = "in"
      ,res = 100
  )
  print(p)
  dev.off()
  
  # saving tiff
  tiff(paste0("../little_sim/figures/",name,".tiff")
       ,width=width, height=height, units = "in"
       ,res = 800
       ,compression = "lzw"
  )
  print(p)
  dev.off()
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Bias & Standard error & Coverage plot----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_function <- function(data){
  # description label
  p_desc <- data %>%
    dplyr::select(bias,mice_std.error,method,facet,structure,true_estimate) %>%
    group_by(method,facet,structure) %>%
    summarise(mean_bias = mean( bias)
              ,sd_bias  = sd( bias) 
              ,mean_se = mean( mice_std.error)
              ,sd_se = sd( mice_std.error)
              ,upper75qu_bias = quantile(bias, 0.75)
              ,upper75qu_se = quantile(mice_std.error, 0.75)
    ) %>%
    group_by(structure) %>%
    mutate(upper75qu_bias = max(upper75qu_bias)
           ,upper75qu_se = max(upper75qu_se))

p_desc$label_bias <- paste0(
    paste0(format(round(p_desc$mean_bias,2), nsmall = 2)," ")
    ,paste0("(",format(round(p_desc$sd_bias,2), nsmall = 2),")")
)
p_desc$label_se <- paste0(
    paste0(format(round(p_desc$mean_se,2), nsmall = 2)," ")
    ,paste0("(",format(round(p_desc$sd_se,2), nsmall = 2),")")
)

# bias plot
p_bias <- data %>%
    ggplot(aes( y = method
                ,x = bias
                ,fill = method
    ) )+
    labs(y = "", x = ""
         ,title = bquote("Bias: "~hat(mu)[s]-hat(mu)[s]^"(full)" )
         ,subtitle = bquote(m==.(unique(data$m))~","~Nsim==.(unique(data$Nsim))~"(outliers and whiskers removed for illustration)")
         ,fill ="")+
    geom_label(data = p_desc,aes(x = upper75qu_bias, y = method, label = label_bias),
               hjust = "outward", vjust = 0.5,fill = "white",label.size = NA,size = 3)+
    geom_boxplot(outliers = F,coef = 0)+
    theme_bw()+
    geom_vline(xintercept = 0,lty = 3)+
    facet_grid( facet ~ structure
                ,scales = "free_x"
                ,labeller = labeller(structure = label_both)
    )+
    theme(
        legend.position = "bottom"
        ,panel.spacing = ggplot2::unit(0, "lines")
        ,strip.text.y = element_text(size = 8, colour = "black", angle = 0)
    )+
    guides(fill = guide_legend(nrow = 1))+
    labs(x = NULL)+
    scale_x_continuous(expand = expansion(mult=c(0,0.8)))

# standard error
p_se <- data %>%
  ggplot(aes( y = method
              ,x = mice_std.error
              ,fill = method
  ) )+
  labs(y = "",title = bquote("Standard error: "~SE(hat(mu))[s]), x = ""
       ,subtitle = bquote(m==.(unique(data$m))~","~Nsim==.(unique(data$Nsim))~"(outliers and whiskers removed for illustration)")
       ,fill ="")+
  geom_label(data = p_desc,aes(x = upper75qu_se, y = method, label = label_se),
             hjust = "outward", vjust = 0.5,fill = "white",label.size = NA,size = 3)+
  geom_boxplot(outliers = F,coef = 0)+
  theme_bw()+
  facet_grid( facet ~ structure
              ,scales = "free_x"
              ,labeller = labeller(structure = label_both)
  )+
  theme(
    legend.position = "bottom"
    ,panel.spacing = ggplot2::unit(0, "lines")
    ,strip.text.y = element_text(size = 8, colour = "black", angle = 0)
  )+
  guides(fill = guide_legend(nrow = 1))+
  labs(x = NULL)+
  scale_x_continuous(expand = expansion(mult=c(0,0.6)))

# Coverage
mycov <- function(x){
  out <- prop.test(x = sum(x,na.rm=TRUE), n = length(x))
  out <- data.frame(out$estimate, t(out$conf.int), row.names = NULL)
  names(out) <- c("cov_estimate", "cov_lower", "cov_upper")
  return(out)
}

res_coverage <- data %>%
  mutate(di = abs(bias)/mice_std.error < qnorm(0.975)) %>%
  group_by(facet,method,structure) %>%
  group_modify(~ mycov(.x$di))
res_coverage$label_cov <- paste0(
  paste0(format(round(res_coverage$cov_estimate,2), nsmall = 2)," ")
  ,paste0("[",format(round(res_coverage$cov_lower,2), nsmall = 2)
          ,"-",format(round(res_coverage$cov_upper,2), nsmall = 2)
          ,"]")
)

base_size = 11
ggtext_size <- function(base_size, ratio = 0.8) {
  ratio * base_size / ggplot2::.pt
}

p_coverage <- res_coverage %>%
  ggplot(aes( y = method
              ,x = cov_estimate
              ,colour = method
  ) )+
  labs(y = ""
       ,x = ""
       ,title = bquote("Coverage: "~abs(hat(mu)[s]-hat(mu)[s]^"(full)")/SE(hat(mu))[s]<z[0.975])
       ,subtitle = bquote(m==.(unique(data$m))~","~Nsim==.(unique(data$Nsim))~", Coverage probability with 95% Wilson confidence interval")
       ,col = "")+
  geom_point()+
  geom_errorbar(aes(xmin = cov_lower, xmax = cov_upper),show.legend = F)+
  geom_label(data = res_coverage,aes(x = 1.1, y = method, label = label_cov),
             hjust = 0, vjust = 0.5,col= "black",fill = "white",label.size = NA,size = ggtext_size(base_size) )+
  theme_bw(base_size = base_size)+
  coord_cartesian(clip = "off")+
  geom_vline(xintercept = 1-0.05)+
  facet_grid( facet ~ structure
              # ,scales = "fixed"
              # ,ncol = 4
              ,labeller = labeller(structure = label_both))+
  theme(
    legend.position = "bottom"
    ,panel.spacing = ggplot2::unit(0, "lines")
    ,strip.text.y = element_text(size = 8, colour = "black", angle = 0)
  )+
  guides(col = guide_legend(nrow = 1))+
  scale_x_continuous(breaks = c(seq(0,1,length.out = 5)[-5],0.95)
                     ,labels= c(seq(0,1,length.out = 5)[-5],0.95)
                     ,expand = expansion(mult=c(0,0.95))
  )+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5))+
  labs(x = NULL)
  return(list("bias" = p_bias, "se" = p_se, "coverage" = p_coverage))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Execute ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p <- plot_function(data = res %>% filter(n==min(res$n)))
# save primary simulation plots
save_stuff(p=p$bias,width=10, height=12)
save_stuff(p=p$se,width=10, height=12)
save_stuff(p=p$coverage,width=10, height=12)

p_secondary <- plot_function(data = res %>%
                                filter(missing == "non-additive" &
                                         structure == "non-additive" &
                                         p == 100 
                                         ))

# combine and save secondary simulation plots
p_secondary_combined <- ggarrange(p_secondary$bias+theme(strip.text.y = element_blank())+labs(subtitle = bquote(m==.(unique(res$m))~","~Nsim==.(unique(res$Nsim))))
          , p_secondary$se+theme(strip.text.y = element_blank(), axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+labs(subtitle = "(outliers and whiskers removed for illustration)")+scale_x_continuous(expand = expansion(mult=c(0,0.3)))
          , p_secondary$coverage+theme( axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+labs(subtitle = "Coverage probability with 95% confidence interval")
, ncol=3, nrow=1, common.legend = TRUE, legend="bottom",align = "h")

p_secondary_combined
save_stuff(p=p_secondary_combined,width=10, height=8)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Time ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# description label
p_desc <- res %>%
  filter(!(method %in% c("full","cc"))) %>%
  dplyr::select(time,method,facet,n,p,m) %>%
  group_by(method,n,p) %>%
  summarise(mean_time = mean(time)
            ,sd_time  = sd( time)
            ,upper75qu_time = quantile(time, 0.75)
  )
p_desc$upper75qu_time <- max(p_desc$upper75qu_time)

p_desc$label_time <- paste0(
  paste0(format(round(p_desc$mean_time,2), nsmall = 2)," ")
  ,paste0("(",format(round(p_desc$sd_time,2), nsmall = 2),")")
)

# same colors 
cols <- hue_pal()(length(unique(res$method)))

time_factors <- res %>% dplyr::select(p,n,m) %>% distinct() %>% arrange(p,n,m) %>%
  mutate(facet = paste0("p: ",p, "\nn: ",n, "\nm: ",m) ) %>%
  dplyr::select(facet) %>% pull

# time plot
p_time <- res %>%
  filter(!(method %in% c("full","cc"))) %>%
  mutate(facet = factor(paste0("p: ",p, "\nn: ",n, "\nm: ",m) , levels = time_factors)
         )%>%
  ggplot(aes( y= method
                ,x = abs(time)
                ,fill = method
    ) )+
  #geom_boxplot()+
  geom_boxplot(outliers = F,coef = 0)+
    labs(y = "", x = "Time in seconds", title = "Imputation time"
         ,subtitle = bquote(m==.(unique(res$m))~","~Nsim==.(unique(res$Nsim))~"(outliers and whiskers removed for illustration)")
         ,fill = "")+
  geom_label(data = p_desc,aes(x = upper75qu_time+1, y = method, label = label_time),
             hjust = "outward", vjust = 0.5,fill = "white",label.size = NA,size = 3)+
    theme_bw()+
    facet_grid(  facet~.
                ,scales = "free"
                ,labeller = labeller(structure = label_both)
                )+
    theme(
        legend.position = "bottom"
        ,panel.spacing = ggplot2::unit(0, "lines")
        ,strip.text.y = element_text(size = 8, colour = "black", angle = 0)
    )+
    guides(fill = guide_legend(nrow = 1))+
    scale_x_continuous(expand = expansion(mult=c(0,0.15)))+
    scale_fill_manual(values = cols[-c(1:2)] )

p_time
save_stuff(p=p_time,width=8, height=8) 

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Variance components ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var_comp <- function(res){
  
# same colors 
cols <- hue_pal()(length(unique(res$method)))

# description label
p_desc <- res %>%
  dplyr::select(mice_ubar,mice_b,method,facet,structure,n,p
                , full_estimate, full_std.error) %>%
  filter(!method %in% c("cc","full")) %>%
  group_by(method,facet,structure,n,p) %>%
  summarise(mean_ubar = mean(mice_ubar)
            ,sd_ubar  = sd(mice_ubar)
            ,upper75qu_ubar = quantile(mice_ubar, 0.75)
            ,mean_b = mean(mice_b)
            ,sd_b  = sd(mice_b)
            ,upper75qu_b = quantile(mice_b, 0.75)
  )
p_desc$upper75qu_ubar <- max(p_desc$upper75qu_ubar)
p_desc$upper75qu_b <- max(p_desc$upper75qu_b)


p_desc$label_ubar <- paste0(
  paste0(format(round(p_desc$mean_ubar,2), nsmall = 2)," ")
  ,paste0("(",format(round(p_desc$sd_ubar,2), nsmall = 2),")")
)
p_desc$label_b <- paste0(
  paste0(format(round(p_desc$mean_b,2), nsmall = 2)," ")
  ,paste0("(",format(round(p_desc$sd_b,2), nsmall = 2),")")
)

# within variance plot (mice_ubar/mice_b)
p_withinvar <- res %>%
  filter(!method %in% c("cc","full")) %>%
  ggplot(aes( y = method
              ,x = mice_ubar
              ,fill = method
  ) )+
  labs(y = "",title = bquote("Within-variance: "~bar(U)==frac(1,m)~sum(bar(U)[l], l==1, m) ), x = ""
       ,subtitle = bquote(m==.(unique(res$m))~","~Nsim==.(unique(res$Nsim))~"(outliers and whiskers removed for illustration)")
       ,fill ="")+
  geom_boxplot(outliers = F,coef = 0)+
  #geom_boxplot()+
  # geom_label(data = p_desc,aes(x = upper75qu_ubar, y = method, label = label_ubar),
  #            hjust = "outward", vjust = 0.5,fill = "white",label.size = NA,size = 3)+
  theme_bw()+
  facet_grid( facet ~ structure
              ,scales = "free_x"
              ,labeller = labeller(structure = label_both)
              #,labeller = label_both
  )+
  theme(
    legend.position = "bottom"
    ,panel.spacing = ggplot2::unit(0, "lines")
    ,strip.text.y = element_text(size = 8, colour = "black", angle = 0)
  )+
  guides(fill = guide_legend(nrow = 1))+
  scale_fill_manual(values = cols[-c(1:2)] )+
  labs(x = NULL)

p_betweenvar <- res %>%
  filter(!method %in% c("cc","full")) %>%
  ggplot(aes( y = method
              ,x = mice_b
              ,fill = method
  ) )+
  labs(y = "",title = bquote("Between-variance: "~B==Var(hat(mu)[s])~"(impact can to some extent be reduced by increasing m)")
       ,x = ""
       ,subtitle = bquote(m==.(unique(res$m))~","~Nsim==.(unique(res$Nsim))~"(outliers and whiskers removed for illustration)")
       ,fill ="")+
  geom_boxplot(outliers = F,coef = 0)+
  #geom_boxplot()+
  # geom_label(data = p_desc,aes(x = upper75qu_b, y = method, label = label_b),
  #           hjust = "outward", vjust = 0.5,fill = "white",label.size = NA,size = 3)+
  theme_bw()+
  facet_grid( facet ~ structure
              ,scales = "free_x"
              ,labeller = labeller(structure = label_both)
  )+
  theme(
    legend.position = "bottom"
    ,panel.spacing = ggplot2::unit(0, "lines")
    ,strip.text.y = element_text(size = 8, colour = "black", angle = 0)
  )+
  guides(fill = guide_legend(nrow = 1))+
  scale_fill_manual(values = cols[-c(1:2)] )+
  labs(x = NULL)

return(list("within" = p_withinvar, "between"=p_betweenvar))

}

# dev.off()
# ggarrange(p_withinvar, p_betweenvar)

p_varcomp <- var_comp(res %>% filter(n==100))
save_stuff(p=p_varcomp$within,width=10, height=8)
save_stuff(p=p_varcomp$between,width=10, height=8)

p_varcomp_secondary <- var_comp(res %>%
                        filter(missing == "non-additive" &
                                 structure == "non-additive" &
                                 p == 100 
                        ))
save_stuff(p=p_varcomp_secondary$within,width=10, height=8)
save_stuff(p=p_varcomp_secondary$between,width=10, height=8)


