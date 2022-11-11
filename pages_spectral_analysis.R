#pages spectral analysis
#July 14, 2020

#put PAGES_crs_cleaned.rds into folder with script
require(dplyr)
require(tidyr)
require(ggplot2)
require(extrafont) 
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(scales)
require(graphics)
require(utils)
require(dplR)

trw<-readRDS(paste("PAGES_crns_cleaned.rds",sep=''))

#get a list of unique site names
sites<-lapply(trw, function(x) x$dataSetName)
lengths<-lapply(trw, function(x) length(x$paleoData_values))

#remove versions with variance stabilization, as these are not useful
trw<-trw[-which(sapply(trw, function(x) (x$paleoData_detrendingMethod == "AgeDependentStdCrnStb")))]
trw<-trw[-which(sapply(trw, function(x) (x$paleoData_detrendingMethod == "SsfCrnStb")))]

#AgeDep<-trw[sapply(trw, function(x) all(x$paleoData_detrendingMethod == "AgeDependentStdCrn"))]
#Ssf<-trw[sapply(trw, function(x) all(x$paleoData_detrendingMethod == "SsfCrn"))]
#Negex<-trw[sapply(trw, function(x) all(x$paleoData_detrendingMethod == "NegEx"))]
#RCS<-trw[sapply(trw, function(x) all(x$paleoData_detrendingMethod == "RCS"))]

#####tests on smoothing level using modified Daniell fitler (spans command)
#test_ts<-AgeDep[1]
#spectrum <- spectrum(test_ts[[1]]$paleoData_values, spans=c(2,2), log="yes", plot=T)
####

#calculate spectra for all series represnted
all_spec_list<-list()

for(i in 1:length(trw)) {
  rwi <- trw[[i]]$paleoData_values
  name <- trw[[i]]$dataSetName
  proxy <- trw[[i]]$paleoData_proxy
  detrend<-trw[[i]]$paleoData_detrendingMethod
  spectrum <- spectrum(rwi, spans=c(2,2), log="no", plot=F)
  all_spec_list[[paste0(name, "_" , proxy, "_", detrend)]] <- spectrum
  all_spec_list[[i]]$detrend<-detrend
  all_spec_list[[i]]$proxy<-proxy
  all_spec_list[[i]]$name<-name
}

#########use Raphael's code to smooth each spectrum using linear interpolation 
approx_specs<-list()
test<-trw[[1]]


SpecApprox<-function(spec, xout=NULL,...){
  if(is.null(xout)){
    frq.bnds<-log10(range(spec$freq))
    xout<-10^seq(frq.bnds[1],frq.bnds[2],0.05)
  }
  int.spec<-approx(x = spec$freq,y = spec$spec,xout = xout,...)
  rtrn.spec<-list(freq=int.spec$x,
                  spec=int.spec$y)
  class(rtrn.spec)<-"spec"
  return(rtrn.spec)
}

test_spec<-spectrum(test$paleoData_values, spans = c(2,2), log = "no", plot = T)
test_spec<-SpecApprox(test_spec)

plot(test_spec)

for(i in 1:length(all_spec_list)) {
  spc<-all_spec_list[[i]]
  spc.apprx<-SpecApprox(spc)
  name <- all_spec_list[[i]]$name
  proxy <- all_spec_list[[i]]$proxy
  detrend<-all_spec_list[[i]]$detrend
  approx_specs[[paste0(name, "_" , proxy, "_", detrend)]] <- spc.apprx
  approx_specs[[i]]$detrend<-detrend
  approx_specs[[i]]$name<-name
  approx_specs[[i]]$proxy<-proxy
}

#take series out of list and transform into one long dataframe for plotting and stats. 
#Takes a minute

spec_df<-setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("freq","spec","detrending_type","proxy","name"))

for(i in 1:length(all_spec_list)) {
  freq <- approx_specs[[i]]$freq
  spec <- approx_specs[[i]]$spec
  detrending_type<-approx_specs[[i]]$detrend
  proxy<-approx_specs[[i]]$proxy
  name <- replicate(length(freq), names(approx_specs[i]))
  spec_bind<-cbind(freq, spec, detrending_type, proxy, name)
  spec_bind<-as.data.frame(spec_bind)
  spec_df<-rbind(spec_df, spec_bind)
}


#preserve original 
spec_df_original<-spec_df
spec_df<-spec_df_original

#######option to bin by frequency####
#spec_df$freq[spec_df$freq <= 0.000] <- 0.0005
#spec_df$freq[spec_df$freq <= 0.002 & spec_df$freq > 0.001]<-"0.002" #500-1000
#spec_df$freq[spec_df$freq <= 0.004 & spec_df$freq > 0.002]<-"0.004" #250-500
#spec_df$freq[spec_df$freq <= 0.01 & spec_df$freq > 0.004]<-"0.01" #100-250
#spec_df$freq[spec_df$freq <= 0.02 & spec_df$freq > 0.01]<-"0.02" #50-100
#spec_df$freq[spec_df$freq <= 0.05 & spec_df$freq > 0.02]<-"0.05" #20-50
#spec_df$freq[spec_df$freq <= 0.1 & spec_df$freq > 0.05]<-"0.1" #10-20
#spec_df$freq[spec_df$freq > 0.1 ]<-"0.1" #0-10
#unique(spec_df$freq)
#######
#########

#for some reason the values output from the previous functions as character vectors. 
spec_df$freq<-as.numeric(spec_df$freq)
spec_df$spec<-as.numeric(spec_df$spec)

#calculate beta values on full dataframe by detrending type
spec_betas <- spec_df %>%
  group_by(detrending_type) %>%
  do(mod = lm(log(spec)~log(freq), data = .))%>%
  mutate(Slope = summary(mod)$coefficients[2]) %>%
  dplyr::select(-mod)              


#do a bit of rounding, but be careful because the lowest frequencies can get rounded to zero. 
spec_df$freq<-round(spec_df$freq, 3)
spec_df$spec<-round(spec_df$spec, 4)

#if needed you can make everything less than 0.001 equal to 0.001, there arent' that many series that go so low. 
spec_df$freq[spec_df$freq <= 0.001] <- 0.001
spec_df$freq<-format(spec_df$freq, scientific=F)

#make frequency into a factor to calculate confidence intervals. 
spec_df$freq<-as.factor(spec_df$freq)

spec_sums <- spec_df %>%
  group_by(detrending_type, freq) %>%
  summarise(mean.spec = mean(spec, na.rm = T),
            sd.spec = sd(spec, na.rm = T),
            var.spec = var(spec, na.rm = T),
            max.spec = max(spec, na.rm = T),
            min.spec = min(spec, na.rm = T),
            n.spec = n()) %>%
  mutate(se.spec = sd.spec / sqrt(n.spec),
         lower.ci.spec = mean.spec - qt(1 - (0.01 / 2), n.spec - 1) * se.spec,
         upper.ci.spec = mean.spec + qt(1 - (0.01 / 2), n.spec - 1) * se.spec)


spec_sums$freq<-as.numeric(as.character(spec_sums$freq))

spec_betas <- spec_sums %>%
  group_by(detrending_type) %>%
  do(mod = lm(log(mean.spec)~log(freq), data = .))%>%
  mutate(Slope = summary(mod)$coefficients[2]) %>%
  dplyr::select(-mod)  


spec_sums<-spec_sums %>%
  filter(freq > 0.0001)

#d a little smoothing on the mean and confidence bands to reduce noise at high frequencies
spec_smooths <- spec_sums %>%
  group_by(detrending_type) %>%
  summarise(smoothed.mean = ffcsaps(mean.spec, nyrs = 10, f=0.5),
            smoothed.upper.ci = ffcsaps(upper.ci.spec, nyrs = 10, f=0.5),
            smoothed.lower.ci = ffcsaps(lower.ci.spec, nyrs = 10, f=0.5))

spec_smooths$freq<-spec_sums$freq
spec_smooths$smoothed.upper.ci<-round(spec_smooths$smoothed.upper.ci, 5)
spec_smooths$smoothed.lower.ci<-round(spec_smooths$smoothed.lower.ci, 5)

spec_smooths$smoothed.lower.ci[spec_smooths$smoothed.lower.ci < 0] <-0.00001

theme1<-theme(legend.key = element_blank(),
                text = element_text(size = 12),
                title = element_text(size = 12),
                strip.text.x = element_text(size = 12),
                strip.background.x = element_blank(),
                strip.text.y = element_text(size = 12),
                #panel.border=element_blank(), 
                #axis.line=element_line())
                #axis.line = element_blank(),
                legend.position=c(.3,.75))

require(lemon)

spec_smooths<-as.data.frame(spec_smooths)

spec_sums<-spec_sums%>%
  filter(lower.ci.spec != "NaN")

p1<-ggplot()+
  #geom_line()+ 
  geom_ribbon(data = spec_sums, aes(x = freq, ymin = lower.ci.spec, ymax = upper.ci.spec, fill = detrending_type), 
              colour = NA, alpha = 0.4)+
  geom_line(aes(x = freq, y = mean.spec, color=detrending_type), data = spec_sums)+
  theme_bw()+
  theme(panel.border=element_blank(), axis.line=element_line())+
  ylab(expression(PSD~(K^2~yr^-1)))+
  xlab(expression(Delta~"t"~"(year)"))+
  scale_fill_manual(name = "Method", breaks = c("AgeDependentStdCrn", "NegEx", "RCS","SsfCrn"), 
                    labels = c("Age-dependent spline","Negative exponential curve","Regional curve standardization", "Signal-free detrending"), 
                    values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"))+
  scale_color_manual(name = "Method", breaks = c("AgeDependentStdCrn", "NegEx", "RCS","SsfCrn"), 
                    labels = c("Age-dependent spline","Negative exponential curve","Regional curve standardization", "Signal-free detrending"), 
                    values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"))+
  scale_y_continuous(trans=c('log10'), limits = c(0.005, 10)) +
  scale_x_continuous(trans=c("log10" , "reverse"), limits = c(0.47,0.001), breaks = c(0.33, 0.1, 0.033, 0.01, 0.0033,0.001), 
                     labels = c("3", "10","30","100","300","1000"))+
  theme1+
  coord_capped_cart(bottom="both", left="both")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("PAGES spectra")
p1

ggsave(p1, file = "PAGES.png", width = 5.1, height = 5, dpi = 300, units = "in")
