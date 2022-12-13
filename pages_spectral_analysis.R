#pages spectral analysis
#November 2022

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
require(lemon)
require(ggforce)
#install.packages("PaleoSpec")
require(PaleoSpec)
#require(devtools)
#devtools::install_github("EarthSystemDiagnostics/paleospec")
###########
##this version uses the version of the PAGES database that I created which has had chronologies
##truncated below an expressed population signal of 0.85.
############
trw<-readRDS(paste("PAGES_crns_cleaned.rds",sep=''))
cru<-readRDS(paste("ExtractedMeanJJAtempsfromCRUData.rds", sep=''))
hadcru_raw<-readRDS(paste("ExtractedMeanJJAtempsfromHADCRUData_NOT_INFILLED.rds", sep=''))
hadcru<-readRDS(paste("ExtractedMeanJJAtempsfromHADCRUData_INFILLED.rds", sep=''))
berear<-readRDS(paste("ExtractedMeanJJAtempsfromBerkeleyEarth 2.rds", sep=''))


####
#get a list of unique site name and chronology metadatata
sites<-lapply(trw, function(x) x$dataSetName)
#lengths<-lapply(trw, function(x) length(x$paleoData_values)) #mean of 577 unfiltered, 479 filtered
#mean(unlist(lengths))
lat<-lapply(trw, function(x) x$geo_latitude)
lon<-lapply(trw, function(x) x$geo_longitude)
meta<-data.frame(cbind(unlist(sites), unlist(lat), unlist(lon)))
colnames(meta)<-c("sites","lat","lon")
#remove duplicates
meta<-distinct(meta)


#remove versions with variance stabilization, as these are not useful
trw<-trw[-which(sapply(trw, function(x) (x$paleoData_detrendingMethod == "AgeDependentStdCrnStb")))]
trw<-trw[-which(sapply(trw, function(x) (x$paleoData_detrendingMethod == "SsfCrnStb")))]

#option to analyze just TRW, MXD or both
#trw<-trw[-which(sapply(trw, function(x) (x$paleoData_proxy == "MXD")))]

#create indivudal data objects for each type for addtional exploration if needed
AgeDep<-trw[sapply(trw, function(x) all(x$paleoData_detrendingMethod == "AgeDependentStdCrn"))]
Ssf<-trw[sapply(trw, function(x) all(x$paleoData_detrendingMethod == "SsfCrn"))]
Negex<-trw[sapply(trw, function(x) all(x$paleoData_detrendingMethod == "NegEx"))]
RCS<-trw[sapply(trw, function(x) all(x$paleoData_detrendingMethod == "RCS"))]

#get the mean length of each type of proxy+detrending combiniation
#lengths<-lapply(AgeDep, function(x) length(x$paleoData_values))#474(both), #523(mxd), #331
#lengths<-lapply(Ssf, function(x) length(x$paleoData_values)) #486 (both), #536(mxd), #
#lengths<-lapply(Negex, function(x) length(x$paleoData_values)) #494(both) #541(mxd), #
#lengths<-lapply(RCS, function(x) length(x$paleoData_values)) #495 (both) #540(mxd), #


#####tests on smoothing level using modified Daniell fitler (spans command)
#test_ts<-AgeDep[1]
#spectrum <- spectrum(test_ts[[1]]$paleoData_values, spans=c(2,2), log="yes", plot=T)
####

#calculate spectra for all series in database using a light smoothing kernel 
all_spec_list<-list()

for(i in 1:length(trw)) {
  rwi <- trw[[i]]$paleoData_values
  name <- trw[[i]]$dataSetName
  proxy <- trw[[i]]$paleoData_proxy
  detrend<-trw[[i]]$paleoData_detrendingMethod
  spectrum <- spectrum(rwi, spans=c(2,2), log="no", plot=F)
  spectrum$freq<-spectrum$freq[-(1:2)]
  spectrum$spec<-spectrum$spec[-(1:2)]
  all_spec_list[[paste0(name, "_" , proxy, "_", detrend)]] <- spectrum
  all_spec_list[[i]]$detrend<-detrend
  all_spec_list[[i]]$proxy<-proxy
  all_spec_list[[i]]$name<-name
}

#########use Raphael's code to interpolate each spectrum using linear interpolation 
approx_specs<-list()
#test<-trw[[1]]

SpecApprox<-function(spec, xout=NULL,...){
  if(is.null(xout)){
    #frq.bnds<-log10(range(spec$freq))
    frq.bnds<-log10(c(1/1000,1/2)) #bin output by regular frequency bands
    xout<-10^seq(frq.bnds[1],frq.bnds[2], 0.01)
  }
  int.spec<-approx(x = spec$freq,y = spec$spec,xout = xout,...)
  rtrn.spec<-list(freq=int.spec$x,
                  spec=int.spec$y)
  class(rtrn.spec)<-"spec"
  return(rtrn.spec)
}

#test out different log smoothing options
#test_spec<-spectrum(test$paleoData_values, spans = c(2,2), log = "no", plot = T)
#test_smooth<-LogSmooth(test_spec, df.log = 0.05, removeFirst = 0, removeLast =  0)
#plot(test_spec)
#plot(test_smooth)
#test3<-SpecApprox(test_smooth,xout = 10^seq(log10(1/5000),log10(1/2),0.001))
#plot(test3)


#apply smoothing and/or linear interpolation to each spectra using SpecApprox
#option to smooth additionally using the LogSmooth command. I found this worked okay, but I still had a lot of 
#noise in the mean curve due to averaging so many spectra together. 
for(i in 1:length(all_spec_list)) {
  spc<-all_spec_list[[i]]
  #smooth<-LogSmooth(spc, df.log = 0.05, removeFirst = 0, removeLast =  0)
  spc.apprx<-SpecApprox(spc)
  name <- all_spec_list[[i]]$name
  proxy <- all_spec_list[[i]]$proxy
  detrend<-all_spec_list[[i]]$detrend
  approx_specs[[paste0(name, "_" , proxy, "_", detrend)]] <- spc.apprx
  approx_specs[[i]]$detrend<-detrend
  approx_specs[[i]]$name<-name
  approx_specs[[i]]$proxy<-proxy
}

##################process CRUTS and HADCRU/HADCRU raw data and BEST data
cru_years<-cru$Year
cru$Year<-NULL

hadcru_years<-hadcru$Year
hadcru$Year<-NULL

hadcru_raw_years<-hadcru_raw$Year
hadcru_raw$Year<-NULL

berear_years<-berear$Year
berear$Year<-NULL

#filter raw values from hadcru raw data to remove the empty pixels
#I've already interpolated across NA values up to 4 years in length when curating the data from HadCRU, so this is now removing all series prior to 
#that point, which are currently contained as series with NA values. 
#this is important because the spectra command will not operate if there are NA values present
hadcru_raw_filter<-list()

for(i in hadcru_raw) {
  filt<-if(all(is.na(i))) NA else((na.contiguous(i)))
  hadcru_raw_filter[[length(hadcru_raw_filter)+1]]<-unlist(filt)
}

#how many years are on average included in the analysis after filtering?
n<-lapply(hadcru_raw_filter, function(x) length(x))
n<-unlist(n)
mean(n)
#131 years on average

hadcru_filter<-list()

for(i in hadcru) {
  filt<-if(all(is.na(i))) NA else((na.contiguous(i)))
  hadcru_filter[[length(hadcru_filter)+1]]<-unlist(filt)
}

n<-lapply(hadcru_filter, function(x) length(x))
n<-unlist(n)
mean(n)
#172 years on average (good)          


#add names and combine climate datasets for efficiency
names<-names(cru)

cru_names<-paste0(names, "_cru")
hadcru_names<-paste0(names, "_hadcru")
hadcru_raw_names<-paste0(names, "_hadcru.raw")
bear_names<-paste0(names, "_bear")

names(cru)<-cru_names
names(hadcru_filter)<-hadcru_names
names(hadcru_raw_filter)<-hadcru_raw_names
names(berear)<-bear_names
 
#bidn all dataframes together
clim<-c(cru, hadcru_filter, hadcru_raw_filter, berear)


#perform the spectral analysis of CRU database of extracted grid cells,
#apply spec approx function to interpolate. 

#NOTE I'm clipping the lowest two frequency bands to where the spectra is unreliable due to series length
clim_specs<-list()

for(i in 1:length(clim)) {
  tmp<-clim[[i]]
  spectrum <- spectrum(tmp, spans=c(2,2), log="no", plot=F)
  spectrum$freq<-spectrum$freq[-(1:2)]
  spectrum$spec<-spectrum$spec[-(1:2)]
  #smooth<-LogSmooth(spectrum, df.log = 0.05, removeFirst = 0, removeLast =  0)
  spc.apprx<-SpecApprox(spectrum)
  name <- names(clim[i])
  clim_specs[[paste0(name)]] <- spc.apprx
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

###make data frame of CRU data
clim_df<-setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("freq","spec","data_type","proxy", "name"))

for(i in 1:length(clim_specs)) {
  freq <- clim_specs[[i]]$freq
  spec <- clim_specs[[i]]$spec
  name <- replicate(length(freq), names(clim_specs[i]))
  data_type<-replicate(length(freq), strsplit(name, '_')[[1]][2])
  proxy<-replicate(length(freq), "instrumental")
  spec_bind<-cbind(freq, spec, data_type,proxy, name)
  spec_bind<-as.data.frame(spec_bind)
  clim_df<-rbind(clim_df, spec_bind)
}



#for some reason the values output from the previous functions as character vectors. 
spec_df$freq<-as.numeric(spec_df$freq)
spec_df$spec<-as.numeric(spec_df$spec)

clim_df$freq<-as.numeric(clim_df$freq)
clim_df$spec<-as.numeric(clim_df$spec)


#do some filtering of the lower frequencies where the spectra become clearly dominated by a few long records
spec_df <- spec_df %>%
  filter(freq > 0.00666)

#clim_bear<-clim_df %>%
#  group_by(data_type)%>%
#  filter(data_type == "bear" & freq > 0.03)

# do some filtering of the hadcru raw (noisy at low frequencies) and cru (short series)
# this is a roundabout way, again there should be a dplyr filter option for this but i couldn't figure it out quickly
clim_hadcru_raw<-clim_df %>%
  group_by(data_type)%>%
  filter(data_type == "hadcru.raw")%>%
  filter(freq > 0.03)

clim_hadcru_raw<-as.data.frame(clim_hadcru_raw)
clim_hadcru_raw$freq<-as.numeric(clim_hadcru_raw$freq)
clim_hadcru_raw$spec<-as.numeric(clim_hadcru_raw$spec)

clim_cru<-clim_df %>%
  group_by(data_type)%>%
  filter(data_type == "cru")%>%
  filter(freq > 0.016)

clim_cru<-as.data.frame(clim_cru)
clim_crufreq<-as.numeric(clim_cru$freq)
clim_cru$spec<-as.numeric(clim_cru$spec)


clim_df <- clim_df%>%
  filter( data_type !="hadcru.raw" & data_type != "cru")

#add everything back together
clim_df<-rbind(clim_df, clim_hadcru_raw, clim_cru)


colnames(spec_df)[3]<-"data_type"
#make into one dataset for plotting
all_specs<-rbind(clim_df, spec_df)

all_specs <- all_specs%>%
  filter(!is.na(spec))

#calculate beta values on full dataframe by detrending type
betas <- all_specs %>% 
  group_by(data_type) %>%
  do(mod = lm(log(spec)~log(freq), data = .))%>%
  mutate(Slope = summary(mod)$coefficients[2]) %>%
  dplyr::select(-mod)              


#calculate a mean and variance for ribbon plot, with 95% confidence intervals
#note that this is NOT a good method of calculating CI intervals for these curves. 
spec_sums <- all_specs %>%
  group_by(data_type, freq) %>%
  summarise(mean.spec = mean(spec, na.rm = T),
            sd.spec = sd(spec, na.rm = T),
            var.spec = var(spec, na.rm = T),
            max.spec = max(spec, na.rm = T),
            min.spec = min(spec, na.rm = T),
            n.spec = n()) %>%
  mutate(se.spec = sd.spec / sqrt(n.spec),
         lower.ci.spec = mean.spec - qt(1 - (0.01 / 2), n.spec - 1) * se.spec,
         upper.ci.spec = mean.spec + qt(1 - (0.01 / 2), n.spec - 1) * se.spec)



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

#####normalize to a consistent starting value using the highest frequency bands of 2-8 years
#mean at high frequencies = 0.26
spec_high<-spec_sums %>%
  group_by(data_type)%>%
  filter(freq > 0.125)%>%
  summarise(mean_high = mean(mean.spec),
            ratio = mean_high/0.26,
            scaled = 1/ratio) 
  
  
#normalize each mean curve by it's respective scaling value
#there should be a better way to do this using a mutate function, but I couldn't figure it out quickly. 
#Also, this might allow you to modify the individual scaling values if you wanted to avoid overplotting, for example,
#dropping the tree-ring slightly to show difference in variability at high frequencies
spec_sums$mean.spec[spec_sums$data_type == "RCS"]<-spec_sums$mean.spec[spec_sums$data_type == "RCS"]*spec_high$scaled[spec_high$data_type=="RCS"]
spec_sums$mean.spec[spec_sums$data_type == "AgeDependentStdCrn"]<-spec_sums$mean.spec[spec_sums$data_type == "AgeDependentStdCrn"]*spec_high$scaled[spec_high$data_type=="AgeDependentStdCrn"]
spec_sums$mean.spec[spec_sums$data_type == "SsfCrn"]<-spec_sums$mean.spec[spec_sums$data_type == "SsfCrn"]*spec_high$scaled[spec_high$data_type=="SsfCrn"]
spec_sums$mean.spec[spec_sums$data_type == "NegEx"]<-spec_sums$mean.spec[spec_sums$data_type == "NegEx"]*spec_high$scaled[spec_high$data_type=="NegEx"]
spec_sums$mean.spec[spec_sums$data_type == "cru"]<-spec_sums$mean.spec[spec_sums$data_type == "cru"]*spec_high$scaled[spec_high$data_type=="cru"]
spec_sums$mean.spec[spec_sums$data_type == "hadcru"]<-spec_sums$mean.spec[spec_sums$data_type == "hadcru"]*spec_high$scaled[spec_high$data_type=="hadcru"]
spec_sums$mean.spec[spec_sums$data_type == "hadcru.raw"]<-spec_sums$mean.spec[spec_sums$data_type == "hadcru.raw"]*spec_high$scaled[spec_high$data_type=="hadcru.raw"]
spec_sums$mean.spec[spec_sums$data_type == "bear"]<-spec_sums$mean.spec[spec_sums$data_type == "bear"]*spec_high$scaled[spec_high$data_type=="bear"]

#I'm plotting using a loess smoothing function instead of the original line. This could also be done higher up in the workflow, but I was having trouble using specApprox and LogSmooth together
p1<-ggplot()+
  #geom_line(aes(x = freq, y = mean.spec, color = data_type), data = spec_sums)+
  theme_bw()+
  geom_smooth(data = spec_sums, aes(x=freq, y=mean.spec, color = data_type),span=0.1, method = "loess", se=FALSE)+
  #geom_ribbon(data = spec_sums, aes(x = freq, ymin = lower.ci.spec, ymax = upper.ci.spec, fill = data_type), alpha = 0.5)+
  theme(panel.border=element_blank(), axis.line=element_line())+
  #annotation_custom(tableGrob(betas), xmin=0.001, ymin=0.2)+
  #ylab(expression(PSD~(TRW^2~yr^-1)))+
  ylab(expression(PSD ^2~yr^-1))+
  xlab(expression(Delta~"t"~"(year)"))+
  scale_fill_manual(name = "Method", breaks = c("AgeDependentStdCrn", "NegEx", "RCS","SsfCrn","cru","hadcru", "hadcru.raw","bear"), 
                    labels = c("Age-dependent spline (n=523","Negative exponential curve (n=541)","Regional curve standardization (n=540)", "Signal-free detrending (n=536)", "CRU Summer Temp (C) (1901-2020)", "HadCRU Summer Temp Anomalies (C) (1850-2021)","HadCRU Summer Temp Anomalies (raw data) (C) (131 year average)", "Berkeley Earth Summer Temp Anomalies (1850-2021) (C)"), 
                    values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#e6ab02", "#a6761d", "#666666"))+
  scale_color_manual(name = "Method", breaks = c("AgeDependentStdCrn", "NegEx", "RCS","SsfCrn","cru","hadcru", "hadcru.raw","bear"), 
                    labels = c("Age-dependent spline (n=523)","Negative exponential curve (n=541)","Regional curve standardization (n=540)", "Signal-free detrending (n=536)", "CRU Summer Temp (C) (1901-2020)", "HadCRU Summer Temp Anomalies (C) (1850-2021)","HadCRU Summer Temp Anomalies (raw data) (C) (131 year average)", "Berkeley Earth Summer Temp Anomalies (1850-2021) (C)"), 
                    values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#e6ab02", "#a6761d", "#666666"))+
  scale_x_continuous(trans = trans_reverser('log10'), limits = c(0.5,0.003), 
                    breaks = c(0.3, 0.1, 0.033, 0.01, 0.0033,0.001), 
                    labels = c("3", "10","30","100","300","1000"))+
  #scale_x_continuous(trans=trans_reverser('log10'))+
  #scale_y_continuous(trans=c('log10'))+
  scale_y_continuous(trans=c('log10'), limits = c(0.1, 15)) +
  #scale_x_continuous(trans = trans_reverser('log10'), limits = c(0.5,0.001))+
  #theme1+
  coord_capped_cart(bottom="both", left="both")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("PAGES MXD and Climate temperture spectra")
p1

#ggsave(p1, file = "PAGES_MXD_and_clim_spectra_allcurves.png", width =12, height = 6, dpi = 300, units = "in")
