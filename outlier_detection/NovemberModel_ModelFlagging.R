
sample.data<-read.csv("outputs/sample.data.csv")
stations<-read.csv("../Desktop/stations.csv")
stations$LogWSA<-log10(stations$AREA_SQKM)
stations$Log_P_MEAN<-  log10(stations$P_MEAN + 0.0001)
stations$Log_N_MEAN<-  log10(stations$N_MEAN + 0.00001)

library(ggplot2)

predictors<-c("LogWSA","New_Long","New_Lat","SITE_ELEV","ELEV_RANGE",
                               "TEMP_00_09","PPT_00_09","SumAve_P","KFCT_AVE","BDH_AVE","MgO_Mean","Log_P_MEAN",
                               "CaO_Mean","PRMH_AVE","S_Mean","PCT_SEDIM","LPREM_mean","Log_N_MEAN")

  
  
stations.pred<-as.matrix(stations[,predictors])



Biotic.Group<-sort(unique(na.omit(stations$BioticGroup)))

sapply(1:length(stations$StationCode), function (site)
{
  sapply(1:length(Biotic.Group), function(BG)
  {
    RefBG<-as.matrix(stations[stations$BioticGroup==Biotic.Group[BG],predictors])
    mahalanobis(x=as.matrix(subset(stations[site,predictors])),center=colMeans(RefBG),cov(RefBG))
  })
}   )     
       





       BG<-1
       site<-1
       RefBG<-as.matrix(na.omit(stations[stations$BioticGroup==BG,predictors]))
       
       
       mahalanobis(x=as.matrix(subset(stations[1,predictors])),center=colMeans(RefBG),cov(RefBG))


mahal.dist<-mahalanobis(x=stations.pred,center=colMeans(stations.pred),cov(stations.pred))

stations$Mahal<-mahal.dist

ggplot(data=stations,aes(x=SiteSet,y=Mahal))+geom_boxplot()+
  scale_y_continuous(trans="log")


#####



Biotic.Group<-sort(unique(na.omit(stations$BioticGroup)))

junk<-sapply(split(stations[,c("StationCode",predictors)],stations$StationCode), function(site){

  sapply(1:length(Biotic.Group),function(BG){
#     print(BG)
    mahalanobis(x=as.matrix(site[1,predictors]),
              center=colMeans(stations.pred[which(stations$BioticGroup==BG),]), 
              cov=cov(stations.pred[which(stations$BioticGroup==BG),]),inverted=F)}) })

junk<-data.frame(junk,stations$BioticGroup,stations$SiteSet)

ggplot(data=na.omit(junk),aes(x=as.factor(stations.BioticGroup),y=junk))+
  geom_boxplot()+
  scale_y_continuous(trans="log")



#####
###########################;
#STEP 1 -- # Obtain the needed R objects;
# first, obtain the five R objects which define a single model (see model.build.v4.1.r);
#   "bugcal.pa"= Matrix of observed presence/absence (1/0) at calibration sites, for all taxa.
#
#   "grps.final" = Corresponding vector identifying the group ID's of calibration sites in bugcal.pa;
#
grps.final<-stations[which(stations$SiteSet=="RefCal"),"BioticGroup"]

#   "preds.final" = Vector with names of the chosen predictor variables, all of which must be available;
#
#   "grpmns" = Matrix of group (cluster) means at calibration sites;


datmat<-as.matrix(stations[which(stations$SiteSet=="RefCal"),predictors])
grpmns<-apply(datmat,2,function(x)tapply(x,grps.final,mean))
#
#   "covpinv" = Inverse of pooled covariance matrix among final predictor variables at calibration sites;

# 4.5) Calculate the inverse pooled covariance matrix(covpinv) of the preds.final variables at the calibration sites;
#First, calculate a list of covariance matrices for each group;
covlist<-lapply(split.data.frame(datmat,grps.final),cov);
#pooled cov matrix is weighted average of group matrices, weighted by group size. Johnson & Wichern, 11-64;
grpsiz<-table(grps.final)
ngrps<-length(grpsiz)
npreds<-length(predictors)
#zero out an initial matrix for pooled covariance;
covpool<-matrix(rep(0,npreds*npreds),nrow=npreds,dimnames=dimnames(covlist[[1]]));
#weighted sum of covariance matrices;
covlist.df <- lapply(1:11, function(i) (grpsiz[i]-1)*covlist[[i]])
covpool2 <- Reduce(`+`, covlist.df)
covpool2 <- covpool2/(sum(grpsiz) - ngrps)
covpinv <- solve(covpool2)
# for(i in 1:ngrps){covpool<-covpool+(grpsiz[i]-1)*covlist[[i]]};
# covpool<-covpool/(sum(grpsiz)-ngrps);#renormalize;
# covpinv<-solve(covpool); #inverse of pooled cov matrix;

## 3. -- predict the group (cluster) membership for all new sites. ;
# Follow RIVPACS assumption of weighting ;
# the membership probabilities by Calibration group size, as a prior probability;
# Also, flag any outlier sites, using the chi-squared statistic;
#Store probs in matrix, sites are rows, columns are groups;
#Uses mahalanobis function, where new vector is taken as the 'center', mu,;
#and matrix of means is taken as the 'data matrix', x;

# dmat<-as.matrix(prednew[,preds.final]); #matrix of predictor data for new samples, include only the predictor variables;
dmat<-stations.pred

#3.1 -- compute the critical chi-squared values for flagging outlier samples;
# df = the MINIMUM of (a)(number of groups-1), and (b) number of predictor variables;
# will flag each sample at P-value =.05 and also P-value =.01 level;
dff<-min(c(length(predictors),(length(Biotic.Group)-1)));
crit.01<-qchisq(0.99,df=dff);
crit.05<-qchisq(0.95,df=dff);





# 3.2 - construct empty matrix for predicted membership probabilities;


nsit.new<-dim(dmat)[[1]]; #number of new samples;
grpprobs<-matrix(rep(0,nsit.new*length(Biotic.Group)),nrow=nsit.new,
                 dimnames=list(dimnames(dmat)[[1]],dimnames(grpmns)[[1]]));

#Also construct empty data.frame for outlier flag;
# include minimum (squared)distance;
# Each sample is either a PASS (denote by 0) or FAIL (denote by 1) for the outlier test;
outlier.flag<-data.frame(outlier.05=rep(0,nsit.new),outlier.01=rep(0,nsit.new),dismin=rep(0,nsit.new),row.names=dimnames(dmat)[[1]]);

# 3.3 - Loop over ALL samples, compute vector of group membership probs for each sample and set its outlier flags;
#execute the following code piece as a single block;
##;
for(i in 1:nsit.new){;
                     #vector of squared Mahal. dist from current sample to each group mean;
                     dist<-mahalanobis(grpmns,dmat[i,],covpinv,inverted=T); #vector of distances;
                     grpprobs[i,]<-grpsiz*exp(-0.5*dist); # see Clarke et al. (2000);
                     grpprobs[i,]<-grpprobs[i,]/sum(grpprobs[i,]);
                     #check for outlier;
                     outlier.flag$dismin[i]<-min(dist); #save minimum distance;
                     if(outlier.flag$dismin[i]>crit.05)outlier.flag[i,'outlier.05']<-1;
                     if(outlier.flag$dismin[i]>crit.01)outlier.flag[i,'outlier.01']<-1;
} #finish sample loop;
#### sample membership probabilities complete for new samples;
               
outlier.flag$SiteSet<-stations$SiteSet

ggplot(data=outlier.flag, aes(x=SiteSet, y=log10(dismin)))+geom_boxplot()

table(as.factor(outlier.flag$outlier.01), outlier.flag$SiteSet)
               
               
               
               
               