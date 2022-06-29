rm(list = ls()) # delete parameters in the workspace
library(scales)

mean_dist <- 50
nsd_dist <- 15
sd_dist <- 40
distances <- seq(0,100,1)
nsamples <- 1000
xlim <- 100
hist.resolution <- 1
alpha.col <- 0.75

# sample from single value
sdisp_vals <- runif(nsamples, min=0, max=xlim) #
scsum <- cumsum(xtabs(~round(sdisp_vals)))
#plot(scsum)
sprobs_disp <- distances/max(distances)
sdprobs_disp <- rep(1,length(distances))/length(distances)
smedian <- median(sdisp_vals)

# sample from normal distribution
ndisp_vals <- rnorm(nsamples,mean=mean_dist,sd = nsd_dist) #
nprobs_disp <- pnorm(distances,mean=mean_dist,sd = nsd_dist)
ndprobs_disp <- dnorm(distances,mean=mean_dist,sd = nsd_dist)
#ndprobs_disp <- ndprobs_disp/max(ndprobs_disp)
nmedian <- median(ndisp_vals)

# hn <- hist(ndisp_vals,xlim = c(0,xlim),n=max(distances)/hist.resolution/(sd_dist/nsd_dist),
#            col = alpha("lightgray",alpha.col),border=alpha("lightgray",alpha.col),probability = TRUE)
# nmax_freq <- max(hn$counts)
# lines(distances,ndprobs_disp*nmax_freq,lwd=2)
# abline(v=mean_dist,col='blue',lty=1,lwd=2)
# abline(v=nmedian,col='red',lty=2,lwd=2)


# sample from negative binomial
cv_dist = sd_dist/mean_dist # coefficient of variation
p_dist = 1/cv_dist/sd_dist # probability of sampling success to represent SD
r_dist = mean_dist * p_dist / (1-p_dist); # change of mean to sample in order to maintain actual mean 
probs_disp <- pnbinom(distances,size=r_dist,prob = p_dist) # calculate probabilities of larval dispersal
dprobs_disp <- dnbinom(distances,size=r_dist,prob = p_dist)
#dprobs_disp <- dprobs_disp/max(dprobs_disp)
drops_disp_rev <- rev(dprobs_disp)
#quant_disp <- pnbinom(distances,size=r_dist,prob = p_dist) # calculate probabilities of larval dispersal

disp_vals <- rnbinom(nsamples,size=r_dist,prob = p_dist) # calculate probabilities of larval dispersal
median_dist <- median(disp_vals)
disp_vals_rev <- abs(disp_vals-xlim)
dmax <- max(disp_vals)
mode_disp_vals <- which.max(tabulate(match(disp_vals, unique(disp_vals))))
#max_freq <- max(tabulate(match(disp_vals, unique(disp_vals))))

# hnb <- hist(disp_vals,xlim = c(0,xlim),n=max(distances)/hist.resolution,
#             col = alpha("lightgray",alpha.col),border=alpha("lightgray",alpha.col))
# max_freq <- max(hnb$counts)
# lines(distances,dprobs_disp*max_freq,lwd=2)
# abline(v=mean_dist,col='blue',lty=1,lwd=2)
# abline(v=median,col='red',lty=2,lwd=2)
# 
# hist(disp_vals_rev,xlim = c(0,xlim),n=max(distances)/hist.resolution,
#      col = alpha("lightgray",alpha.col),border=alpha("lightgray",alpha.col))
# lines(rev(distances),dprobs_disp*max_freq,lwd=2)
# abline(v=mean_dist,col='blue',lty=1,lwd=2)
# abline(v=xlim-median,col='red',lty=2,lwd=2)


# Protection plots relative to MPA size
# ncsum <- cumsum(xtabs(~round(sdisp_vals)))
# ncsum.dists <- as.numeric(names(ncsum))
# #nbar <- barplot(ncsum/nsamples,ncsum.dists,col = alpha("lightgray",alpha.col),border=alpha("lightgray",alpha.col),xlim=c(0,xlim))
# plot(ncsum.dists,ncsum/nsamples,xlim=c(0,xlim),type='l',lwd=2)

ncsum <- cumsum(xtabs(~round(ndisp_vals)))
ncsum.dists <- as.numeric(names(ncsum))
#nbar <- barplot(ncsum/nsamples,ncsum.dists,col = alpha("lightgray",alpha.col),border=alpha("lightgray",alpha.col),xlim=c(0,xlim))
# plot(distances,nprobs_disp,xlim=c(0,xlim),type='l',lwd=2,
#      ylab='Individuals protected (%)',xlab="MPA size (km)")
# lines(distances,seq(0,1,1/xlim),lty=1,lwd=2)
# abline(v=mean_dist,col='blue',lty=1,lwd=2)
# abline(v=mean_dist,col='red',lty=2,lwd=2)

csum <- cumsum(xtabs(~round(disp_vals)))
csum.dists <- as.numeric(names(csum))
# plot(csum.dists,csum/nsamples,type='l',lwd=2,xlim=c(0,xlim),
#      ylab='Individuals protected (%)',xlab="MPA size (km)")
# 
# plot(distances,probs_disp,type='l',lwd=2,xlim=c(0,xlim),
#      ylab='Individuals protected (%)',xlab="MPA size (km)")
# #lines(distances,seq(0,1,1/xlim),lty=3)
# abline(v=mean_dist,col='blue',lty=1,lwd=2)
# abline(v=median,col='red',lty=2,lwd=2)


# plot everything in one figure ####

par(mfrow=c(4,3))
#> par("mar")
#[1] 5.1 4.1 4.1 2.1
par(mar=c(4.1,4.1,1.1,2.1)) # c(bottom, left, top, right)


hist.metrics <- rbind(sdisp_vals,ndisp_vals,disp_vals,disp_vals_rev)/10
fit.metrics <- rbind(sdprobs_disp,ndprobs_disp,dprobs_disp,rev(dprobs_disp))*10
median.metrics <- c(mean_dist,mean_dist,median_dist,xlim-median_dist)
protection.metrics <- round(rbind(sprobs_disp,nprobs_disp,probs_disp,rev(1-probs_disp))*100)
difference.metrics <- round(rbind(rep(0,length(distances)),nprobs_disp-sprobs_disp,probs_disp-sprobs_disp,rev(1-probs_disp)-sprobs_disp)*100)

xlim.plot<-10
nrows <- 4
xlabs <- rep("",nrows) #,nrow=nrows,ncol=ncols)
xlabs[nrows] <- "MPA size (km)"
xlabs1 <- xlabs
xlabs1[nrows] <- "Movement distance (km)"

count <- 0
  
#for(yplot in 1:ncols){
for(xplot in 1:nrows){
    count <- count + 1
    # plot histograms (xplot == 1)
    
     #if(xplot==1){
      hnb <- hist(hist.metrics[xplot,],xlim = c(0,xlim.plot),n=max(distances)/hist.resolution,
      col = alpha("lightgray",alpha.col),border=alpha("lightgray",alpha.col),freq = FALSE,
      ylab = 'Relative Frequency',xlab=xlabs1[xplot],main=NULL)
      lines(distances/10,fit.metrics[xplot,],lwd=2)
      abline(v=mean_dist/10,col='blue',lty=1,lwd=2)
      abline(v=median.metrics[xplot]/10,col='red',lty=2,lwd=2)
      
      #} else if (xplot==2){
      plin <- plot(distances/10,protection.metrics[xplot,],type='l',lwd=2,xlim=c(0,xlim.plot),
              ylab = 'Survival (% individuals)',xlab=xlabs[xplot],main=NULL)
      abline(v=median.metrics[xplot]/10,col='red',lty=2,lwd=2)
      
      #} else if (xplot==3){
      if(xplot==1){
        plot.new()
        legend("topleft",legend=c("Mean","Median"),col=c("blue","red"),lty=c(1,2),lwd=2,bty="n")
      } else {
      dlin <- plot(distances/10,difference.metrics[xplot,],type='l',lwd=2,xlim=c(0,xlim.plot),
              ylab = 'Difference in survival',xlab=xlabs[xplot],main=NULL)
      } 
      

   #  } else { # level 1 condition
   # 
   #    #if(xplot==1){
   #    hnb <- hist(hist.metrics[xplot,],xlim = c(0,xlim),n=max(distances)/hist.resolution,
   #    col = alpha("lightgray",alpha.col),border=alpha("lightgray",alpha.col),freq = FALSE,
   #    ylab = 'Relative Frequency',xlab=NULL,main=NULL)
   #    lines(distances,fit.metrics[xplot,],lwd=2)
   #    abline(v=mean_dist,col='blue',lty=1,lwd=2)
   #    abline(v=median.metrics[xplot],col='red',lty=2,lwd=2)
   #    
   #    #} else if (xplot==2){
   #    plin <- plot(distances,protection.metrics[xplot,],type='l',lwd=2,xlim=c(0,xlim),
   #    ylab='Individual survival (%)',xlab=NULL)
   #    abline(v=mean_dist,col='blue',lty=1,lwd=2)
   #    abline(v=median.metrics[xplot],col='red',lty=2,lwd=2)
   #    abline(v=median.metrics[xplot],col='red',lty=2,lwd=2)
   #    
   #    #} else if (xplot==3){
   #    dlin <- plot(distances,difference.metrics[xplot,],type='l',lwd=2,xlim=c(0,xlim),
   #    ylab='Difference in survival (%)',xlab=NULL)
   #    #} # level 2 condition
   #    
   #          #} # level 2 condition
   # } # level 1 condition
} # xplot loop
#} # yplot loop


# save hist metrics to file
hist.metrics.save<-hist.metrics
hist.metrics.save<-hist.metrics.save[1:3,]
saveRDS(hist.metrics.save,'hist.metrics.RDS')
