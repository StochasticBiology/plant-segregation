### plot data summaries and posteriors for the test bed datasets
# sparsely commented here -- the code is the same as in parse-data.R and plot-posteriors-data.R (with fewer plots generated)

library(ggplot2)
library(gridExtra)

# transform effective cell division counts, given an old and a new population size
# simulations are run with n=50. we might want, for example, n=200 for mitos or n=10 for plastids
mytrans = function(n, ne.orig=50, ne.new=100) {
  return(n*log(1-1/ne.orig)/log(1-1/ne.new))
}

# initialise set of plots
subplots0 = subplots1 = subplots2 = list()

mtplot.ne = 50
ptplot.ne = 7

gsf = 3
ressf = 4.5

defaultplot.max = 40
mtplot.max = 150

# set of vectors encoding different experimental setups -- used for experiment-specific tweaks to plotting style
# vector of filenames
expts = c("data-bed-1", "data-bed-2", "data-bed-3", "data-bed-4", "data-bed-5")

var.plot = data.plot = list()

gsf = 3
ressf = 4.5
for(eref in 1:length(expts)) {
  
  expt = expts[eref]
  filename = paste(c(expt, ".csv"), collapse="")
  
  df = read.csv(filename, header=T)
  colnames(df) = c("id", "family", "time", "h")
  
  df$h = df$h/100
  
  lines = data.frame()
  families = unique(df$family)
  for(family in families) {
    sub = df[df$family == family,]
    ref = sub$h[which(sub$time == 0)]
     for(offspring in unique(sub$id)) {
      meanh = mean(sub$h[sub$id == offspring & sub$time == 1])
      if(!is.na(meanh)) {
        lines = rbind(lines, data.frame(family=family,id=offspring,h0=ref,h=meanh))
      }
    }
  }
  
  stats = data.frame()
  for(family in families) {
    sub = lines[lines$family == family,]
    hs = sub$h
    hdiffs = (hs-mean(hs))**2
    vest = sum(hdiffs)/(length(hs)-1)
    #   stats = rbind(stats, data.frame(family=family, stage=0, vhp = var(hs)/(mean(hs)*(1-mean(hs)))))
    stats = rbind(stats, data.frame(family=family, stage=0, vhp = vest/(mean(hs)*(1-mean(hs)))))
  }
  
  points = data.frame()
  for(id in 1:max(df$id)) {
    sub = df[df$id == id,]
    if(nrow(sub) > 1) {
      points = rbind(points, data.frame(family=sub$family, id=id, time=sub$time+1, h=sub$h))
    }
  }
  
  olines = data.frame()
  for(id in 1:max(df$id)) {
    sub = df[df$id == id,]
    if(nrow(sub) > 1) {
      ref = mean(sub$h[sub$time==1])
      olines = rbind(olines, data.frame(family=sub$family, h0=ref, time=sub$time, h=sub$h))
      subsub = sub[sub$time==1,]
      vest = sum((subsub$h-ref)**2)/(length(subsub$h)-1)
      if(nrow(subsub) > 1) {stats = rbind(stats, data.frame(family=sub$family[1], stage=1, vhp = vest/(ref*(1-ref)))) }
      subsub = sub[sub$time==2,]
      vest = sum((subsub$h-ref)**2)/(length(subsub$h)-1)
      if(nrow(subsub) > 1) {stats = rbind(stats, data.frame(family=sub$family[1], stage=2, vhp = vest/(ref*(1-ref)))) }
      subsub = sub[sub$time==3,]
      vest = sum((subsub$h-ref)**2)/(length(subsub$h)-1)
      if(nrow(subsub) > 1) {stats = rbind(stats, data.frame(family=sub$family[1], stage=3, vhp = vest/(ref*(1-ref)))) }
    }
  }
  
  if(nrow(olines) > 0) {
    data.plot[[eref]] = ggplot() +
      geom_segment(data=lines, aes(x=0, y=h0, xend=1, yend=h, color=factor(family))) +
      geom_point(data=lines, aes(x=1, y=h, color=factor(family))) +
      geom_segment(data=olines, aes(x=1.25, y=h0, xend=time/2+1, yend=h, color=factor(family))) +
      geom_point(data=olines, aes(x=time/2+1, y=h, color=factor(family))) +
      scale_x_continuous(name="stage", breaks = c(0, 1, 1.25, 1.5, 2, 2.5), labels=c("Mother EL", "Mean child EL", "Estimate child 0", "Child EL", "Child LL", "Child INF"), limits = c(0,2.5)) +
      scale_y_continuous(name="h", limits = c(0,1)) + theme_classic()+
      theme(axis.text.x = element_text(angle = 90), legend.position="none") +
      coord_flip() 
    
  }       else {
    data.plot[[eref]] = ggplot() +
      geom_segment(data=lines, aes(x=0, y=h0, xend=1, yend=h, color=factor(family))) +
      geom_point(data=lines, aes(x=1, y=h, color=factor(family))) +
      scale_x_continuous(name="stage", breaks = c(0, 1, 1.25, 1.5, 2, 2.5), labels=c("Mother EL", "Mean child EL", "Estimate child 0", "Child EL", "Child LL", "Child INF"), limits = c(0,2.5)) +
      scale_y_continuous(name="h", limits = c(0,1)) + theme_classic() +
      theme(axis.text.x = element_text(angle = 90), legend.position="none") +
      coord_flip() 
  }
  
}


# MCMC characteristics
steps = 100000; burnin = 5000; burnskip = 10;

fulldf = data.frame()
# loop through experiments
for(i in 1:length(expts)) {
  expt = expts[i]
  prior = priors[i]
  
  # read data from two parallel output files
  df1 = read.csv(paste(c("outtrjswitch-", expt, ".csv-1-50-100-", prior, "-100000.csv"), collapse=""))
  df1$seed = 1
  df1$index = 1:nrow(df1)

  # construct overall data frame
  df = rbind(df1[seq(burnin,nrow(df1),burnskip),])
  nref = which(colnames(df)=="n0")
  df$expt = expt

  # set experiment-specific plot parameters
  # effective population size (pt 7, mt 50)
  ne.new = mtplot.ne 
  # x-axis (broader for wildtype)
  xmin= -5; xmax = defaultplot.max 
  # model slice to focus on for first plot set
  modelslice = 0
  # number of bins for histograms
  nbins = 15
  
  # transform inferred division count to suit this population size
  df$n0p = mytrans(df$n0, ne.new=ne.new)
  df$n1p = mytrans(df$n1, ne.new=ne.new)
  df$n2p = mytrans(df$n2, ne.new=ne.new)
  df$n3p = mytrans(df$n3, ne.new=ne.new)
  

  # gather effective division counts for each developmental stage (respecting model index)
  df$eff1 = df$eff2 = df$eff3 = df$eff4 = 0
  for(j in 1:nrow(df)) {
    if(df$model[j] == 0) { df$eff1[j] = df$n0p[j]; df$eff2[j] = df$n0p[j]+df$n1p[j]; df$eff3[j] = df$n0p[j]+df$n1p[j]+df$n2p[j]; df$eff4[j] = df$n0p[j]+df$n1p[j]+df$n2p[j]+df$n3p[j] }
    if(df$model[j] == 1) { df$eff1[j] = df$n0p[j]; df$eff2[j] = df$n0p[j]+df$n1p[j]; df$eff3[j] = df$n2p[j]; df$eff4[j] = df$n2p[j]+df$n3p[j] }
    if(df$model[j] == 2) { df$eff1[j] = df$n0p[j]; df$eff2[j] = df$n1p[j]; df$eff3[j] = df$n2p[j]; df$eff4[j] = df$n2p[j]+df$n3p[j] }
  }
  
  nbins = 10
  # construct second plot set (all model indices)
  subplots1[[i]] = list()
  # individual effective division counts
  subplots1[[i]][[1]] = data.plot[[i]]
  subplots1[[i]][[2]] = ggplot(df, aes(x=eff1, fill=factor(seed))) + geom_histogram(aes(y=..density..), position="dodge",bins=nbins) + xlim(xmin,xmax) +theme_classic() + theme(legend.position="none") 
  subplots1[[i]][[3]] = ggplot(df, aes(x=eff2, fill=factor(seed))) + geom_histogram(aes(y=..density..), position="dodge",bins=nbins) + xlim(xmin,xmax) +theme_classic() + theme(legend.position="none") 
  subplots1[[i]][[4]] = ggplot(df, aes(x=eff3, fill=factor(seed))) + geom_histogram(aes(y=..density..), position="dodge",bins=nbins) + xlim(xmin,xmax) +theme_classic() + theme(legend.position="none") 
  subplots1[[i]][[5]] = ggplot(df, aes(x=eff4, fill=factor(seed))) + geom_histogram(aes(y=..density..), position="dodge",bins=nbins) + xlim(xmin,xmax) +theme_classic() + theme(legend.position="none") 
  # model index
  subplots1[[i]][[6]] = ggplot(df, aes(x=model, fill=factor(seed))) + geom_histogram(position="dodge") + scale_x_continuous(breaks = c(0,1,2)) + theme_classic() + theme(legend.position="none") 
  subplots1[[i]][[7]] = ggplot(df, aes(x=index, y=lik)) + geom_line()  + theme_classic()
  
  
  df$exptindex = i
  fulldf = rbind(fulldf, subset(df, select=c(eff1, eff2, eff3, eff4, seed, exptindex)))
}

g.empty = ggplot() + theme_void()


# plot summaries and posteriors for each test bed set
png("fig-s1-test-all.png", width=1400*gsf, height=800*gsf, res=72*ressf)
grid.arrange(arrangeGrob(grobs = subplots1[[1]], left = "Test 1", nrow=1),
             arrangeGrob(grobs = subplots1[[2]], left = "Test 2", nrow=1),
                          arrangeGrob(grobs = subplots1[[3]], left = "Test 3", nrow=1),
                          arrangeGrob(grobs = subplots1[[4]], left = "Test 4", nrow=1),
                          arrangeGrob(grobs = subplots1[[5]], left = "Test 5", nrow=1),
             nrow = 5)
dev.off()

