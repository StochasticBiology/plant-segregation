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

defaultplot.max = 35
mtplot.max = 150

# set of vectors encoding different experimental setups -- used for experiment-specific tweaks to plotting style
# vector of filenames
expts = c("old-mito-MSH1",
          "old-mito-WILD", "new-mito-WILD", "all-mito-WILD",
          "old-plastid-MSH1", "new-plastid-MSH1", "all-plastid-MSH1", 
          "new-mito-WILD",
          "new-plastid-MSH1")
# vector of priors on developmental model (0 -- uniform, 1 -- linear dev only)
priors = c(0,
           0, 0, 0,
           0, 0, 0,
           -1, 
           -1)
# vector of organelle being studied (0 -- mt, 1 -- pt)
organelles = c(0,
               0, 0, 0,
               1, 1, 1,
               0,
               1)
# vector of genotypes (0 -- msh1, 1 -- wild)
wildtype = c(0,
             1, 1, 1,
             0, 0, 0,
             1,
             0)

# MCMC characteristics
steps = 100000; burnin = 1000; burnskip = 10;

fulldf = data.frame()
# loop through experiments
for(i in 1:length(expts)) {
  expt = expts[i]
  prior = priors[i]
  
  # read data from two parallel output files
  df1 = read.csv(paste(c("outtrj-newest-data-", expt, ".csv-1-50-100-", prior, "-100000.csv"), collapse=""))
  df1$seed = 1
  df1$index = 1:nrow(df1)
  df2 = read.csv(paste(c("outtrj-newest-data-", expt, ".csv-2-50-100-", prior, "-100000.csv"), collapse=""))
  df2$seed = 2
  df2$index = 1:nrow(df2)
  
  # construct overall data frame
  df = rbind(df1[seq(burnin,nrow(df1),burnskip),], df2[seq(burnin,nrow(df2),burnskip),])
  nref = which(colnames(df)=="n0")
  df$expt = expt
  df$prior = prior
  df$organelle = organelles[i]
  df$wildtype = wildtype[i]
  
  # set experiment-specific plot parameters
  # effective population size (pt 7, mt 50)
  if(organelles[i] == 1) { ne.new = ptplot.ne } else { ne.new = mtplot.ne }
  # x-axis (broader for wildtype)
  if(wildtype[i] == 1) { xmin = -25; xmax = mtplot.max } else { xmin= -5; xmax = defaultplot.max }
  # model slice to focus on for first plot set
  modelslice = 0
  # number of bins for histograms
  nbins = 15
  
  # transform inferred division count to suit this population size
  df$n0p = mytrans(df$n0, ne.new=ne.new)
  df$n1p = mytrans(df$n1, ne.new=ne.new)
  df$n2p = mytrans(df$n2, ne.new=ne.new)
  df$n3p = mytrans(df$n3, ne.new=ne.new)
  
  # construct first plot set (focussed on modelslice)
  subplots0[[i]] = list()
  # individual division counts
  subplots0[[i]][[1]] = ggplot(df[df$model==modelslice,], aes(x=n0p, fill=factor(seed))) + geom_histogram(aes(y=..density..), position="dodge",bins=nbins) + xlim(xmin,xmax) +theme_classic() + theme(legend.position="none")
  subplots0[[i]][[2]] = ggplot(df[df$model==modelslice,], aes(x=n1p, fill=factor(seed))) + geom_histogram(aes(y=..density..), position="dodge",bins=nbins) + xlim(xmin,xmax) +theme_classic() + theme(legend.position="none")
  subplots0[[i]][[3]] = ggplot(df[df$model==modelslice,], aes(x=n2p, fill=factor(seed))) + geom_histogram(aes(y=..density..), position="dodge",bins=nbins) + xlim(xmin,xmax) +theme_classic() + theme(legend.position="none")
  subplots0[[i]][[4]] = ggplot(df[df$model==modelslice,], aes(x=n3p, fill=factor(seed))) + geom_histogram(aes(y=..density..), position="dodge",bins=nbins) + xlim(xmin,xmax) +theme_classic() + theme(legend.position="none")
  # summed division count
  subplots0[[i]][[5]] = ggplot(df[df$model==modelslice,], aes(x=n0p+n1p+n2p+n3p, fill=factor(seed))) + geom_histogram(aes(y=..density..), position="dodge",bins=nbins) + xlim(2*xmin,2*xmax) +theme(legend.position="none")  + theme_classic() 
  # model index
  subplots0[[i]][[6]] = ggplot(df, aes(x=model, fill=factor(seed))) + geom_histogram(position="dodge") +theme(legend.position="none")  + theme_classic() 
  # likelihood trace as proxy for mixing
  subplots0[[i]][[7]] = ggplot(df, aes(x=index, y=lik)) + geom_line()  + theme_classic()
  
  # gather effective division counts for each developmental stage (respecting model index)
  df$eff1 = df$eff2 = df$eff3 = df$eff4 = 0
  for(j in 1:nrow(df)) {
    if(df$model[j] == 0) { df$eff1[j] = df$n0p[j]; df$eff2[j] = df$n0p[j]+df$n1p[j]; df$eff3[j] = df$n0p[j]+df$n1p[j]+df$n2p[j]; df$eff4[j] = df$n0p[j]+df$n1p[j]+df$n2p[j]+df$n3p[j] }
    if(df$model[j] == 1) { df$eff1[j] = df$n0p[j]; df$eff2[j] = df$n0p[j]+df$n1p[j]; df$eff3[j] = df$n2p[j]; df$eff4[j] = df$n2p[j]+df$n3p[j] }
    if(df$model[j] == 2) { df$eff1[j] = df$n0p[j]; df$eff2[j] = df$n1p[j]; df$eff3[j] = df$n2p[j]; df$eff4[j] = df$n2p[j]+df$n3p[j] }
  }
  
  nbins = 20
  # construct second plot set (all model indices)
  subplots1[[i]] = list()
  # individual effective division counts
  subplots1[[i]][[1]] = ggplot(df, aes(x=eff1, fill=factor(seed))) + geom_histogram(aes(y=..density..), position="dodge",bins=nbins) + xlim(xmin,xmax) +theme_classic() + theme(legend.position="none") 
  subplots1[[i]][[2]] = ggplot(df, aes(x=eff2, fill=factor(seed))) + geom_histogram(aes(y=..density..), position="dodge",bins=nbins) + xlim(xmin,xmax) +theme_classic() + theme(legend.position="none") 
  subplots1[[i]][[3]] = ggplot(df, aes(x=eff3, fill=factor(seed))) + geom_histogram(aes(y=..density..), position="dodge",bins=nbins) + xlim(xmin,xmax) +theme_classic() + theme(legend.position="none") 
  subplots1[[i]][[4]] = ggplot(df, aes(x=eff4, fill=factor(seed))) + geom_histogram(aes(y=..density..), position="dodge",bins=nbins) + xlim(xmin,xmax) +theme_classic() + theme(legend.position="none") 
  # model index
  subplots1[[i]][[5]] = ggplot(df, aes(x=model, fill=factor(seed))) + geom_histogram(position="dodge") + theme_classic() + theme(legend.position="none") 

  # construct third plot set (all model indices, for just one seed)
  subplots2[[i]] = list()
  # individual effective division counts
  subplots2[[i]][[1]] = ggplot(df[df$seed==1,], aes(x=eff1)) + geom_histogram(aes(y=..density..), position="dodge",bins=nbins) + xlim(xmin,xmax) +theme_classic() + theme(legend.position="none") + xlab("O -> EL")
  subplots2[[i]][[2]] = ggplot(df[df$seed==1,], aes(x=eff2)) + geom_histogram(aes(y=..density..), position="dodge",bins=nbins) + xlim(xmin,xmax) +theme_classic() + theme(legend.position="none") + xlab("O -> LL")
  subplots2[[i]][[3]] = ggplot(df[df$seed==1,], aes(x=eff3)) + geom_histogram(aes(y=..density..), position="dodge",bins=nbins) + xlim(xmin,xmax) +theme_classic() + theme(legend.position="none") + xlab("O -> INF")
  subplots2[[i]][[4]] = ggplot(df[df$seed==1,], aes(x=eff4)) + geom_histogram(aes(y=..density..), position="dodge",bins=nbins) + xlim(xmin,xmax) +theme_classic() + theme(legend.position="none") + xlab("OM -> O")
  # model index
  
  df$exptindex = i
  fulldf = rbind(fulldf, subset(df, select=c(eff1, eff2, eff3, eff4, seed, exptindex)))
}

g.empty = ggplot() + theme_void()

# plot just those posteriors constrained to the linear development model
#grid.arrange(grobs = c(subplots0[[1]], subplots0[[8]], subplots0[[9]]), nrow=3)

# plot old-new-all data
m.g1 = ggplot(fulldf[fulldf$exptindex %in% c(2,3,4) & fulldf$seed == 1,], aes(x = eff1, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none") + xlim(0,mtplot.max) + xlab("O -> EL")
m.g2 = ggplot(fulldf[fulldf$exptindex %in% c(2,3,4) & fulldf$seed == 1,], aes(x = eff2, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none") + xlim(0,mtplot.max) + xlab("O -> LL")
m.g3 = ggplot(fulldf[fulldf$exptindex %in% c(2,3,4) & fulldf$seed == 1,], aes(x = eff3, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none") + xlim(0,mtplot.max) + xlab("O -> INF")
m.g4 = ggplot(fulldf[fulldf$exptindex %in% c(2,3,4) & fulldf$seed == 1,], aes(x = eff4, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none") + xlim(0,mtplot.max) + xlab("OM -> O")
p.g1 = ggplot(fulldf[fulldf$exptindex %in% c(5,6,7) & fulldf$seed == 1,], aes(x = eff1, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none") + xlim(0,10) + xlab("O -> EL")
p.g2 = ggplot(fulldf[fulldf$exptindex %in% c(5,6,7) & fulldf$seed == 1,], aes(x = eff2, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none")+ xlim(0,10) + xlab("O -> LL")
p.g3 = ggplot(fulldf[fulldf$exptindex %in% c(5,6,7) & fulldf$seed == 1,], aes(x = eff3, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none")+ xlim(0,10) + xlab("O -> INF")
p.g4 = ggplot(fulldf[fulldf$exptindex %in% c(5,6,7) & fulldf$seed == 1,], aes(x = eff4, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none")+ xlim(0,50) + xlab("OM -> O")

png("oldnewdata-alt.png", width=1000, height=400)
grid.arrange(m.g1, m.g2, m.g3, m.g4, p.g1, p.g2, p.g3, p.g4, nrow = 2)
dev.off()

png("oldnewdata-alt-pt.png", width=1000, height=250)
grid.arrange( p.g1, p.g2, p.g3, p.g4, nrow = 1)
dev.off()

# plot effective division counts across models
png("all9.png", width=1000, height=1500)
grid.arrange(arrangeGrob(grobs = subplots1[[1]], left = "old-mito-msh1", nrow=1),
             arrangeGrob(grobs = subplots1[[2]], left = "old-mito-wild", nrow=1),
             arrangeGrob(grobs = subplots1[[3]], left = "new-mito-wild", nrow=1),
             arrangeGrob(grobs = subplots1[[4]], left = "all-mito-wild", nrow=1),
             arrangeGrob(grobs = subplots1[[5]], left = "old-pltd-msh1", nrow=1),
             arrangeGrob(grobs = subplots1[[6]], left = "new-pltd-msh1", nrow=1),
             arrangeGrob(grobs = subplots1[[7]], left = "all-pltd-msh1", nrow=1),
             arrangeGrob(grobs = subplots1[[8]], left = "new-mito-wild*", nrow=1),
             arrangeGrob(grobs = subplots1[[9]], left = "new-pltd-msh1*", nrow=1),
             nrow = 9)
dev.off()

subplots3 = list()
subplots3[[1]] = list()
subplots3[[1]][[1]] = subplots2[[1]][[4]]
subplots3[[1]][[2]] = subplots2[[1]][[1]]
subplots3[[1]][[3]] = subplots2[[1]][[2]]
subplots3[[1]][[4]] = subplots2[[1]][[3]]
subplots3[[2]] = list()
subplots3[[2]][[1]] = subplots2[[2]][[4]]
subplots3[[2]][[2]] = g.empty
subplots3[[2]][[3]] = g.empty
subplots3[[2]][[4]] = g.empty
subplots3[[3]] = list()
subplots3[[3]][[1]] = subplots2[[5]][[4]]
subplots3[[3]][[2]] = subplots2[[5]][[1]]
subplots3[[3]][[3]] = subplots2[[5]][[2]]
subplots3[[3]][[4]] = subplots2[[5]][[3]]

# plot effective division counts across models (just one seed)
png("olddata-old.png", width=1000, height=500)
grid.arrange(arrangeGrob(grobs = subplots2[[1]], left = "MT msh1", nrow=1),
             arrangeGrob(grobs = subplots2[[2]], left = "MT wildtype", nrow=1),
             arrangeGrob(grobs = subplots2[[5]], left = "PT msh1", nrow=1),
             nrow = 3)
dev.off()
png("olddata.png", width=1000, height=500)
grid.arrange(arrangeGrob(grobs = subplots3[[1]], left = "MT msh1", nrow=1),
             arrangeGrob(grobs = subplots3[[3]], left = "PT msh1", nrow=1),
             arrangeGrob(grobs = subplots3[[2]], left = "MT wildtype", nrow=1),
             
             nrow = 3)
dev.off()

# plot effective division counts across models (just one seed)
png("oldandnewdata.png", width=1000, height=650)
grid.arrange(arrangeGrob(grobs = subplots2[[2]], left = "MT wildtype", nrow=1),
             arrangeGrob(grobs = subplots2[[3]], left = "MT wildtype (new)", nrow=1),
             arrangeGrob(grobs = subplots2[[5]], left = "PT msh1", nrow=1),
             arrangeGrob(grobs = subplots2[[6]], left = "PT msh1 (new)", nrow=1),
             nrow = 4)
dev.off()

# try MT predictions based on scaling from MSH1 to WT
# plot old-new-all data
sf = 7
pred.g1 = ggplot(fulldf[fulldf$exptindex == 1 & fulldf$seed == 1,], aes(x = eff1*sf, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none") + xlim(0,mtplot.max)
pred.g2 = ggplot(fulldf[fulldf$exptindex == 1 & fulldf$seed == 1,], aes(x = eff2*sf, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none") + xlim(0,mtplot.max)
pred.g3 = ggplot(fulldf[fulldf$exptindex == 1 & fulldf$seed == 1,], aes(x = eff3*sf, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none") + xlim(0,mtplot.max)
pred.g4 = ggplot(fulldf[fulldf$exptindex == 1 & fulldf$seed == 1,], aes(x = eff4*sf, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none") + xlim(0,mtplot.max)

pred1.g1 = ggplot(fulldf[fulldf$exptindex == 1 & fulldf$seed == 1,], aes(x = eff1*2, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none") + xlim(0,mtplot.max)
pred1.g2 = ggplot(fulldf[fulldf$exptindex == 1 & fulldf$seed == 1,], aes(x = eff2*2, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none") + xlim(0,mtplot.max)
pred1.g3 = ggplot(fulldf[fulldf$exptindex == 1 & fulldf$seed == 1,], aes(x = eff3*2, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none") + xlim(0,mtplot.max)
pred1.g4 = ggplot(fulldf[fulldf$exptindex == 1 & fulldf$seed == 1,], aes(x = eff4*2, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none") + xlim(0,mtplot.max)

true.g1 = ggplot(fulldf[fulldf$exptindex == 3 & fulldf$seed == 1,], aes(x = eff1, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none") + xlim(0,mtplot.max)
true.g2 = ggplot(fulldf[fulldf$exptindex == 3 & fulldf$seed == 1,], aes(x = eff2, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none") + xlim(0,mtplot.max)
true.g3 = ggplot(fulldf[fulldf$exptindex == 3 & fulldf$seed == 1,], aes(x = eff3, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none") + xlim(0,mtplot.max)
true.g4 = ggplot(fulldf[fulldf$exptindex == 3 & fulldf$seed == 1,], aes(x = eff4, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none") + xlim(0,mtplot.max)

png("mt-pred.png", width=1000, height=400)
#grid.arrange(pred.g1, pred.g2, pred.g3, pred.g4, pred1.g1, pred1.g2, pred1.g3, pred1.g4, true.g1, true.g2, true.g3, true.g4, nrow = 3)
grid.arrange(pred.g1, pred.g2, pred.g3, pred.g4, true.g1, true.g2, true.g3, true.g4, nrow = 2)
dev.off()

both = fulldf[fulldf$exptindex == 1 & fulldf$seed == 1,]
both$eff1 = both$eff1*sf
both$eff2 = both$eff2*sf
both$eff3 = both$eff3*sf
both$eff4 = both$eff4*sf
both = rbind(both, fulldf[fulldf$exptindex == 3 & fulldf$seed == 1,])
both.g1 = ggplot(both, aes(x = eff1, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none") + xlim(0,mtplot.max) + xlab("O -> EL")
both.g2 = ggplot(both, aes(x = eff2, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none") + xlim(0,mtplot.max) + xlab("O -> LL")
both.g3 = ggplot(both, aes(x = eff3, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none") + xlim(0,mtplot.max) + xlab("O -> INF")
both.g4 = ggplot(both, aes(x = eff4, fill=factor(exptindex))) + geom_histogram(aes(y=..density..), alpha=0.5,position="identity",bins=nbins) + theme_classic() + theme(legend.position="none") + xlim(0,mtplot.max) + xlab("OM -> O")

png("mt-pred-alt.png", width=1000, height=250)
grid.arrange(both.g1, both.g2, both.g3, both.g4, nrow = 1)
dev.off()

png("new-pred-both.png", width=1000, height=300)
grid.arrange(both.g1, both.g2, both.g3, both.g4, p.g1, p.g2, p.g3, p.g4, nrow = 2)
dev.off()
