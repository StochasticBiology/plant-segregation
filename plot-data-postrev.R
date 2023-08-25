### plot parsed data with "fan" diagrams and summary statistics

library(ggplot2)
library(gridExtra)

# set of the individual datafiles
expts = c("newest-data-old-mito-MSH1", "newest-data-old-mito-WILD", "newest-data-old-plastid-MSH1", "newest-data-new-mito-WILD", "newest-data-new-plastid-MSH1", "newest-data-all-mito-WILD", "newest-data-all-plastid-MSH1")

# initialise lists storing plots
var.plot = data.plot = mean.plot = vartraj.plot = list()

# dimension and resolution scale factors for output images
gsf = 3
ressf = 4.5

# loop through datafiles
for(eref in 1:length(expts)) {
  # read and normalise this datafile  
  expt = expts[eref]
  filename = paste(c(expt, ".csv"), collapse="")
  df = read.csv(filename, header=T)
  df$h = df$h/100

  # first pull the family-wise observations from this dataset
  # !!! NOTE -- for eref == 5 the mother's initial heteroplasmy is added manually -- to keep the original data intact
  lines = data.frame()
  families = unique(df$family)
  for(family in families) {
    sub = df[df$family == family,]
    ref = sub$h[which(sub$time == 0)]
    if(eref == 5) { ref = 0.14 }  # manually include mother reference for family examined in new plastid data
    for(offspring in unique(sub$id)) {
      meanh = mean(sub$h[sub$id == offspring & sub$time == 1])
      if(!is.na(meanh)) {
        lines = rbind(lines, data.frame(family=family,id=offspring,h0=ref,h=meanh))
      }
    }
  }

  # next build a dataframe storing heteroplasmy mean and normalised variance for these observations
  stats = data.frame()
  for(family in families) {
    sub = lines[lines$family == family,]
    hs = sub$h
    hdiffs = (hs-mean(hs))**2
    vest = sum(hdiffs)/(length(hs)-1)
    stats = rbind(stats, data.frame(family=family, id=-1, stage=0, mest=mean(hs), vhp = vest/(mean(hs)*(1-mean(hs)))))
  }

  # dataframe of individual heteroplasmy points
  points = data.frame()
  for(id in 1:max(df$id)) {
    sub = df[df$id == id,]
    if(nrow(sub) > 1) {
      points = rbind(points, data.frame(family=sub$family, id=id, time=sub$time+1, h=sub$h))
    }
  }

  # rather awkward code gathering stats and drawing developmental lines between individual points
  # loop through individual IDs, pulling mean and variance estimates for all samples at a given developmental stage
  olines = data.frame()
  for(id in 1:max(df$id)) {
    sub = df[df$id == id,]
    if(nrow(sub) > 1) {
      # reference is average over all tissues
      ref = mean(sub$h[sub$time>=1])
      olines = rbind(olines, data.frame(family=sub$family, id=id, h0=ref, time=sub$time, h=sub$h))
      # process early leaves
      subsub = sub[sub$time==1,]
      vest = sum((subsub$h-ref)**2)/(length(subsub$h)-1)
      mest = mean(subsub$h)
      if(nrow(subsub) > 1) {stats = rbind(stats, data.frame(family=sub$family[1], id=id, stage=1, mest=mest, vhp = vest/(ref*(1-ref)))) }
      # process late leaves
      subsub = sub[sub$time==2,]
      vest = sum((subsub$h-ref)**2)/(length(subsub$h)-1)
      mest = mean(subsub$h)
      if(nrow(subsub) > 1) {stats = rbind(stats, data.frame(family=sub$family[1], id=id, stage=2, mest=mest, vhp = vest/(ref*(1-ref)))) }
      # process inflorescences
      subsub = sub[sub$time==3,]
      vest = sum((subsub$h-ref)**2)/(length(subsub$h)-1)
      mest = mean(subsub$h)
      if(nrow(subsub) > 1) {stats = rbind(stats, data.frame(family=sub$family[1], id=id, stage=3, mest=mest, vhp = vest/(ref*(1-ref)))) }
    }
  }

  # now build up the "fan" plots using sample points and connecting lines
  # other faff labels and structures axes
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
    
  } else {
    data.plot[[eref]] = ggplot() +
      geom_segment(data=lines, aes(x=0, y=h0, xend=1, yend=h, color=factor(family))) +
      geom_point(data=lines, aes(x=1, y=h, color=factor(family))) +
      scale_x_continuous(name="stage", breaks = c(0, 1, 1.25, 1.5, 2, 2.5), labels=c("Mother EL", "Mean child EL", "Estimate child 0", "Child EL", "Child LL", "Child INF"), limits = c(0,2.5)) +
      scale_y_continuous(name="h", limits = c(0,1)) + theme_classic() +
      theme(axis.text.x = element_text(angle = 90), legend.position="none") +
      coord_flip() 
  }

  # quick data frame storing mean heteroplasmy variances across samples and for each stage
  stages = unique(stats$stage)
  means = data.frame()
  for(stage in stages) {
    means = rbind(means, data.frame(stage = stage, meanvh = mean(stats$vhp[stats$stage == stage], na.rm=T)))
  }

  # variance plots: boxplots and points for normalised heteroplasmy variances by stage
  var.plot[[eref]] = ggplot(data=stats, aes(x=factor(stage), y=vhp)) +
    geom_boxplot(alpha=0.3, col="#CCCCCC") + geom_point(aes(color=factor(family))) + 
    geom_point(data=means, aes(x=factor(stage), y=meanvh)) +
    scale_x_discrete(name="stage", breaks = c(0, 1, 2, 3), labels=c("M→C", "EL", "LL", "INF")) + theme_classic() +
    scale_y_continuous(name="V'(h)") +
    theme(axis.text.x = element_text(angle = 90), legend.position="none") + ylab("V'(h)") +
    labs(color = "family") + coord_flip()

  # mean trajectory plot: lines and points for mean heteroplasmy by stage, for each individual
  mean.plot[[eref]] = ggplot(stats[stats$id != -1,], aes(x=mest,y=stage,colour=factor(id))) + 
    geom_line() + geom_point() + xlim(0,1) + 
    scale_y_continuous(name="stage", breaks = c(1, 2, 3), labels=c("EL", "LL", "INF")) +
    theme_classic() + xlab("Individual mean h") + theme(axis.text.x = element_text(angle = 90), legend.position="none")

  # variance trajectory plot: lines and points for variance heteroplasmy by stage, for each individual
  vartraj.plot[[eref]] = ggplot(stats[stats$id != -1,], aes(x=vhp,y=stage,colour=factor(id))) + 
    geom_line() + geom_point() + xlim(0,1) + 
    scale_y_continuous(name="stage", breaks = c(1, 2, 3), labels=c("EL", "LL", "INF")) +
    theme_classic() + xlab("Individual V'(h)") + theme(axis.text.x = element_text(angle = 90), legend.position="none")
}

g.empty = ggplot() + theme_void()
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette = c(cbPalette[-(0:6)], cbPalette[0:6])
#data.plot[[1]] + scale_colour_manual(values=cbPalette)

png("fig-3ab-new-postrev.png", width=820*gsf, height=600*gsf, res=72*ressf)
grid.arrange(
             arrangeGrob(grobs = list(data.plot[[2]] + scale_colour_brewer(palette="Set2")), left = "Old, MT wild", nrow=1), 
             arrangeGrob(grobs = list(data.plot[[4]] + scale_colour_brewer(palette="Set2")), left = "New, MT wild", nrow=1), 
             arrangeGrob(grobs = list(data.plot[[3]] + scale_colour_brewer(palette="Set2")), left = "Old, PT msh1", nrow=1), 
             arrangeGrob(grobs = list(data.plot[[5]] + scale_colour_brewer(palette="Set2")), left = "New, PT msh1", nrow=1), 
             nrow=2)
dev.off()

# figure 1 -- existing heteroplasmy measurements, fans only
png("fig-1-postrev.png", width=820*gsf, height=400*gsf, res=72*ressf)
grid.arrange(arrangeGrob(grobs = list(data.plot[[1]] + scale_colour_brewer(palette="Set2")), left = "MT msh1", nrow=1), 
             arrangeGrob(grobs = list(data.plot[[2]] + scale_colour_brewer(palette="Set2")), left = "MT wild", nrow=1), 
             arrangeGrob(grobs = list(data.plot[[3]] + scale_colour_brewer(palette="Set2")), left = "PT msh1", nrow=1), 
             arrangeGrob(grobs = list(g.empty), left = "PT wild", nrow=1),
             nrow=2)
dev.off()

# figure 3A-B -- new heteroplasmy measurements, fans and variances
png("fig-3ab-postrev.png", width=820*gsf, height=400*gsf, res=72*ressf)
grid.arrange(arrangeGrob(grobs = list(data.plot[[2]]+ scale_colour_brewer(palette="Set2")), left = "MT wild, old", nrow=1),
             arrangeGrob(grobs = list(data.plot[[4]]+ scale_colour_brewer(palette="Set2")), left = "MT wild, new", nrow=1),
             arrangeGrob(grobs = list(data.plot[[3]]+ scale_colour_brewer(palette="Set2")), left = "PT msh1, old", nrow=1),
             arrangeGrob(grobs = list(data.plot[[5]]+ scale_colour_brewer(palette="Set2")), left = "PT msh1, new", nrow=1),
             nrow=2)
dev.off()

# figure 1 -- existing heteroplasmy measurements, fans and variances
png("old-fig-1-old-data-plots-var-newest.png", width=820*gsf, height=600*gsf, res=72*ressf)
grid.arrange(arrangeGrob(grobs = list(data.plot[[1]], var.plot[[1]]), left = "MT msh1", nrow=1), 
             arrangeGrob(grobs = list(data.plot[[2]], var.plot[[2]]), left = "MT wild", nrow=1), 
             arrangeGrob(grobs = list(data.plot[[3]], var.plot[[3]]), left = "PT msh1", nrow=1), 
             nrow=3)
dev.off()

# figure 3A-B -- new heteroplasmy measurements, fans and variances
png("old-fig-3ab-new-data-plots-var-newest.png", width=820*gsf, height=200*gsf, res=72*ressf)
grid.arrange(arrangeGrob(grobs = list(data.plot[[4]], var.plot[[4]]), left = "MT wild", nrow=1),
             arrangeGrob(grobs = list(data.plot[[5]], var.plot[[5]]), left = "PT msh1", nrow=1),
             nrow=2)
dev.off()

# figure 4 -- combined measurements, fans and statistic trajectories
png("fig-4-all-data-plots-trajs-newest.png", width=820*gsf, height=600*gsf, res=72*ressf)
grid.arrange(arrangeGrob(grobs = list(data.plot[[1]], vartraj.plot[[1]], mean.plot[[1]]), left = "MT msh1", widths = c(2, 1, 1), nrow=1),
            arrangeGrob(grobs = list(data.plot[[6]], vartraj.plot[[6]], mean.plot[[6]]), left = "MT wild", widths = c(2, 1, 1), nrow=1), 
            arrangeGrob(grobs = list(data.plot[[7]], vartraj.plot[[7]], mean.plot[[7]]), left = "PT msh1", widths = c(2, 1, 1), nrow=1), 
             nrow=3)
dev.off()

# old plots for legacy reasons
png("all-data-plots-newest.png", width=820*gsf, height=1800*gsf, res=72*ressf)
grid.arrange(data.plot[[1]], data.plot[[2]], data.plot[[3]], data.plot[[4]], nrow=4)
dev.off()

png("all-data-plots-var-newest.png", width=820*gsf, height=600*gsf, res=72*ressf)
grid.arrange(data.plot[[1]], var.plot[[1]], data.plot[[6]], var.plot[[6]], data.plot[[7]], var.plot[[7]], nrow=3)
dev.off()

png("all-data-plots-var-mean-newest.png", width=820*gsf, height=600*gsf, res=72*ressf)
grid.arrange(data.plot[[1]], var.plot[[1]], mean.plot[[1]], data.plot[[6]], var.plot[[6]], mean.plot[[6]], data.plot[[7]], var.plot[[7]], mean.plot[[7]], nrow=3, widths = c(2, 1, 1))
dev.off()

