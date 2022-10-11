library(ggplot2)
library(gridExtra)

expts = c("newest-data-old-mito-MSH1", "newest-data-old-mito-WILD", "newest-data-old-plastid-MSH1", "newest-data-new-mito-WILD", "newest-data-new-plastid-MSH1", "newest-data-all-mito-WILD", "newest-data-all-plastid-MSH1")

var.plot = data.plot = list()

gsf = 3
ressf = 4.5
for(eref in 1:length(expts)) {

  expt = expts[eref]
  filename = paste(c(expt, ".csv"), collapse="")

    df = read.csv(filename, header=T)

  df$h = df$h/100
  
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
 	  
#  var.plot[[eref]] = ggplot() +
#      geom_boxplot(data=stats, aes(group=stage, y=vhp)) +
#                     geom_point(data=stats, aes(x=factor(stage), y=vhp, color=factor(family))) +
#                     scale_x_continuous(name="stage", breaks = c(0, 1, 2, 3), labels=c("Cross-gen", "Early", "Late", "Inflor")) +
#                     scale_y_continuous(name="V'(h)") +
#                     theme(axis.text.x = element_text(angle = 90)) +
#		     labs(color = "Germline\nmodel")
  stages = unique(stats$stage)
  means = data.frame()
  for(stage in stages) {
    means = rbind(means, data.frame(stage = stage, meanvh = mean(stats$vhp[stats$stage == stage])))
  }
  var.plot[[eref]] = ggplot(data=stats, aes(x=factor(stage), y=vhp)) + geom_boxplot(alpha=0.3, col="#CCCCCC") + geom_point(size=3,aes(color=factor(family))) + 
    geom_point(data=means, aes(x=factor(stage), y=meanvh)) +
    scale_x_discrete(name="stage", breaks = c(0, 1, 2, 3), labels=c("Cross-gen", "EL", "LL", "INF")) + theme_classic() +
                        scale_y_continuous(name="V'(h)") + ylim(0,2) +
                         theme(axis.text.x = element_text(angle = 90)) +
    		     labs(color = "family") + coord_flip()
    
}

empty.plot = ggplot() + theme_void()

png("all-data-plots-newest.png", width=1000*gsf, height=1800*gsf, res=72*ressf)

#grid.arrange(data.plot[[1]], data.plot[[2]], data.plot[[3]], data.plot[[4]], data.plot[[5]], nrow=5)
grid.arrange(data.plot[[1]], data.plot[[2]], data.plot[[3]], data.plot[[4]], nrow=4)
dev.off()

png("old-data-plots-var-newest.png", width=1000*gsf, height=600*gsf, res=72*ressf)
grid.arrange(data.plot[[1]], var.plot[[1]], data.plot[[2]], var.plot[[2]], data.plot[[3]], var.plot[[3]], nrow=3)
dev.off()

png("new-data-plots-var-newest.png", width=1000*gsf, height=400*gsf, res=72*ressf)
grid.arrange(data.plot[[4]], var.plot[[4]], data.plot[[5]], var.plot[[5]], nrow=2)
dev.off()

png("all-data-plots-var-newest.png", width=1000*gsf, height=600*gsf, res=72*ressf)
grid.arrange(data.plot[[1]], var.plot[[1]], data.plot[[6]], var.plot[[6]], data.plot[[7]], var.plot[[7]], nrow=3)
dev.off()