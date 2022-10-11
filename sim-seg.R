library(ggplot2)
library(gridExtra)
library(see)

df = read.csv("sim-seg.csv")

dfr = data.frame()
for(ngc in unique(df$ngc)) {
  for(div in unique(df$div)) {
    sub = df[df$ngc == ngc & df$div == div,]
    dfr = rbind(dfr, data.frame(ngc = ngc, div=div, vh = var(sub$h), meangc = mean(sub$realgc)/200/div))
  }
}

ga1 = ggplot(dfr, aes(x = div, y = vh, col = factor(ngc))) + geom_line() + xlab("Divisions") + ylab("V(h)") + labs(col = "GC rate\n(AU)") + ylim(0.20,0.252) + xlim(0,100) + theme_classic()
ga2 = ggplot(dfr, aes(x = div, y = meangc, col = factor(ngc))) + geom_line() + xlab("Divisions") + ylab("Mean GC events\nper oDNA") + labs(col = "GC rate\n(AU)") + xlim(0,100) + theme_classic()
sf = 3
png("sim-seg.png", width=400*sf, height=150*sf, res=72*sf)
grid.arrange(ga1, ga2, nrow=1)
dev.off()

sub = df[(df$div == 1 | df$div == 10 | df$div == 20 | df$div == 50 ),] 
png("sim-pred.png", width=300*sf, height=250*sf, res=72*sf)
ggplot(sub[sub$ngc %in% c(0, 400),], aes(x=h)) + geom_histogram(binwidth=0.1) + facet_grid(div~ngc) + theme_classic()
dev.off()


