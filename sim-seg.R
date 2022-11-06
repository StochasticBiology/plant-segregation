### produce supplementary plots relating to influence of gene conversion on heteroplasmy distributions

library(ggplot2)
library(gridExtra)
library(see)

# read output from simulation
df = read.csv("sim-seg.csv")

# oDNA population size in model
npop = 50

# subset out content, normalising gene conversion counts by oDNA population size
dfr = data.frame()
for(ngc in unique(df$ngc)) {
  for(div in unique(df$div)) {
    sub = df[df$ngc == ngc & df$div == div,]
    dfr = rbind(dfr, data.frame(ngc = ngc, div=div, vh = var(sub$h), meangc = mean(sub$realgc)/npop/div))
  }
}

# first plot -- heteroplasmy variance with divisions and gene conversion rates
ga1 = ggplot(dfr, aes(x = div, y = vh, col = factor(ngc))) + geom_line() + xlab("Divisions") + ylab("V(h)") + labs(col = "GC rate\n(AU)") + ylim(0.20,0.252) + xlim(0,100) + theme_classic()
ga2 = ggplot(dfr, aes(x = div, y = meangc, col = factor(ngc))) + geom_line() + xlab("Divisions") + ylab("Mean GC events\nper oDNA") + labs(col = "GC rate\n(AU)") + xlim(0,100) + theme_classic()
sf = 3
png("fig-s4-sim-seg.png", width=400*sf, height=150*sf, res=72*sf)
grid.arrange(ga1, ga2, nrow=1)
dev.off()

# second plot -- predicted heteroplasmy distributions for divisions and gene conversion rates
sub = df[(df$div == 1 | df$div == 10 | df$div == 20 | df$div == 50 ),] 
png("fig-s5-sim-pred.png", width=300*sf, height=250*sf, res=72*sf)
ggplot(sub[sub$ngc %in% c(0, 400),], aes(x=h)) + geom_histogram(binwidth=0.1) + facet_grid(div~ngc) + theme_classic()
dev.off()

# data frame to illustate simplified conversion factor between effective population sizes
x = 1:200
df = data.frame(new.n = x, old.n = 50, scale.factor = log(1-1/50)/log(1-1/x))
df = rbind(df, data.frame(new.n = x, old.n = 7, scale.factor = log(1-1/7)/log(1-1/x)))

df$ratio = df$scale.factor / (df$new.n/df$old.n)

png("fig-s6-scalepopsize.png", width=200*sf, height=200*sf, res=72*sf)
ggplot(df, aes(x=new.n, y=scale.factor, col=factor(old.n))) + geom_line() + xlab("New Ne'") + ylab("Scale") + labs(col="Old Ne") + theme_light()
dev.off()


