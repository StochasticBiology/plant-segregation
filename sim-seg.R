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

ga1 = ggplot(dfr, aes(x = div, y = vh, col = factor(ngc))) + geom_line() + ylim(0.20,0.252) + xlim(0,100)
ga2 = ggplot(dfr, aes(x = div, y = meangc, col = factor(ngc))) + geom_line() + xlim(0,100)
sf = 3
png("sim-seg.png", width=400*sf, height=150*sf, res=72*sf)
grid.arrange(ga1, ga2, nrow=1)
dev.off()

sub = df[(df$div == 1 | df$div == 10 | df$div == 20 | df$div == 50 | df$div == 100),] 
png("sim-pred.png", width=400*sf, height=350*sf, res=72*sf)
ggplot(sub[sub$ngc %in% c(0, 400),], aes(x=factor(div), y=h)) + geom_violinhalf() + facet_wrap(~ngc, ncol=1)
dev.off()

sub = df[df$div == 39 | df$div == 299,]
g1 = ggplot(sub, aes(x=h)) + geom_histogram() + facet_grid(div ~ ngc)
g2 = ggplot(sub, aes(x=div, y=realgc/200/div)) + geom_jitter() + facet_wrap(~ngc)
grid.arrange(g1, g2)

ggplot(sub, aes(x = div, y=h)) + geom_jitter() + facet_wrap(~ngc)
