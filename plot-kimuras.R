library(kimura)

nbin = 10
x = (0:nbin)/nbin
results = data.frame()
for(nb in 1:14) {
  pdf = dkimura(x, 0.5, 1-1/nb)
  pdf[2:nbin] = pdf[2:nbin]*(1/(nbin-1))
  results = rbind(results, data.frame(nb=nb,x=x,pdf=pdf))
}
png("kimura-dists.png", width=1000, height=1000)
ggplot(results, aes(x=x,y=pdf)) + geom_col() + facet_wrap(~nb) + theme_void()
dev.off()


df = read.csv("dbtrj-50-100.csv")
ggplot(df[df$h0 == 0.5 & df$n %in% c(100,1,2,5),], aes(x=h,y=Ph)) + geom_col() + facet_wrap(~ n)
