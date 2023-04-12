# transform effective cell division counts, given an old and a new population size
# simulations are run with n=50. we might want, for example, n=200 for mitos or n=10 for plastids
mytrans = function(n, ne.orig=50, ne.new=100) {
  return(n*log(1-1/ne.orig)/log(1-1/ne.new))
}

# read in pre-existing PDF database (make sure simulation code has been run to provide this file)
df = read.csv("dbtrj-50-100.csv")
png("kimuras-mt.png", width=600, height=400)
# plot distributions relevant for MT behaviour
ggplot(df[df$h0 == 0.5 & df$n %in% c(4,6,8,12,25,50,80),], aes(x=h,y=Ph)) + 
  geom_col(width=0.05) + facet_wrap(~ n, scales="free") + theme_void()
dev.off()

mytrans(c(1.5, 1, 0.5, 18), ne.orig=7, ne.new=50)
png("kimuras-pt.png", width=400, height=400)
# plot distributions relevant for PT behaviour
ggplot(df[df$h0 == 0.5 & df$n %in% c(4,8,12,137),], aes(x=h,y=Ph)) + 
  geom_col(width=0.05) + facet_wrap(~ n, scales="free") + theme_void()
dev.off()
