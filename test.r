require(gaston, lib.loc="/env/export/cngstkprd1/v_home/q_unix/dandine/R/x86_64-pc-linux-gnu-library/3.3")
require(Hmisc, lib.loc="/env/export/cngstkprd1/v_home/q_unix/dandine/R/x86_64-pc-linux-gnu-library/3.3")
require(htmlTable, lib.loc="/env/export/cngstkprd1/v_home/q_unix/dandine/R/x86_64-pc-linux-gnu-library/3.3")
require(haplo.stats, lib.loc="/env/export/cngstkprd1/v_home/q_unix/dandine/R/x86_64-pc-linux-gnu-library/3.3")
require(gaston.pop, lib.loc="/env/export/cngstkprd1/v_home/q_unix/dandine/R/x86_64-pc-linux-gnu-library/3.3")
source('~/Programmes/mating.r')

n <- 50
K <- random.pm(n=100)$K
sex <- rep(c(1,2), n)
couple <- rep(1:n, each=2)

mean_couple(K, couple)
app.pairs(K, couple, mean)

sum( app_nonmates(K, couple)!= app.unpairs(K, couple) )

sum( app.pairs(K, couple)!=app_couple(K, couple) )

sum( app.unpairs(K, list(couple, sex))!=app_couple(K, couple, sex, opposite=T) )

z_score(K, couple, sex, method='all')
zscore(app.pairs(K, couple), m=app.unpairs(K, couple, mean), sd=app.unpairs(K, couple, sd))

z_score(K, couple, sex, method='opp_comp')
zscore(app.pairs(K, couple), m=app.unpairs(K, list(couple, sex), mean), sd=app.unpairs(K, list(couple, sex), sd))

z_score(K, couple, sex, method='opp_int')
zscore(app.unpairs(K, list(couple, sex)), m=app.pairs(K, sex, mean), sd=app.pairs(K, sex, sd))



source('~/Mixed Model/Package R Gaston/gaston.pop/R/permut_pairs.r')

app_permut(K, sex, couple, 10, method='mean', thread=1)
app.permut(K, couple, group=sex, mean, unpairs=FALSE, fixed.pairs=NULL, B=10, thread=1)
par(mfrow=c(1,2))
hist(app_permut(K, sex, couple, 1e3, method='mean', thread=1)$permut, breaks=30)
hist( app.permut(K, couple, group=sex, mean, unpairs=FALSE, fixed.pairs=NULL, B=1e3, thread=1)$permut, breaks=30)

app_permut(K, sex, couple, 10, method='mean', thread=2)
app.permut(K, couple, group=sex, mean, unpairs=FALSE, fixed.pairs=NULL, B=10, thread=2)
par(mfrow=c(1,2))
hist(app_permut(K, sex, couple, 1e3, method='mean', thread=2)$permut, breaks=30)
hist( app.permut(K, couple, group=sex, mean, unpairs=FALSE, fixed.pairs=NULL, B=1e3, thread=2)$permut, breaks=30)


app_permut(K, sex, couple, 10, method='zscore', thread=1)
app.permut(K, couple, group=sex, function(x) zscore(x, m=app.unpairs(K, couple, mean), sd=app.unpairs(K, couple, sd)), unpairs=FALSE, fixed.pairs=NULL, B=10, thread=1)
par(mfrow=c(1,2))
hist(app_permut(K, sex, couple, 1e3, method='zscore', thread=1)$permut, breaks=30)
hist( app.permut(K, couple, group=sex, function(x) zscore(x, m=app.unpairs(K, couple, mean), sd=app.unpairs(K, couple, sd)), unpairs=FALSE, fixed.pairs=NULL, B=1e3, thread=1)$permut, breaks=30)

app_permut(K, sex, couple, 10, method='zscore', thread=2)
app.permut(K, couple, group=sex, function(x) zscore(x, m=app.unpairs(K, couple, mean), sd=app.unpairs(K, couple, sd)), unpairs=FALSE, fixed.pairs=NULL, B=10, thread=2)
par(mfrow=c(1,2))
hist(app_permut(K, sex, couple, 1e3, method='zscore', thread=2)$permut, breaks=30)
hist( app.permut(K, couple, group=sex, function(x) zscore(x, m=app.unpairs(K, couple, mean), sd=app.unpairs(K, couple, sd)), unpairs=FALSE, fixed.pairs=NULL, B=1e3, thread=2)$permut, breaks=30)


app_permut(K, sex, couple, 10, method='zscore_coupleVSopposite', thread=1)
app.permut(K, couple, group=sex, FUN=function(x) zscore(x, m=app.unpairs(K, list(couple, sex), mean), sd=app.unpairs(K, list(couple, sex), sd)), B=10, thread=1)
par(mfrow=c(1,2))
hist(app_permut(K, sex, couple, 1e3, method='zscore_coupleVSopposite', thread=1)$permut, breaks=30)
hist( app.permut(K, couple, group=sex, FUN=function(x) zscore(x, m=app.unpairs(K, list(couple, sex), mean), sd=app.unpairs(K, list(couple, sex), sd)), B=1e3, thread=1)$permut, breaks=30)

app_permut(K, sex, couple, 10, method='zscore_coupleVSopposite', thread=2)
app.permut(K, couple, group=sex, FUN=function(x) zscore(x, m=app.unpairs(K, list(couple, sex), mean), sd=app.unpairs(K, list(couple, sex), sd)), B=10, thread=2)
par(mfrow=c(1,2))
hist(app_permut(K, sex, couple, 1e3, method='zscore_coupleVSopposite', thread=2)$permut, breaks=30)
hist( app.permut(K, couple, group=sex, FUN=function(x) zscore(x, m=app.unpairs(K, list(couple, sex), mean), sd=app.unpairs(K, list(couple, sex), sd)), B=1e3, thread=2)$permut, breaks=30)


app_permut(K, sex, couple, 10, method='zscore_oppositeVSall', thread=1)
app.permut(K, couple, group=sex, FUN=function(x) zscore(x, m=app.pairs(K, sex, mean), sd=app.pairs(K, sex, sd)), B=10, thread=1, unpairs=TRUE, fixed.pairs=sex)
par(mfrow=c(1,2))
hist(app_permut(K, sex, couple, 1e3, method='zscore_oppositeVSall', thread=1)$permut, breaks=30)
hist( app.permut(K, couple, group=sex, FUN=function(x) zscore(x, m=app.pairs(K, sex, mean), sd=app.pairs(K, sex, sd)), B=1e3, thread=1, unpairs=TRUE, fixed.pairs=sex)$permut, breaks=30)

app_permut(K, sex, couple, 10, method='zscore_oppositeVSall', thread=2)
app.permut(K, couple, group=sex, FUN=function(x) zscore(x, m=app.pairs(K, sex, mean), sd=app.pairs(K, sex, sd)), B=10, thread=2, unpairs=TRUE, fixed.pairs=sex)
par(mfrow=c(1,2))
hist(app_permut(K, sex, couple, 1e3, method='zscore_oppositeVSall', thread=2)$permut, breaks=30)
hist( app.permut(K, couple, group=sex, FUN=function(x) zscore(x, m=app.pairs(K, sex, mean), sd=app.pairs(K, sex, sd)), B=1e3, thread=2, unpairs=TRUE, fixed.pairs=sex)$permut, breaks=30)



source('~/Programmes/Het_haplo.R')
source('~/Mixed Model/Package R Gaston/gaston.pop/R/Het_haplo.r')
data <- read.bed.matrix('~/math_stats/Projet FIGHT HF/scratch/FIGHTHF_GSA')

i=22
map <- genmapB37$All_1000G[[paste('chr', i, sep='')]]
str(het.haplo.chr(select.snps(data, chr==i), map, block=0.5, minl=2, maxl=+Inf))
str(het_haplo_chr(select.snps(data, chr==i), map, block=0.5, minl=2, maxl=+Inf, method='haplo.stats'))

i=c(21,22)
str(het.haplo(select.snps(data, chr %in% i), 'B37', block=0.5, minl=2, maxl=+Inf, thread=1)->t)
str(het_haplo(select.snps(data, chr %in% i), 'B37', block=0.5, minl=2, maxl=+Inf, method='haplo.stats', thread=1) -> tt)
sum(t$chr21!=tt$chr21)
sum(t$chr22!=tt$chr22)

str(het.haplo(select.snps(data, chr %in% i), 'B37', block=0.5, minl=2, maxl=+Inf, thread=2)->t)
str(het_haplo(select.snps(data, chr %in% i), 'B37', block=0.5, minl=2, maxl=+Inf, method='haplo.stats', thread=2) -> tt)
sum(t$chr21!=tt$chr21)
sum(t$chr22!=tt$chr22)





x <- select.snps(read.bed.matrix("Example"), chr==22 & pos<3600000*10 & maf>0.05)
data(LCT)
x <- select.inds(x, id %in% LCT.fam$id)
x@ped$famid <- LCT.pop[LCT.fam$id %in% x@ped$id]
x <- select.snps(x, maf>0.05)

app_genome( x, relatedness='Rousset', method='mean', thread=1, build='B37' )[[22]] -> t

source('~/Mixed Model/Package R Gaston/gaston.pop/R/across_genome.r')
rouss <- function(x) Rousset(x)$Rousset
tt <- app.genome(x, as.numeric(x@ped$famid), FUN=mean, relatedness=rouss, windows=3600000, unit="base", sliding=300000, thread=1, map='B37', LD.thin=NULL)
ttt <- app.genome(x, as.numeric(x@ped$famid), FUN=mean, relatedness=rouss, windows=3600000, unit="base", sliding=300000, thread=2, map='B37', LD.thin=NULL)

sum(t$start!=tt$start)
sum(t$end!=tt$end)
sum(t$centro!=tt$centro)
sum(t$num!=tt$num, na.rm=T)
sum(t$recombi_region!=tt$recombi_region, na.rm=T)
t[which(tt$recombi_region!=t$recombi_region),]
tt[which(tt$recombi_region!=t$recombi_region),]
sum(t$recombi_all!=tt$recombi_snps, na.rm=T)
t[which(t$recombi_all!=tt$recombi_snps),]
tt[which(t$recombi_all!=tt$recombi_snps),]

sum(t$LD!=tt$LD, na.rm=T)
sum(t$LD_sd!=tt$LD_sd, na.rm=T)
sum(t$maf!=tt$maf, na.rm=T)
sum(t$maf_sd!=tt$maf_sd, na.rm=T)

sum(t$rousset_mean!=tt$rouss_mean, na.rm=T)



sum(t$start!=ttt$start)
sum(t$end!=ttt$end)
sum(t$centro!=ttt$centro)
sum(t$num!=ttt$num, na.rm=T)
sum(t$recombi_region!=ttt$recombi_region, na.rm=T)
t[which(ttt$recombi_region!=t$recombi_region),]
ttt[which(ttt$recombi_region!=t$recombi_region),]
sum(t$recombi_all!=ttt$recombi_snps, na.rm=T)
t[which(t$recombi_all!=ttt$recombi_snps),]
ttt[which(t$recombi_all!=ttt$recombi_snps),]

sum(t$LD!=ttt$LD, na.rm=T)
sum(t$LD_sd!=ttt$LD_sd, na.rm=T)
sum(t$maf!=ttt$maf, na.rm=T)
sum(t$maf_sd!=ttt$maf_sd, na.rm=T)

sum(t$rousset_mean!=ttt$rouss_mean, na.rm=T)

