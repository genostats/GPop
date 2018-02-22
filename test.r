source('Programmes/mating.r')
source('Mixed Model/Package R Gaston/gaston.pop/R/app_pairs.r')
require(gaston)

n <- 50
K <- random.pm(n=100)$K
sex <- rep(c(1,2), n)
couple <- rep(1:n, each=2)

mean_couple(K, couple)
app_pairs(K, couple, mean)

sum( app_nonmates(K, couple)!= app_unpairs(K, couple) )

sum( app_pairs(K, couple)!=app_couple(K, couple) )

sum( app_unpairs(K, list(couple, sex))!=app_couple(K, couple, sex, opposite=T) )

z_score(K, couple, sex, method='all')
zscore(app_pairs(K, couple), m=app_unpairs(K, couple, mean), sd=app_unpairs(K, couple, sd))

z_score(K, couple, sex, method='opp_comp')
zscore(app_pairs(K, couple), m=app_unpairs(K, list(couple, sex), mean), sd=app_unpairs(K, list(couple, sex), sd))

z_score(K, couple, sex, method='opp_int')
zscore(app_unpairs(K, list(couple, sex)), m=app_pairs(K, sex, mean), sd=app_pairs(K, sex, sd))



source('Mixed Model/Package R Gaston/gaston.pop/R/permut_pairs.r')

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
app.permut(K, couple, group=sex, function(x) zscore(x, m=app_unpairs(K, couple, mean), sd=app_unpairs(K, couple, sd)), unpairs=FALSE, fixed.pairs=NULL, B=10, thread=1)
par(mfrow=c(1,2))
hist(app_permut(K, sex, couple, 1e3, method='zscore', thread=1)$permut, breaks=30)
hist( app.permut(K, couple, group=sex, function(x) zscore(x, m=app_unpairs(K, couple, mean), sd=app_unpairs(K, couple, sd)), unpairs=FALSE, fixed.pairs=NULL, B=1e3, thread=1)$permut, breaks=30)

app_permut(K, sex, couple, 10, method='zscore', thread=2)
app.permut(K, couple, group=sex, function(x) zscore(x, m=app_unpairs(K, couple, mean), sd=app_unpairs(K, couple, sd)), unpairs=FALSE, fixed.pairs=NULL, B=10, thread=2)
par(mfrow=c(1,2))
hist(app_permut(K, sex, couple, 1e3, method='zscore', thread=2)$permut, breaks=30)
hist( app.permut(K, couple, group=sex, function(x) zscore(x, m=app_unpairs(K, couple, mean), sd=app_unpairs(K, couple, sd)), unpairs=FALSE, fixed.pairs=NULL, B=1e3, thread=2)$permut, breaks=30)


app_permut(K, sex, couple, 10, method='zscore_coupleVSopposite', thread=1)
app.permut(K, couple, group=sex, FUN=function(x) zscore(x, m=app_unpairs(K, list(couple, sex), mean), sd=app_unpairs(K, list(couple, sex), sd)), B=10, thread=1)
par(mfrow=c(1,2))
hist(app_permut(K, sex, couple, 1e3, method='zscore_coupleVSopposite', thread=1)$permut, breaks=30)
hist( app.permut(K, couple, group=sex, FUN=function(x) zscore(x, m=app_unpairs(K, list(couple, sex), mean), sd=app_unpairs(K, list(couple, sex), sd)), B=1e3, thread=1)$permut, breaks=30)

app_permut(K, sex, couple, 10, method='zscore_coupleVSopposite', thread=2)
app.permut(K, couple, group=sex, FUN=function(x) zscore(x, m=app_unpairs(K, list(couple, sex), mean), sd=app_unpairs(K, list(couple, sex), sd)), B=10, thread=2)
par(mfrow=c(1,2))
hist(app_permut(K, sex, couple, 1e3, method='zscore_coupleVSopposite', thread=2)$permut, breaks=30)
hist( app.permut(K, couple, group=sex, FUN=function(x) zscore(x, m=app_unpairs(K, list(couple, sex), mean), sd=app_unpairs(K, list(couple, sex), sd)), B=1e3, thread=2)$permut, breaks=30)




app_permut(K, sex, couple, 10, method='zscore_oppositeVSall', thread=1)
app.permut(K, couple, group=sex, FUN=function(x) zscore(x, m=app_unpairs(K, couple, mean), sd=app_unpairs(K, couple, sd)), B=10, thread=1, unpairs=TRUE, fixed.pairs=sex)








