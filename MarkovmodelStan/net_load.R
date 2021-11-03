

# load csv file
data_R <- read.csv(file="../hole_transitions.csv", header=TRUE, sep=",")

# initial states
# notation: 1 no-use&intact, 2 no-use&holes, 3 use&intact, 4 use&holed, 5 discarded, 6 no-use&(intact|holed), 7 use&(intact|holed)
data_R$S.i[data_R$Status.i=='New' & data_R$PHI.i == 0 ] <- 0 # state N
data_R$S.i[data_R$Status.i=='NotUse' & data_R$PHI.i == 0 ] <- 1 # state S1
data_R$S.i[data_R$Status.i=='NotUse' & data_R$PHI.i > 0 ] <- 2 # state S2
data_R$S.i[(data_R$Status.i=='LastNight' | data_R$Status.i=='NotLastNight') & data_R$PHI.i == 0 ] <- 3 # state S3
data_R$S.i[(data_R$Status.i=='LastNight' | data_R$Status.i=='NotLastNight') & data_R$PHI.i > 0 ] <- 4 # state S4

# final states, 
data_R$S.f[data_R$Status.f=='NotUse' & data_R$PHI.f == 0 ] <- 1 # state S1
data_R$S.f[data_R$Status.f=='NotUse' & data_R$PHI.f > 0 ] <- 2  # state S2
data_R$S.f[(data_R$Status.f=='LastNight' | data_R$Status.f=='NotLastNight') & data_R$PHI.f == 0 ] <- 3  # state S3
data_R$S.f[(data_R$Status.f=='LastNight' | data_R$Status.f=='NotLastNight') & data_R$PHI.f > 0 ] <- 4  # state S$
# data_R$S.f[(data_R$Status.f=='LastNight' | data_R$Status.f=='NotLastNight') & (data_R$PHI.f > 0 | data_R$PHI.i > 0)] <- 4 
data_R$S.f[data_R$Status.f=='Attrited'] <- 5  # state A
data_R$S.f[data_R$Status.f=='NotUse' & is.na(data_R$PHI.f)] <- 6  # state 'S1ORS2'
data_R$S.f[(data_R$Status.f=='LastNight' | data_R$Status.f=='NotLastNight') & is.na(data_R$PHI.f)] <- 7  # state 'S3ORS4'

# check that no NAs left, should give 0
sum(is.na(data_R$S.i))
sum(is.na(data_R$S.f))

# produce matrix of counts
counts <- matrix(data=0,nrow=5,ncol=7)

for (ii in c(1:nrow(data_R))){
 counts[data_R$S.i[ii]+1,data_R$S.f[ii]] <- counts[data_R$S.i[ii]+1,data_R$S.f[ii]] + 1
}

