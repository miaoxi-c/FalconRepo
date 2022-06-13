## Not run:
#Simulate 3 species data
N <- 1000
nspecies <- 3
J <- 5
occ_covs <- as.data.frame(matrix(rnorm(N * 10),ncol=10))
names(occ_covs) <- paste('occ_cov',1:10,sep='')
det_covs <- list()
for (i in 1:nspecies){
  det_covs[[i]] <- matrix(rnorm(N*J),nrow=N)
}
names(det_covs) <- paste('det_cov',1:nspecies,sep='')
#True vals
beta <- c(0.5,0.2,0.4,0.5,-0.1,-0.3,0.2,0.1,-1,0.1)
f1 <- beta[1] + beta[2]*occ_covs$occ_cov1
f2 <- beta[3] + beta[4]*occ_covs$occ_cov2
f3 <- beta[5] + beta[6]*occ_covs$occ_cov3
f4 <- beta[7]
f5 <- beta[8]
f6 <- beta[9]
f7 <- beta[10]
f <- cbind(f1,f2,f3,f4,f5,f6,f7)
z <- expand.grid(rep(list(1:0),nspecies))[,nspecies:1]
colnames(z) <- paste('sp',1:nspecies,sep='')
dm <- model.matrix(as.formula(paste0("~.^",nspecies,"-1")),z)
psi <- exp(f %*% t(dm))
psi <- psi/rowSums(psi)
#True state
ztruth <- matrix(NA,nrow=N,ncol=nspecies)
for (i in 1:N){
  ztruth[i,] <- as.matrix(z[sample(8,1,prob=psi[i,]),])
}
p_true <- c(0.6,0.7,0.5)
# fake y data
y <- list()
for (i in 1:nspecies){
  y[[i]] <- matrix(NA,N,J)
  for (j in 1:N){
    for (k in 1:J){
      y[[i]][j,k] <- rbinom(1,1,ztruth[j,i]*p_true[i])
    } 
  }
}
names(y) <- c('coyote','tiger','bear')
#Create the unmarked data object
data = unmarkedFrameOccuMulti(y=y,siteCovs=occ_covs,obsCovs=det_covs)
#Summary of data object
summary(data)
plot(data)
# Look at f parameter design matrix
data@fDesign
# Formulas for state and detection processes
# Length should match number/order of columns in fDesign
occFormulas <- c('~occ_cov1','~occ_cov2','~occ_cov3','~1','~1','~1','~1')
#Length should match number/order of species in data@ylist
detFormulas <- c('~1','~1','~1')
fit <- occuMulti(detFormulas,occFormulas,data)
#Look at output
fit
plot(fit)
#Compare with known values
cbind(c(beta,log(p_true/(1-p_true))),fit@opt$par)
#predict method
lapply(predict(fit,'state'),head)
lapply(predict(fit,'det'),head)
#marginal occupancy
head(predict(fit,'state',species=2))
head(predict(fit,'state',species='bear'))
head(predict(fit,'det',species='coyote'))
#probability of co-occurrence of two or more species
head(predict(fit, 'state', species=c('coyote','tiger')))
#conditional occupancy
head(predict(fit,'state',species=2,cond=3)) #tiger | bear present
head(predict(fit,'state',species='tiger',cond='bear')) #tiger | bear present
head(predict(fit,'state',species='tiger',cond='-bear')) #bear absent
head(predict(fit,'state',species='tiger',cond=c('coyote','-bear')))
#residuals (by species)
lapply(residuals(fit),head)
#ranef (by species)
ranef(fit, species='coyote')
#parametric bootstrap
bt <- parboot(fit,nsim=30)
#update model
occFormulas <- c('~occ_cov1','~occ_cov2','~occ_cov2+occ_cov3','~1','~1','~1','~1')
fit2 <- update(fit,stateformulas=occFormulas)
#List of fitted models
fl <- fitList(fit,fit2)
coef(fl)
#Model selection
modSel(fl)
#Fit model while forcing some natural parameters to be 0
#For example: fit model with no species interactions
occFormulas <- c('~occ_cov1','~occ_cov2','~occ_cov2+occ_cov3','0','0','0','0')
fit3 <- occuMulti(detFormulas,occFormulas,data)
#Alternatively, you can force all interaction parameters above a certain
#order to be zero with maxOrder. This will be faster.
occFormulas <- c('~occ_cov1','~occ_cov2','~occ_cov2+occ_cov3')
fit4 <- occuMulti(detFormulas,occFormulas,data,maxOrder=1)
#Add Bayes penalty term to likelihood. This is useful if your parameter
#estimates are very large, eg because of separation.
fit5 <- occuMulti(detFormulas, occFormulas, data, penalty=1)
#Find optimal penalty term value from a range of possible values using
#K-fold cross validation, and re-fit the model
fit_opt <- optimizePenalty(fit5, penalties=c(0,1,2))
## End(Not run)

