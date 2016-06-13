setwd('C:/Users/MariaIzabel/Desktop/MASTER/QUANTITATIVE GENOMICS/Genomic_Prediction_Models/')

dat = read.table("http://www.bayz.biz/mouse_allgeno.txt", header = T)
dat <- dat[!is.na(dat$FI),]
dat_train = dat[dat$CV1 == 1, ]
dat_test = dat[dat$CV1 == 2, ]
geno_train = dat_train[, 9:1821]
geno_test = dat_test[, 9:1821]
map = read.table("http://www.bayz.biz/mousesnps.map",header=F)

yteste <- dat_test$FI # Selecting the weigth as phenotypic measure for the test set
ytreino <- dat_train$FI # Selecting the weigth as phenotypic measure for the train set

#replacing zeros by the mean
yteste[is.na(yteste)] <- mean(yteste, na.rm = TRUE)
ytreino[is.na(ytreino)] <- mean(ytreino, na.rm = TRUE)

#Estimating the heritability of the trait based on markers
require(heritability)
geno.vector = seq(1, nrow(dat), by = 1) 

# Calculating the Genetic Relatedness Matrix (G or K):
genotypes = dat[,9:1821]
geno_centered = scale(genotypes,scale=F, center=T) # Z

# Inputting zero in the missing genotypes
geno_centered[is.na(geno_centered)] <- 0 # Z correct for missing values
mean(geno_centered) #almost zero
Z = geno_centered %*% t(geno_centered)
p = as.vector(colMeans(genotypes,na.rm=T)/2)
q = 1- p
var = sum(2*p*q)
K = Z/var
row.names(K) = geno.vector
colnames(K) = geno.vector


weight_heritability = marker_h2(dat$BW, geno.vector, covariates = dat$Sex + dat$batch, K, alpha = 0.05,
                                eps = 1e-06, max.iter = 100, fix.h2 = FALSE, h2 = 0.5)

weight_heritability$h2
FI_heritability = marker_h2(dat$FI, geno.vector, covariates = dat$Sex + dat$batch, K, alpha = 0.05,
                            eps = 1e-06, max.iter = 100, fix.h2 = FALSE, h2 = 0.5)
FI_heritability$h2
# The real variance of the data is:
#The real variance of the data
var(dat$BW, na.rm= T)

require(BGLR)
#install.packages('rrBLUP')
require(rrBLUP) # To be able to use the mixed.solve function
nIter<-20000 ; burnIn<-2000

# It makes a 'geno_noNA' table that can be used in BGLR.
geno_train_noNA = t(geno_train) - colMeans(geno_train, na.rm = T)
geno_train_noNA[is.na(geno_train_noNA)] = 0
geno_train_noNA = t(geno_train_noNA)

# also make "noNA" version of test genotypes; it avoids that predicted GEBV fall out
# due to missing genotypes.
geno_test_noNA = t(geno_test) - colMeans(geno_test, na.rm = T)
geno_test_noNA[is.na(geno_test_noNA)] = 0
geno_test_noNA = t(geno_test_noNA)

fixedmod = model.matrix(~factor(dat_train$Sex) + factor(dat_train$batch))

print(paste('Modelo BayesB'))
ETA = list(list(X = fixedmod, model = "FIXED"), list(X = geno_train_noNA, model='BayesB'))
fmRBB <- BGLR(y=ytreino,ETA=ETA, nIter=nIter, burnIn=burnIn,verbose= FALSE, saveAt="BB")

print(paste('Modelo BayesA'))
ETA[[1]]$model<-'BayesA'
fmRBA <- BGLR(y=ytreino,ETA=ETA, nIter=nIter, burnIn=burnIn,verbose= FALSE, saveAt="BA")

print(paste('Modelo BRR'))
ETA[[1]]$model<-'BRR'
fmRBRR <- BGLR(y=ytreino,ETA=ETA, nIter=nIter, burnIn=burnIn,verbose= FALSE, saveAt="BRR")

print(paste('Modelo BayesC'))
ETA[[1]]$model<-'BayesC'
fmRBC <- BGLR(y=ytreino,ETA=ETA, nIter=nIter, burnIn=burnIn,verbose= FALSE, saveAt="BC")

print(paste('Modelo BayesLASSO'))
ETA[[1]]$model<-'BL'
fmRBL <- BGLR(y=ytreino,ETA=ETA, nIter=nIter, burnIn=burnIn, verbose= FALSE, saveAt="BL")  

print(paste('Modelo RRBlup'))
rrblup = mixed.solve(ytreino, Z=geno_train_noNA)
assign(paste('fmRR'), rrblup)

save.image("Intake.RData")
load("Intake.RData")  

###########################################################
#       PARTE 5: OBtaining GEBVs     
############################################################

GEBVBB <-geno_test_noNA%*%fmRBB$ETA[[2]]$b
GEBVBA <-geno_test_noNA%*%fmRBA$ETA[[2]]$b
GEBVBRR <-geno_test_noNA%*%fmRBRR$ETA[[2]]$b
GEBVBC <-geno_test_noNA%*%fmRBC$ETA[[2]]$b
GEBVBL <-geno_test_noNA%*%fmRBL$ETA[[2]]$b

GEBVRRBlup <- geno_test_noNA%*%rrblup$u #GEBV for individuals of test population
GEBVRRBlup2 <- geno_train_noNA%*%rrblup$u #GEBV for individuals of train population
############################################################
#       PARTE 6: Obtaining Correlations test set
############################################################

models_breeding_values = cbind(GEBVBB, GEBVBA, GEBVBRR, GEBVBC, GEBVBL, GEBVRRBlup)
colnames(models_breeding_values) = c('GEBVBB', 'GEBVBA', 'GEBVBRR', 'GEBVBC', 'GEBVBL', 'GEBVRRBlup')
models_breeding_values = as.data.frame(models_breeding_values)
cor.tst = vector()
for(i in 1:6){
  cor.tst[i] <- cor(models_breeding_values[,i], yteste)
}
cor.tst

seq_h2 = rep(FI_heritability$h2, 6)
test_cor = rbind(cor.tst, seq_h2)
test_cor = as.data.frame(t(test_cor))
test_cor$Model = c('Bayes B', 'Bayes A', 'Bayes RR', 'Bayes C', 'Bayes Lasso', 'RRBlup')

p1 <- ggplot(test_cor, aes(x = seq_h2, y = cor.tst))
p1 + geom_point(aes(color=factor(Model)), size = 4) + geom_point(colour = "grey90", size = 1.5) + ggtitle("Food Intake") +labs(x="Heritability",y="Correlation with test data") 

############################################################
#       PARTE 7: Obtaining Correlations train set
############################################################
predicted_Y = cbind(fmRBB$yHat, fmRBA$yHat, fmRBRR$yHat, fmRBC$yHat, fmRBL$yHat, GEBVRRBlup2)
colnames(predicted_Y) = c('GEBVBB', 'GEBVBA', 'GEBVBRR', 'GEBVBC', 'GEBVBL', 'GEBVRRBlup')

cor.trn = vector()
for(i in 1:6){
  cor.trn[i] <- cor(predicted_Y[,i], ytreino) ## this correlation indicates overfitting
}

train_cor = rbind(cor.trn, seq_h2)
train_cor = as.data.frame(t(train_cor))
train_cor$Model = c('Bayes B', 'Bayes A', 'Bayes RR', 'Bayes C', 'Bayes Lasso', 'RRBlup')

p2 <- ggplot(train_cor, aes(x = seq_h2, y = cor.trn))
p2 + geom_point(aes(color=factor(Model)), size = 4) + geom_point(colour = "grey90", size = 1.5) + ggtitle("Food Intake") +labs(x="Heritability",y="Correlation with train data")  
###############################################


par(mfrow=c(2,3))
plot(fmRBB$yHat, ytreino, ylab="Food Intake",
     xlab="Estimated phenotype", cex=.8, bty="L", main= 'Bayes B')

plot(fmRBA$yHat, ytreino, ylab="Food Intake",
     xlab="Estimated phenotype", cex=.8, bty="L", main= 'Bayes A')

plot(fmRBRR$yHat, ytreino, ylab="Food Intake",
     xlab="Estimated phenotype", cex=.8, bty="L", main= 'Bayes RR')

plot(fmRBC$yHat, ytreino, ylab="Food Intake", xlab="Estimated phenotype", cex=.8, bty="L", main= 'Bayes C')

plot(fmRBL$yHat, ytreino, ylab="Food Intake",
     xlab="Estimated phenotype", cex=.8, bty="L", main= 'Bayes Lasso')

dev.off()
########### test data
par(mfrow=c(2,3))
plot(GEBVBB, yteste, ylab="Food Intake",
     xlab="Pred. Gen. Value", cex=.8, bty="L", main= 'Bayes B')

plot(GEBVBA, yteste, ylab="Food Intake",
     xlab="Pred. Gen. Value", cex=.8, bty="L", main= 'Bayes A')

plot(GEBVBRR, yteste, ylab="Food Intake",
     xlab="Pred. Gen. Value", cex=.8, bty="L", main= 'Bayes RR')

plot(GEBVBC, yteste, ylab="Food Intake",
     xlab="Pred. Gen. Value", cex=.8, bty="L", main= 'Bayes C')

plot(GEBVBL, yteste, ylab="Food Intake",
     xlab="Pred. Gen. Value", cex=.8, bty="L", main= 'Bayes Lasso')

plot(GEBVRRBlup, yteste, ylab="Food Intake",
     xlab="Pred. Gen. Value", cex=.8, bty="L", main= 'RRBlup')

dev.off()
###############################################
#             Gráficos BGLR                   #
###############################################
par(mfrow=c(2,3))

bhat <- fmRBB$ETA[[2]]$b
SD.bhat <- fmRBB$ETA[[2]]$SD.b
plot(bhat^2, ylab='Estimated Squared-Marker Effect', type='o', cex=.5, col=4, main='Bayes B')

bhat <- fmRBA$ETA[[2]]$b
SD.bhat <- fmRBA$ETA[[2]]$SD.b
plot(bhat^2, ylab='', type='o', cex=.5, col=4, main='Bayes A')

bhat <- fmRBC$ETA[[2]]$b
SD.bhat <- fmRBC$ETA[[2]]$SD.b
plot(bhat^2, ylab='', type='o', cex=.5, col=4, main='Bayes C')

bhat <- fmRBL$ETA[[2]]$b
SD.bhat <- fmRBL$ETA[[2]]$SD.b
plot(bhat^2, ylab='Estimated Squared-Marker Effect', type='o', cex=.5, col=4, main='Bayes Lasso')

bhat <- fmRBRR$ETA[[2]]$b
SD.bhat <- fmRBRR$ETA[[2]]$SD.b
plot(bhat^2, ylab='', type='o', cex=.5, col=4, main='Bayes RR')

dev.off()

###################################
# Estimating Genomic prediction usinf FI
######################################


