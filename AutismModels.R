#Comparing models of underlying nature of autism

#Model A - there is a single autism dimension which determines social and nonsocial features

#We will simulate 3 attributes, to include language as well

require(MASS)
jpeg("modelplot.jpg", width = 450, height = 250)
par(mfrow=c(1,2)) #2 plots side by side
#simulate one group with autism
mu <- c(2,2,2) #means on the dimensions - +ve indicates impairment
sigma <-matrix(c(1,.75,.75,
                 .75,1,.75,
                 .75,.75,1),nrow=3) #dimensions correlated in autism group
N1 = 200
adat <- mvrnorm(N1,mu,sigma)

#simulate larger group from gen population
mu <- c(0,0,0)
sigma <-matrix(c(1,.3,.3,
                 .3,1,.3,
                 .3,.3,1),nrow=3) #dimensions weakly correlated in typical group
N2=600
tdat <- mvrnorm(N2,mu,sigma)

alldat <- data.frame(rbind(adat,tdat))
alldat$group<-'1'
alldat$group[1:N1] <- '2'

plot(alldat$X1,alldat$X2,col=alldat$group,pch=1,cex=.8,
     xlab = 'Social features',ylab = 'Nonsocial features',main='Model A')
abline(h=1,lty=2)
abline(v=1,lty=2)

#Model B - the different features are only weakly correlated in population as a whole; 
#the correlation in autism is an artefact of selection

bdat<- as.data.frame(mvrnorm(800,mu,sigma))
plot(bdat$V1,bdat$V2,pch=1,cex=.8,
     xlab = 'Social features',ylab = 'Nonsocial features',
     xlim = c(-2,4),ylim=c(-2,4),main='Model B')
abline(h=1,lty=2)
abline(v=1,lty=2)


dev.off()

