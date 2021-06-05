
library(psych)
library(nFactors)
library(vegan)
library(mice)
library(RVAideMemoire)
library(lavaan)

set.seed(200)

#====================================
# Collect data for autistic and SCT groups
#====================================
readosf <- 0 
if(readosf==1){
dataAutSCT<-read.csv("https://osf.io/dhu3j/download",stringsAsFactors = F)

vocabdata<-read.csv("https://osf.io/7fyre/download",stringsAsFactors = F)
grammardata<-read.csv("https://osf.io/c3bq8/download",stringsAsFactors = F)
impdata<-read.csv("https://osf.io/grafk/download",stringsAsFactors = F)
infdata<-read.csv("https://osf.io/5hr2e/download",stringsAsFactors = F)
overturesdata<-read.csv("https://osf.io/2f8ga/download",stringsAsFactors = F)
matricesdata<-read.csv("https://osf.io/xqjs8/download",stringsAsFactors = F)

CCCdata<-read.csv("https://osf.io/9rfba/download",stringsAsFactors = F)
PSCdata<-read.csv("https://osf.io//zyph6//download",stringsAsFactors = F)
threedidata<-read.csv("https://osf.io/teg8z/download",stringsAsFactors = F)
}
if(readosf ==0){
  dataAutSCT<-read.csv("participant.information.csv",stringsAsFactors = F)
  
  vocabdata<-read.csv("vocab.csv",stringsAsFactors = F)
  grammardata<-read.csv("grammar.csv",stringsAsFactors = F)
  impdata<-read.csv("implicature.csv",stringsAsFactors = F)
  infdata<-read.csv("inference.csv",stringsAsFactors = F)
  overturesdata<-read.csv("overtures.csv",stringsAsFactors = F)
  matricesdata<-read.csv("matrices.csv",stringsAsFactors = F)
  
  CCCdata<-read.csv("cccquestionnaire.csv",stringsAsFactors = F)
  PSCdata<-read.csv("psc17.csv",stringsAsFactors = F)
  threedidata<-read.csv("3di.csv",stringsAsFactors = F)
}

dataAutSCT$vocab<-NA
dataAutSCT$grammar<-NA
dataAutSCT$imp<-NA
dataAutSCT$inf<-NA
dataAutSCT$overtures<-NA
dataAutSCT$matrices<-NA
dataAutSCT$cccStruct<-NA
dataAutSCT$pscTotal<-NA
dataAutSCT$threediSocial<-NA
dataAutSCT$threediComm<-NA
dataAutSCT$threediRRBI<-NA
dataAutSCT$threediAUT<-NA

for(i in 1:length(dataAutSCT$ID)){
  myID<-dataAutSCT$ID[i]
  
  if(myID %in% vocabdata$ID){
    index<-which(myID %in% vocabdata$ID)
    dataAutSCT$vocab[i]<-vocabdata$total[which(vocabdata$ID==myID[index])]
  }
  if(myID %in% grammardata$ID){
    index<-which(myID %in% grammardata$ID)
    dataAutSCT$grammar[i]<-grammardata$total[which(grammardata$ID==myID[index])]
  }
  if(myID %in% impdata$ID){
    index<-which(myID %in% impdata$ID)
    dataAutSCT$imp[i]<-impdata$implicature_total[which(impdata$ID==myID[index])]
  }
  if(myID %in% infdata$ID){
    index<-which(myID %in% infdata$ID)
    dataAutSCT$inf[i]<-infdata$total[which(infdata$ID==myID[index])]
  }
  if(myID %in% overturesdata$ID){
    index<-which(myID %in% overturesdata$ID)
    dataAutSCT$overtures[i]<-overturesdata$total[which(overturesdata$ID==myID[index])]
  }
  if(myID %in% matricesdata$ID){
    index<-which(myID %in% matricesdata$ID)
    dataAutSCT$matrices[i]<-matricesdata$total[which(matricesdata$ID==myID[index])]
  }
  if(myID %in% CCCdata$ID){
    index<-which(myID %in% CCCdata$ID)
    dataAutSCT$cccStruct[i]<-CCCdata$struct.corrected[which(CCCdata$ID==myID[index])]
  }
  if(myID %in% PSCdata$ID){
    index<-which(myID %in% PSCdata$ID)
    dataAutSCT$pscTotal[i]<-PSCdata$PSC17_total[which(PSCdata$ID==myID[index])]
  }
  if(myID %in% threedidata$ID){
    index<-which(myID %in% threedidata$ID)
    dataAutSCT$threediSocial[i]<-threedidata$Social[which(threedidata$ID==myID[index])]
    dataAutSCT$threediComm[i]<-threedidata$Communication[which(threedidata$ID==myID[index])]
    dataAutSCT$threediRRBI[i]<-threedidata$RRBIs[which(threedidata$ID==myID[index])]
    dataAutSCT$threediAUT[i]<-threedidata$ThreeDi.autism.diagnosis[which(threedidata$ID==myID[index])]
  }
}

#drop 7 participants for whom there is no 3di data
dataAutSCT3di<-dataAutSCT[!is.na(dataAutSCT$threediAUT),]

describeBy(dataAutSCT3di,group=dataAutSCT3di$trisomy.diagnosis)
describeBy(dataAutSCT3di,group=dataAutSCT3di$prenatal.trisomy.diagnosis)


#
#====================================
# DB added: convert test scores to age-scaled z scores
#====================================
coeffs <- read.csv('reg_coeffs_lang.csv')
dataAutSCT3di$year<-dataAutSCT3di$UK.school.year.group
w<-which(dataAutSCT3di$year>8)
dataAutSCT3di$year[w]<-8#have to use year 8 norms for older cases

mycols<- c('vocab','grammar','imp','inf','overtures','matrices')
for (i in 1:length(mycols)){
  c<-which(colnames(dataAutSCT3di) ==mycols[i])
   nucol<-length(colnames(dataAutSCT3di))+1 #add column for z-score
  dataAutSCT3di[,nucol]<-NA #initialise
  for (j in 1:nrow(dataAutSCT3di)){
    a <- coeffs[i,2]
    b <- coeffs[i,3]
    se <- coeffs[i,4]
    x <- dataAutSCT3di$year[j]
    y <- dataAutSCT3di[j,c]
    if(is.numeric(y)){
    #transform where necessary
    if (i==2) {    #grammar
      y <- log(50-y)}
      if (i==3) {    #implicature
        y <- log(34-y)}
    if (i==4) {    #inference
      y <- log(21-y)}

    
    predscore <- a+b*x
    zscore<- (y-predscore)/se
    if (i %in% c(2:4)){
      zscore <-(predscore-y)/se #reverse for those where log transform error used
    }
    dataAutSCT3di[j,nucol]<-round(zscore,2)
    }
  }
                
}
colnames(dataAutSCT3di)[(nucol-5):nucol]<-paste(mycols,'z',sep='.')

describeBy(dataAutSCT3di[(nucol-5):nucol],group=dataAutSCT3di$trisomy.diagnosis)
describeBy(dataAutSCT3di[c(3,14:19)],group=dataAutSCT3di$trisomy.diagnosis)

#Check all looks coherent - plot raw vs z
plot(dataAutSCT3di$vocab,dataAutSCT3di$vocab.z,col=(1+dataAutSCT3di$trisomy.diagnosis))
plot(dataAutSCT3di$grammar,dataAutSCT3di$grammar.z,col=(1+dataAutSCT3di$trisomy.diagnosis))
plot(dataAutSCT3di$imp,dataAutSCT3di$imp.z,col=(1+dataAutSCT3di$trisomy.diagnosis))
plot(dataAutSCT3di$inf,dataAutSCT3di$inf.z,col=(1+dataAutSCT3di$trisomy.diagnosis))
plot(dataAutSCT3di$overtures,dataAutSCT3di$overtures.z,col=(1+dataAutSCT3di$trisomy.diagnosis))
plot(dataAutSCT3di$matrices,dataAutSCT3di$matrices.z,col=(1+dataAutSCT3di$trisomy.diagnosis))
#====================================
# Produce a language factor score
#====================================

#extract cases with language data
dataAutSCTLang<-dataAutSCT[dataAutSCT$include.language.test.data==1,]
#impute missing data
dataAutSCTLang.imputed<-complete(mice(dataAutSCTLang[,c(14:18)]))

#extract eigenvalues for correlation matrix of language measures
#and identify number of factors to extract
ev<-eigen(cor(dataAutSCTLang.imputed)) # get eigenvalues
ap<-parallel(subject=nrow(dataAutSCTLang.imputed),var=5,rep=100,cent=.05)
nS<-nScree(x=ev$values,aparallel=ap$eigen$qevpea)
plotnScree(nS)
#indicates one factor should be extracted

#carry out exploratory factor analysis to extract language factor score
lang.efa<-fa(dataAutSCTLang.imputed,1)
dataAutSCTLang$lang.score<-lang.efa$scores


#DB:check age vs 3di measures - no association
plot(dataAutSCTLang$age,dataAutSCTLang$threediSocial)
cor(dataAutSCTLang$age,dataAutSCTLang$threediSocial)
plot(dataAutSCTLang$age,dataAutSCTLang$threediComm)
cor(dataAutSCTLang$age,dataAutSCTLang$threediComm)
plot(dataAutSCTLang$age,dataAutSCTLang$threediRRBI)
cor(dataAutSCTLang$age,dataAutSCTLang$threediRRBI)

#plot language scores against age
plot(dataAutSCTLang$age,dataAutSCTLang$lang.score)
#strong linear association, although there is an outlier (case 73),
#which is in the bottom right corner
#---------------------------------------
#DB>rather than AW correction via regression  would be better to standardize scores against norms.
# How does this compare?
#I shall include matrices as well - I doubt it is different
#---------------------------------------
dataAutSCTLangz<-dataAutSCT3di[dataAutSCT3di$include.language.test.data==1,]

w<-which(colnames(dataAutSCTLangz)=='vocab.z')
dataAutSCTLangz.imputed<-complete(mice(dataAutSCTLangz[,w:(w+5)]))

#extract eigenvalues for correlation matrix of language measures
#and identify number of factors to extract
ev<-eigen(cor(dataAutSCTLangz.imputed)) # get eigenvalues
ap<-parallel(subject=nrow(dataAutSCTLangz.imputed),var=5,rep=100,cent=.05)
nS<-nScree(x=ev$values,aparallel=ap$eigen$qevpea)
plotnScree(nS)
#indicates one factor should be extracted

#carry out exploratory factor analysis to extract language factor score
lang.efa<-fa(dataAutSCTLangz.imputed,1)
dataAutSCTLangz$lang.score<-lang.efa$scores


#DB:check age vs 3di measures - no association
plot(dataAutSCTLang$age,dataAutSCTLang$threediSocial)
cor(dataAutSCTLang$age,dataAutSCTLang$threediSocial)
plot(dataAutSCTLang$age,dataAutSCTLang$threediComm)
cor(dataAutSCTLang$age,dataAutSCTLang$threediComm)
plot(dataAutSCTLang$age,dataAutSCTLang$threediRRBI)
cor(dataAutSCTLang$age,dataAutSCTLang$threediRRBI)

#plot language scores against age
plot(dataAutSCTLangz$age,dataAutSCTLangz$lang.score)
cor(dataAutSCTLangz$age,dataAutSCTLangz$lang.score)



#set up regression to remove effect of age on language.(original factor score)
#exclude case 73 when determining the regression coefficent for the age adjustment
age.adjustment.lang<-lm(lang.score ~ age, data=dataAutSCTLang[-73,])
#now remove the prediction for age from the language score to control for age
dataAutSCTLang$lang.score.age.adjusted<-dataAutSCTLang$lang.score-
  predict(age.adjustment.lang,newdata=dataAutSCTLang)
dataAutSCTLang$lang.score.age.adjusted<-(dataAutSCTLang$lang.score.age.adjusted-mean(dataAutSCTLang$lang.score.age.adjusted))/sd(dataAutSCTLang$lang.score.age.adjusted)

#carry out a similar adjustment for age on nonverbal reasoning
plot(dataAutSCTLang$age,dataAutSCTLang$matrices)
age.adjustment.NV<-lm(matrices ~ age, data=dataAutSCTLang)
dataAutSCTLang$matrices.age.adjusted<-dataAutSCTLang$matrices-
  predict(age.adjustment.NV,newdata=dataAutSCTLang)
dataAutSCTLang$matrices.age.adjusted<-(dataAutSCTLang$matrices.age.adjusted-mean(na.omit(dataAutSCTLang$matrices.age.adjusted)))/sd(na.omit(dataAutSCTLang$matrices.age.adjusted))

#====================================
# Assess factorial invariance of the 3di
# across autistic and SCT groups
#====================================

model.3di<-'
autism=~threediComm+threediSocial+threediRRBI
'
noconstraints<-cfa(model=model.3di,estimator="MLM",data=dataAutSCT,group="trisomy.diagnosis")
fit1<-t(data.frame(fitMeasures(noconstraints,c("cfi.robust","rmsea.robust","rmsea.ci.lower.robust","rmsea.ci.upper.robust"))))
loadings<-cfa(model=model.3di,estimator="MLM",data=dataAutSCT,group="trisomy.diagnosis",group.equal=c("loadings"))
fit2<-t(data.frame(fitMeasures(loadings,c("cfi.robust","rmsea.robust","rmsea.ci.lower.robust","rmsea.ci.upper.robust"))))
anova(noconstraints,loadings)

intercepts<-cfa(model=model.3di,data=dataAutSCT,estimator="MLM",group="trisomy.diagnosis",group.equal=c("loadings","intercepts"))
fit3<-t(data.frame(fitMeasures(intercepts,c("cfi.robust","rmsea.robust","rmsea.ci.lower.robust","rmsea.ci.upper.robust"))))
anova(loadings,intercepts)

intercepts.varyComm<-cfa(model=model.3di,estimator="MLM",data=dataAutSCT,group="trisomy.diagnosis",group.equal=c("loadings","intercepts"),group.partial=c("threediComm~1"))
fit4<-t(data.frame(fitMeasures(intercepts.varyComm,c("cfi.robust","rmsea.robust","rmsea.ci.lower.robust","rmsea.ci.upper.robust"))))
intercepts.varySocial<-cfa(model=model.3di,estimator="MLM",data=dataAutSCT,group="trisomy.diagnosis",group.equal=c("loadings","intercepts"),group.partial=c("threediSocial~1"))
fit5<-t(data.frame(fitMeasures(intercepts.varySocial,c("cfi.robust","rmsea.robust","rmsea.ci.lower.robust","rmsea.ci.upper.robust"))))

intercepts.varyRRBI<-cfa(model=model.3di,estimator="MLM",data=dataAutSCT,group="trisomy.diagnosis",group.equal=c("loadings","intercepts"),group.partial=c("threediRRBI~1"))
fit4<-t(data.frame(fitMeasures(intercepts.varyRRBI,c("cfi.robust","rmsea.robust","rmsea.ci.lower.robust","rmsea.ci.upper.robust"))))
residuals<-cfa(model=model.3di,estimator="MLM",data=dataAutSCT,group="trisomy.diagnosis",group.equal=c("loadings","intercepts","residuals"),group.partial=c("threediRRBI~1"))
fit7<-t(data.frame(fitMeasures(residuals,c("cfi.robust","rmsea.robust","rmsea.ci.lower.robust","rmsea.ci.upper.robust"))))
anova(intercepts.varyRRBI,residuals)

#Analysis indicates partial factorial invariance -- the RRBI intercept
#varies across groups.
#Now set up full group factor analysis, with RRBI intercept varying by group
model.3di.adjustedRRBI<-'
autism=~threediComm+threediSocial+threediRRBI
threediRRBI+autism~trisomy.diagnosis
'
model.3di.fit<-cfa(model=model.3di.adjustedRRBI,estimator="MLM",data=dataAutSCT)
standardizedsolution(model.3di.fit)

#extract factor score for 3di
dataAutSCTLang$threedi.Fscore<-predict(model.3di.fit,newdata=dataAutSCTLang)

#====================================
# Assess relationship between language
# and 3di dimensions across groups
#====================================

manova.lang.aut.1<-adonis.II(cbind(threediComm,threediSocial,threediRRBI) ~ 
                trisomy.diagnosis*lang.score.age.adjusted, data=dataAutSCTLang)
manova.lang.aut.2<-adonis.II(cbind(threediComm,threediSocial,threediRRBI) ~ 
                trisomy.diagnosis*cccStruct, data=dataAutSCT3di)

#====================================
# Assess whether groups differed in
# level of SEN support when controlling
# for level of difficulties
#====================================

SEN.1<-lm(SEN.provision ~ lang.score.age.adjusted + age + pscTotal + threedi.Fscore,
            data=dataAutSCTLang)
SEN.2<-lm(SEN.provision ~ lang.score.age.adjusted + age + pscTotal + threedi.Fscore + trisomy.diagnosis,
       data=dataAutSCTLang)
anova(SEN.1,SEN.2)
