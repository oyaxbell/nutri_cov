#### Dataset loading ####
library(readxl); library(nephro); library(haven); library(survival); library(pROC); library(survminer); library(glmnetUtils); library(ggplotify); library(timeROC)
library(glmnet);library(mice); library(tidyverse); library(caret); library(OptimalCutpoints);library(rms);library(pROC); library(coxphw); library(MAMI)
library(FactoMineR); library(fpc); library(factoextra); library(rms); library(gridExtra); library(smoothROCtime); library (VIM); library(Hmisc)

setwd("~/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/COVID-Nutricion")
setwd("C:/Users/Usuario Lab. Datos/OneDrive - UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/COVID-Nutricion")
setwd("/Users/nefoantonio/UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO/OMAR YAXMEHEN BELLO CHAVOLLA - COVID-Nutricion")
nutri <- read_excel("covid_nutri.xlsx", na = c("", "N/A", "NA","na", "#N/A"))
names(nutri)<-gsub(" ","_", names(nutri))
nutri$defuncion[nutri$Status_actual==1]<-1;nutri$defuncion[nutri$Status_actual==0 | is.na(nutri$Status_actual)]<-0
nutri$edad40[nutri$edad<40]<-1;nutri$edad40[nutri$edad>=40]<-0
nutri$edad65[nutri$edad<65]<-0;nutri$edad65[nutri$edad>=65]<-1
nutri$tac[nutri$`RESULTADO_TAC_(ver_comentario)`=="compatible"]<-1;nutri$tac[!nutri$`RESULTADO_TAC_(ver_comentario)`=="compatible" | is.na(nutri$`RESULTADO_TAC_(ver_comentario)`)]<-0
nutri$sexo[nutri$sexo=="H"]<-1;nutri$sexo[nutri$sexo=="M" | nutri$sexo=="m"]<-0
nutri$nac<-as.numeric((nutri$disnea+nutri$tos+nutri$fiebre+nutri$malestargeneral>2)+nutri$tac>1)
table(nutri$nac)
nutri$msl<-nutri$nac*7+nutri$dm2*nutri$edad40*5+nutri$edad65*3+nutri$irc*3+
  nutri$inmunoupresion*1+nutri$epoc*1+nutri$obesidad*1+
  nutri$dm2*1+nutri$edad40*(-6)
nutri$score_msl[nutri$msl<=0]<-0
nutri$score_msl[nutri$msl>=1 & nutri$msl<=3]<-1;nutri$score_msl[nutri$msl>=4 & nutri$msl<=7]<-2
nutri$score_msl[nutri$msl>=8 & nutri$msl<=11]<-3;nutri$score_msl[nutri$msl>=12]<-4
nutri$fr_cat[nutri$fr<28 & nutri$fr>=20]<-1;nutri$fr_cat[nutri$fr<20]<-0;nutri$fr_cat[nutri$fr>=28]<-2
nutri$o2[!is.na(nutri$fio2)]<-1;nutri$o2[is.na(nutri$fio2)]<-0
nutri$sat[nutri$sao2aa<87]<-2;nutri$sat[nutri$sao2aa>=87 & nutri$sao2aa<93]<-1;nutri$sat[nutri$sao2aa>=93]<-0
nutri$uci<-nutri$UTI;nutri$uci[is.na(nutri$UTI)]<-0
nutri$uci_dias<-nutri$FECHA_DE_INGRESO_A_UTI-nutri$fechainiciosintomas
nutri$uci_dias[is.na(nutri$FECHA_DE_INGRESO_A_UTI)]<-nutri$diasestancia
nutri$uci_dias[is.na(nutri$FECHA_DE_INGRESO_A_UTI)]<-nutri$FU_time
nutri$uci_dias<-as.numeric(nutri$uci_dias)
nutri$intubado<-nutri$VMI;nutri$intubado[is.na(nutri$VMI)]<-0
nutri$defuncion[is.na(nutri$defuncion)]<-0
nutri$critical[(nutri$defuncion+nutri$uci+nutri$intubado)>0]<-1;nutri$critical[(nutri$defuncion+nutri$uci+nutri$intubado)==0]<-0
nutri$fechavaloracioninicio<-as.Date(nutri$fechavaloracioninicio)
nutri$FU_time[nutri$defuncion==0 & critical==1]<-nutri$uci_dias
nutri$FU_time[nutri$FU_time<5 & nutri$critical==0]<-7
nutri$n_comorb<-nutri$obesidad+nutri$hipertension+nutri$dm2+nutri$asma+nutri$epoc+nutri$vih+
  nutri$inmunoupresion+nutri$cardiovascular+nutri$irc+nutri$tabaquismo+nutri$insuficiencia_hepatica
nutri$ecgme0sde15<-na.tools::na.replace(nutri$ecgme0sde15,0)
nutri$comorb<-ifelse(nutri$n_comorb>0, 1, 0);nutri$comorb<-na.tools::na.replace(nutri$comorb, 0)
nutri<-nutri %>% filter(!is.na(msl))

## Función para Cox imputado
cox.imp<-function(x){
  var<-summary(pool(x))[1]
  HR<-exp(summary(pool(x))[2])
  lwr<-exp(summary(pool(x))[2]-1.96*summary(pool(x))[3])
  upr<-exp(summary(pool(x))[2]+1.96*summary(pool(x))[3])
  pval<-summary(pool(x))[6]
  print(cbind(var,HR, lwr, upr, pval))
}

#### Multiple imputation ####
nutri2<- nutri %>% filter(!is.na(msl)) %>% dplyr::select(1,4,24:25,5:6,9,28:44,48:51,54:55,211,64,206,214:219, 221)
nutri2$fio2[nutri2$pin %in% nutri2$pin]<-nutri$fio2
nutri2$Charlson[nutri2$pin %in% nutri2$pin]<-nutri$Charlson
nutri2$o2[nutri2$pin %in% nutri2$pin]<-nutri$o2
nutri2$qsofa[nutri2$pin %in% nutri2$pin]<-nutri$qsofa
nutri2$news[nutri2$pin %in% nutri2$pin]<-nutri$news
nutri2$o2[!is.na(nutri2$fio2)]<-1;nutri2$o2[is.na(nutri2$fio2)]<-0
nutri2$fio2[is.na(nutri2$fio2)]<-0.21
nutri2$sexo[is.na(nutri2$sexo)]<-0
nutri2$hospitalizacion<-ifelse(nutri2$hospitalizacion==1, 1, 0)
nutri2$FU_time<-nutri2$FU_time+1
nutri2$pin[is.na(nutri2$pin)]<-c("cov1", "cov2")
aggr_plot <- aggr(nutri_score1, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE, 
                  labels=names(aggr), cex.axis=.7, gap=3, ylab=c("Histogram of missing data","Pattern"))
set.seed(123);nutri_imp1<-mice::mice(nutri2, m=5, maxit=5)
nutri_score1<-complete(nutri_imp1, "long", include = T)
nutri_score1$score_msl[nutri_score1$msl<=0]<-0
nutri_score1$score_msl[nutri_score1$msl>=1 & nutri_score1$msl<=3]<-1;nutri_score1$score_msl[nutri_score1$msl>=4 & nutri_score1$msl<=7]<-2
nutri_score1$score_msl[nutri_score1$msl>=8 & nutri_score1$msl<=11]<-3;nutri_score1$score_msl[nutri_score1$msl>=12]<-4
nutri_score1$fr_cat[nutri_score1$fr>=25 & nutri_score1$fr<30]<-1;nutri_score1$fr_cat[nutri_score1$fr<25]<-0; nutri_score1$fr_cat[nutri_score1$fr>=30]<-2
nutri_score1$sat[nutri_score1$sao2aa<85]<-2;nutri_score1$sat[nutri_score1$sao2aa>=85 & nutri_score1$sao2aa<92]<-1;nutri_score1$sat[nutri_score1$sao2aa>=92]<-0
nutri_imp1<-as.mids(nutri_score1)
nutri0<-complete(nutri_imp1, 1)

#### Preliminary model selection ####
nutri_score1<-complete(nutri_imp1, "long", include = T)
nutri3_pre <- nutri_score1[nutri_score1$fechavaloracioninicio<="2020-06-04",]
nutri_imp3<-as.mids(nutri3_pre)

set.seed(123);m1 <- mami(nutri_imp3,
                         model="cox", outcome=c("FU_time","critical"),
                         method="LASSO", alpha=1, kfold=10,report.exp=FALSE,
                         var.remove=c("pin", "hospitalizacion", "diasdesdeiniciosintomasavaloracion",
                                      "uci", "uci_dias", "intubado", "critical", "comorb", "fio2", "Charlson",
                                      "qsofa", "news", "defuncion", "sat", "fr_cat", "score_msl", "fechavaloracioninicio", "edad"))
k1<-summary(m1)

set.seed(123);m1 <- mami(nutri_imp3,
                         model="cox", outcome=c("FU_time","critical"),
                         method="LASSO", alpha=1, kfold=10,report.exp=FALSE,
                         var.remove=c("pin", "hospitalizacion", "diasdesdeiniciosintomasavaloracion",
                                      "uci", "uci_dias", "intubado", "critical", "comorb", "fio2", "Charlson",
                                      "qsofa", "news", "defuncion", "sat", "fr_cat", "score_msl", "fechavaloracioninicio",
                                      "fiebre", "imc", "tos", "diarrea", "odinofagia", "malestargeneral", "rinorrea", "vomito",
                                      "conjuntivitis", "fc", "tad", "edad"))
summary(m1)


### Respiratory Rate, SpO2 and MSL-COVID-19 remain with significant with LASSO CIs

#md.pattern(nutri_score1)
#densityplot(nutri_imp1)

#stripplot(tempData, pch = 20, cex = 1.2)
#### Severity models ####
m0<-with(nutri_imp1, coxph(Surv(FU_time, critical)~edad))
cox.imp(m0)

m0<-with(nutri_imp1, coxph(Surv(FU_time, critical)~diasdesdeiniciosintomasavaloracion))
cox.imp(m0)

m0<-with(nutri_imp1,coxph(Surv(FU_time, critical)~sexo))
cox.imp(m0)

m0<-with(nutri_imp1,coxph(Surv(FU_time, critical)~imc))
cox.imp(m0)

m0<-with(nutri_imp1,coxph(Surv(FU_time, critical)~rcs(sao2aa,3)))
cox.imp(m0)

m0<-with(nutri_imp1,coxph(Surv(FU_time, critical)~dm2))
cox.imp(m0)

m0<-with(nutri_imp1,coxph(Surv(FU_time, critical)~obesidad))
cox.imp(m0)

m0<-with(nutri_imp1,coxph(Surv(FU_time, critical)~Charlson))
cox.imp(m0)

m0<-with(nutri_imp1,coxph(Surv(FU_time, critical)~hipertension))
cox.imp(m0)

m0<-with(nutri_imp1,coxph(Surv(FU_time, critical)~cefalea))
cox.imp(m0)

m0<-with(nutri_imp1,coxph(Surv(FU_time, critical)~disnea))
cox.imp(m0)

m0<-with(nutri_imp1,coxph(Surv(FU_time, critical)~rcs(fr,3)))
cox.imp(m0)

m0<-with(nutri_imp1,coxph(Surv(FU_time, critical)~tad))
cox.imp(m0)

### Combinacion ###
library(simPH)
m1<-with(nutri_imp1,coxph(Surv(FU_time, critical)~msl+poly(sao2aa,3)*fr))
cox.imp(m1)

m1<-coxph(Surv(FU_time, critical)~msl, data=nutri0)
summary(m1)
cbind(m1$concordance[6],m1$concordance[6]-(m1$concordance[7]*1.96),
      m1$concordance[6]+(m1$concordance[7]*1.96))

m1<-coxph(Surv(FU_time, critical)~edad+sao2aa+fr, data=nutri0)
summary(m1)
cbind(m1$concordance[6],m1$concordance[6]-(m1$concordance[7]*1.96),
      m1$concordance[6]+(m1$concordance[7]*1.96))

m1<-coxph(Surv(FU_time, critical)~msl+sao2aa+fr, data=nutri0)
summary(m1)
g1<-ggcoxdiagnostics(m1, type = "score")
g2<-ggcoxdiagnostics(m1)
g3<-ggarrange(g1, g2)
m1<-coxph(Surv(FU_time, critical)~pspline(sao2aa,2), data=nutri0)
summary(m1)

library(Greg)
m0 <- cph(Surv(FU_time, critical) ~ rcs(sao2aa, 3), data = nutri0, x = TRUE, y = TRUE)
BIC(m0)
plotHR(m0,
       term = "sao2aa", 
       xlab = "SpO2 (%)",
       col.term = c("#08519C", "#77777799"),
       lty.term = c(1, 2),
       plot.bty = "l", xlim = c(0,100))
m1 <- cph(Surv(FU_time, critical) ~ rcs(score_crit,3), data = nutri0, x = TRUE, y = TRUE)
BIC(m1)
plotHR(m1,term = "score_crit", 
       xlab = "Respiratory rate (bpm)",
       col.term = c("#08519C", "#77777799"),
       lty.term = c(1, 2),
       plot.bty = "l", xlim = c(0,40))


#### Sample size ####

## Assuming 15% of cases develop severe COVID-19 during a mean of 15 days of follow-up
## Event rate of 10 cases per 1,000 person-days and Cox-Snell of 0.167
## Modeling up to 36 parameters

pmsampsize::pmsampsize(type = "s", rsquared = 0.167, parameters = 35, rate = 0.01,
               timepoint = 15, meanfup = 11.72)

#### Validation MSL ####
### Mortality ###
m1<-coxph(Surv(FU_time, defuncion)~msl, data=nutri2)
summary(m1);cox.zph(m1)
m1<-ggcoxdiagnostics(m1)
cbind(m1$concordance[6],m1$concordance[6]-(m1$concordance[7]*1.96),
      m1$concordance[6]+(m1$concordance[7]*1.96))

set.seed(123)
m1<-cph(Surv(FU_time, defuncion)~msl, data=nutri0, x=T, y=T)
rms::validate(fit = m1, method="boot", B=1000)

### UCI ###
m1<-coxph(Surv(uci_dias, uci)~msl, data=nutri0)
summary(m1)
cbind(m1$concordance[6],m1$concordance[6]-(m1$concordance[7]*1.96),
      m1$concordance[6]+(m1$concordance[7]*1.96))


set.seed(123)
m1<-cph(Surv(uci_dias, uci)~msl, data=nutri0, x=T, y=T)
rms::validate(fit = m1, method="boot", B=1000)

### Intubacion ###
m1<-coxph(Surv(uci_dias, intubado)~msl, data=nutri0)
summary(m1)
cbind(m1$concordance[6],m1$concordance[6]-(m1$concordance[7]*1.96),
      m1$concordance[6]+(m1$concordance[7]*1.96))

set.seed(123)
m1<-cph(Surv(uci_dias, intubado)~msl, data=nutri0, x=T, y=T)
rms::validate(fit = m1, method="boot", B=1000)

### Severe disease ###
m1<-coxph(Surv(FU_time, critical)~msl, data=nutri0)
summary(m1)
cbind(m1$concordance[6],m1$concordance[6]-(m1$concordance[7]*1.96),
      m1$concordance[6]+(m1$concordance[7]*1.96))

set.seed(123)
m1<-cph(Surv(FU_time, critical)~msl, data=nutri0, x=T, y=T)
rms::validate(fit = m1, method="boot", B=1000)

set.seed(123)
m1<-cph(Surv(FU_time, critical)~msl+fr+sao2aa+sexo, data=nutri2, x=T, y=T)
rms::validate(fit = m1, method="boot", B=1000)

m1<-coxph(Surv(FU_time, critical)~msl+fr+sao2aa+sexo, data=nutri2)
cbind(m1$concordance[6],m1$concordance[6]-(m1$concordance[7]*1.96),
      m1$concordance[6]+(m1$concordance[7]*1.96))

#### Nutri-CoV ####
nutri_score<-nutri2 %>% filter(!is.na(msl))%>%dplyr::select(pin, fechavaloracioninicio,FU_time, critical,edad,sexo,fr,sao2aa,msl)
nutri3 <- nutri_score[nutri_score$fechavaloracioninicio<="2020-06-04",]
set.seed(123);nutri3_imp<-mice(nutri3, m=5, maxit=5)
nutri3_1<-complete(nutri3_imp, "long", include = TRUE)
nutri3_1$score_msl[nutri3_1$msl<=0]<-0
nutri3_1$score_msl[nutri3_1$msl>=1 & nutri3_1$msl<=3]<-1;nutri3_1$score_msl[nutri3_1$msl>=4 & nutri3_1$msl<=7]<-2
nutri3_1$score_msl[nutri3_1$msl>=8 & nutri3_1$msl<=11]<-3;nutri3_1$score_msl[nutri3_1$msl>=12]<-4
nutri3_1$fr_cat[nutri3_1$fr>=24 & nutri3_1$fr<30]<-1;nutri3_1$fr_cat[nutri3_1$fr<24]<-0; nutri3_1$fr_cat[nutri3_1$fr>=30]<-2
nutri3_1$sat[nutri3_1$sao2aa<85]<-2;nutri3_1$sat[nutri3_1$sao2aa>=85 & nutri3_1$sao2aa<92]<-1;nutri3_1$sat[nutri3_1$sao2aa>=92]<-0
nutri3_imp<-as.mids(nutri3_1)

nutri3$score_msl[nutri3$msl<=0]<-0
nutri3$score_msl[nutri3$msl>=1 & nutri3$msl<=3]<-1;nutri3$score_msl[nutri3$msl>=4 & nutri3$msl<=7]<-2
nutri3$score_msl[nutri3$msl>=8 & nutri3$msl<=11]<-3;nutri3$score_msl[nutri3$msl>=12]<-4
nutri3$fr_cat[nutri3$fr>=24 & nutri3$fr<30]<-1;nutri3$fr_cat[nutri3$fr<24]<-0; nutri3$fr_cat[nutri3$fr>=30]<-2
nutri3$sat[nutri3$sao2aa<85]<-2;nutri3$sat[nutri3$sao2aa>=85 & nutri3$sao2aa<92]<-1;nutri3$sat[nutri3$sao2aa>=92]<-0

### For k=5 ###

dfs <- lapply(1:5, function(i) complete(nutri3_imp, action = i))
x<-NULL;y<-NULL
for (i in 1:5) {
  x[[i]] <- as.matrix(dfs[[i]][, c("sat","fr_cat","score_msl")])
  y[[i]] <- Surv(dfs[[i]]$FU_time, dfs[[i]]$critical)
}

alpha0<-NULL;mod_cva<-NULL;alpha<-NULL
set.seed(123)
for (i in 1:5) {
  ALPHA <- seq(0,1, by=0.01)
  mod_cva[[i]] <- cva.glmnet(x[[i]], y[[i]],family="cox",alpha = ALPHA, k=5)
  alpha[[i]] <- ALPHA[which.min(sapply(mod_cva[[i]]$modlist, function(mod) min(mod$cvm)))];alpha
}

alpha0<-mean(alpha)

set.seed(123);m1 <- mami(nutri3_imp,
                         model="cox", outcome=c("FU_time","critical"),
                         method="LASSO", alpha=alpha0, kfold=5,report.exp=FALSE,
                         add.factor = c("sat", "fr_cat"), B=100,
                         var.remove=c("fr","msl", "sao2aa", "pin", "edad",
                                      "sexo","fechavaloracioninicio"))
summary(m1)

coef<-round(m1$coefficients.s[,1]/min(abs(m1$coefficients.s[,1])))
coef

### For k=10 ###

dfs <- lapply(1:5, function(i) complete(nutri3_imp, action = i))
x<-NULL;y<-NULL
for (i in 1:5) {
  x[[i]] <- as.matrix(dfs[[i]][, c("sat","fr_cat","score_msl")])
  y[[i]] <- Surv(dfs[[i]]$FU_time, dfs[[i]]$critical)
}

alpha0<-NULL;mod_cva<-NULL;alpha<-NULL
set.seed(123)
for (i in 1:5) {
  ALPHA <- seq(0,1, by=0.01)
  mod_cva[[i]] <- cva.glmnet(x[[i]], y[[i]],family="cox",alpha = ALPHA, k=10)
  alpha[[i]] <- ALPHA[which.min(sapply(mod_cva[[i]]$modlist, function(mod) min(mod$cvm)))];alpha
}
plot(mod_cva[[5]])
alpha0<-mean(alpha)

set.seed(123);m1 <- mami(nutri3_imp,
           model="cox", outcome=c("FU_time","critical"),
           method="LASSO", alpha=alpha0, kfold=10,report.exp=FALSE,
           add.factor = c("sat", "fr_cat"), B=100,
           var.remove=c("fr","msl", "sao2aa", "pin", "edad",
                        "sexo","fechavaloracioninicio"))
summary(m1)
plot.mami(m1)
coef<-round(m1$coefficients.s[,1]/min(abs(m1$coefficients.s[,1])))
coef
#### Score estimation in both cohorts ####
nutri3_1<-complete(nutri3_imp, 1)
nutri3_1$score_crit<-(nutri3_1$fr_cat==2)*2+(nutri3_1$fr_cat==1)+nutri3_1$score_msl+(nutri3_1$sat==1)*5+(nutri3_1$sat==2)*9

m0<-coxph(Surv(FU_time, critical)~score_crit, data=nutri3_1)
summary(m0)
cbind(m0$concordance[6],m0$concordance[6]-(m0$concordance[7]*1.96),
      m0$concordance[6]+(m0$concordance[7]*1.96))

m0<-cph(Surv(FU_time, critical)~score_crit, data=nutri3_1, x=T, y=T, surv=T,time.inc=15)
set.seed(123);rms::validate(fit = m0, method="boot", B=1000)

c1<-calibrate(m0, cmethod=c('hare'),
          method="boot", u=15,B=1000, 
          bw=FALSE, rule="aic", type="residual", estimates=TRUE,
          pr=FALSE, what="observed-predicted")
g1<-as.ggplot(~plot(c1))

### Validation ###
nutri4_1  <- nutri0[nutri0$fechavaloracioninicio>"2020-06-04",]
nutri4_1$score_msl[nutri4_1$msl<=0]<-0
nutri4_1$score_msl[nutri4_1$msl>=1 & nutri4_1$msl<=3]<-1;nutri4_1$score_msl[nutri4_1$msl>=4 & nutri4_1$msl<=7]<-2
nutri4_1$score_msl[nutri4_1$msl>=8 & nutri4_1$msl<=11]<-3;nutri4_1$score_msl[nutri4_1$msl>=12]<-4
nutri4_1$fr_cat[nutri4_1$fr>=20 & nutri4_1$fr<30]<-1;nutri4_1$fr_cat[nutri4_1$fr<20]<-0; nutri4_1$fr_cat[nutri4_1$fr>=30]<-2
nutri4_1$sat[nutri4_1$sao2aa<85]<-2;nutri4_1$sat[nutri4_1$sao2aa>=85 & nutri4_1$sao2aa<92]<-1;nutri4_1$sat[nutri4_1$sao2aa>=92]<-0
nutri4_1$score_crit<-(nutri4_1$fr_cat==2)*2+nutri4_1$score_msl+(nutri4_1$sat==1)*5+(nutri4_1$sat==2)*9

m0<-coxph(Surv(FU_time, critical)~score_crit, data=nutri4_1)
cbind(m0$concordance[6],m0$concordance[6]-(m0$concordance[7]*1.96),
      m0$concordance[6]+(m0$concordance[7]*1.96))

m0<-cph(Surv(FU_time, critical)~score_crit, data=nutri4_1, x=T, y=T)
set.seed(123);rms::validate(fit = m0, method="boot", B=1000)
c1<-calibrate(m0, cmethod=c('hare'),
              method="boot", u=15,B=100, 
              bw=FALSE, rule="aic", type="residual", estimates=TRUE,
              pr=FALSE, what="observed-predicted")
g2<-as.ggplot(~plot(c1))
g2
fig2<-ggarrange(g1, g2, labels = LETTERS[1:2])

#### Time-dependent cut-offs ####
nutri3_1<-as.data.frame(nutri3_1)
m1 <- optimal.cutpoints(X = "score_crit", status = "critical", tag.healthy = 0, 
                        methods = "Youden", data = nutri3_1, pop.prev = NULL,
                        control = control.cutpoints(), ci.fit = TRUE, 
                        conf.level = 0.95, trace = FALSE)
summary(m1)

table(nutri3_1$critical, (nutri3_1$score_crit>4))
table(nutri4_1$critical, (nutri4_1$score_crit>4))

m1<- SeSpPPVNPV(cutpoint=8, T=nutri3_1$FU_time,
                delta=nutri3_1$critical,marker=nutri3_1$score_crit,
                cause=1,weighting="marginal",
                times=c(7,10,15,20,30), iid=T)
print(m1)

m1<- SeSpPPVNPV(cutpoint=8, T=nutri4_1$FU_time,
                delta=nutri4_1$critical,marker=nutri4_1$score_crit,
                cause=1,weighting="marginal",
                times=c(7,10,15,20,30), iid=T)
print(m1)
#### Estimate overall ####
hist(nutri3_1$score_crit)
nutri3_1$crit_cat[nutri3_1$score_crit<5]<-1;nutri3_1$crit_cat[nutri3_1$score_crit>=5 & nutri3_1$score_crit<9]<-2
nutri3_1$crit_cat[nutri3_1$score_crit>=9 & nutri3_1$score_crit<12]<-3;nutri3_1$crit_cat[nutri3_1$score_crit>=12]<-4
nutri4_1$crit_cat[nutri4_1$score_crit<5]<-1;nutri4_1$crit_cat[nutri4_1$score_crit>=5 & nutri4_1$score_crit<9]<-2
nutri4_1$crit_cat[nutri4_1$score_crit>=9 & nutri4_1$score_crit<12]<-3;nutri4_1$crit_cat[nutri4_1$score_crit>=12]<-4
table(nutri4_1$crit_cat, nutri4_1$critical)

## Kaplan-Meier trainint ###
mod1_km<-survfit(Surv(FU_time, critical) ~ crit_cat, data = nutri3_1)
mod1_km
summary(mod1_km)
fig1<-ggsurvplot(mod1_km, data = nutri3_1, size = 1,palette = "bw",conf.int = T,
                 risk.table = T,pval = TRUE,ggtheme = theme_classic(),xlab="Time (Days)",
                 ylab="Survival probability",
                 legend.labs = c("Mild-Risk", 
                                 "Moderate-Risk",
                                 "High-Risk",
                                 "Very High Risk"),
                 ylim= c(0.3,1.0),
                 xlim=c(0, 21),
                 break.y.by= c(0.1),
                 break.x.by= c(3),
                 pval.coord = c(0, 0.4))+theme_survminer(base_size = 9,
                                                         base_family = "Arial")

fig1a<-ggarrange(fig1$plot, fig1$table, heights = c(2, 0.7),
                 ncol = 1, nrow = 2)

mod1_km<-survfit(Surv(FU_time, critical) ~ crit_cat, data = nutri4_1)
mod1_km
summary(mod1_km)
fig1<-ggsurvplot(mod1_km, data = nutri4_1, size = 1,palette = "bw",conf.int = T,
                 risk.table = T,pval = TRUE,ggtheme = theme_classic(),xlab="Time (Days)",
                 ylab="Survival probability",
                 legend.labs = c("Mild-Risk", 
                                 "Moderate-Risk",
                                 "High-Risk",
                                 "Very High Risk"),
                 ylim= c(0.3,1.0),
                 xlim=c(0, 21),
                 break.y.by= c(0.1),
                 break.x.by= c(3),
                 pval.coord = c(0, 0.4))+theme_survminer(base_size = 9,
                                                          base_family = "Arial")
hist(m0$linear.predictors)

fig1b<-ggarrange(fig1$plot, fig1$table, heights = c(2, 0.7),
                 ncol = 1, nrow = 2)

### Risk prediction ####
nutri0$score_msl[nutri0$msl<=0]<-0
nutri0$score_msl[nutri0$msl>=1 & nutri0$msl<=3]<-1;nutri0$score_msl[nutri0$msl>=4 & nutri0$msl<=7]<-2
nutri0$score_msl[nutri0$msl>=8 & nutri0$msl<=11]<-3;nutri0$score_msl[nutri0$msl>=12]<-4
nutri0$fr_cat[nutri0$fr>24 & nutri0$fr<30]<-1;nutri0$fr_cat[nutri0$fr<=24]<-0; nutri0$fr_cat[nutri0$fr>=30]<-2
nutri0$sat[nutri0$sao2aa<85]<-2;nutri0$sat[nutri0$sao2aa>=85 & nutri0$sao2aa<92]<-1;nutri0$sat[nutri0$sao2aa>=92]<-0
nutri0$score_crit<-(nutri0$fr_cat==2)*2+nutri0$score_msl+(nutri0$sat==1)*5+(nutri0$sat==2)*9

nutri0$critical1<-factor(nutri$critical, labels = c("Non-severe", "Severe"))
fig1c<-ggplot(nutri0, aes(x=score_crit, fill=critical1, color=critical1)) +
  geom_histogram(aes(y=..density..), position="identity", alpha=0.5)+
  geom_vline(xintercept = 4, linetype="dotted")+
  geom_density(alpha=0.6)+labs(fill="Severe COVID-19", col="Severe COVID-19")+
  theme_classic()+ylab("Density")+xlab("Nutri-CoV score")

#### Comparison with qSOFA, ROX, NEWS y ABC-GOALS ####
nutri0$spo2_fio2<-nutri0$sao2aa/nutri0$fio2
nutri0$rr_news2[nutri0$fr<=8 | nutri0$fr>=25]<-3;nutri0$rr_news2[nutri0$fr>=9 & nutri0$fr<=11]<-1
nutri0$rr_news2[nutri0$fr>=12 & nutri0$fr<=20]<-0;nutri0$rr_news2[nutri0$fr>=21 & nutri0$fr<=24]<-2
nutri0$sp_news2[nutri0$sao2aa>=96]<-0;nutri0$sp_news2[nutri0$sao2aa>=94 & nutri0$sao2aa<=95]<-1
nutri0$sp_news2[nutri0$sao2aa>=92 & nutri0$sao2aa<=93]<-2;nutri0$sp_news2[nutri0$sao2aa<=91]<-3
nutri0$ta_news2[nutri0$tas<=90 | nutri0$tas>=220]<-3;nutri0$ta_news2[nutri0$tas>=91 & nutri0$tas<=100]<-2
nutri0$ta_news2[nutri0$tas>=101 & nutri0$tas<=110]<-1; nutri0$ta_news2[nutri0$tas>=111 & nutri0$tas<=219]<-0
nutri0$hr_news2[nutri0$fc<=40 | nutri0$fc>=131]<-3; nutri0$hr_news2[(nutri0$fc<=50 & nutri0$fc>=41) | (nutri0$fc<=110 & nutri0$fc>=91)]<-1
nutri0$hr_news2[nutri0$fc>=111 & nutri0$fc<=130]<-1; nutri0$hr_news2[nutri0$fc>=51 & nutri0$fc<=90]<-0
nutri0$t_news2[nutri0$temp<=35]<-3;nutri0$t_news2[nutri0$temp>=39.1]<-2
nutri0$t_news2[(nutri0$temp<=36.0 & nutri0$temp>=35.1) | (nutri0$temp>=38.1 & nutri0$temp<=39.0)]<-1
nutri0$t_news2[nutri0$temp>=36.1 & nutri0$temp<=38.0]<-0
nutri0$char_cat[nutri0$Charlson==0]<-0;nutri0$char_cat[nutri0$Charlson>0 & nutri0$Charlson<3]<-1
nutri0$fr_cat2[nutri0$fr<24]<-0;nutri0$fr_cat2[nutri0$fr>=24 & nutri0$fr<29]<-1
nutri0$fr_cat2[nutri0$fr>28]<-3
nutri0$char_cat[nutri0$Charlson>=3]<-3
nutri0$tas_cat[nutri0$tas>=100]<-0;nutri0$tas_cat[nutri0$tas<100]<-1
nutri0$news2<-nutri0$t_news2+nutri0$hr_news2+nutri0$rr_news2+nutri0$ta_news2+nutri0$sp_news2+(nutri$ecgme0sde15<3)*3+nutri0$o2 *2
nutri0$abc_goals<-as.numeric(nutri0$sexo)+nutri0$tas_cat*4+nutri0$disnea+nutri0$char_cat+nutri0$fr_cat2+(nutri0$imc>30)*2

nutri5  <- nutri0[nutri$fechavaloracioninicio>"2020-06-04",]
roc1<-timeROC(T=nutri5$FU_time,
              delta=nutri5$critical,marker=nutri5$score_crit,
              cause=1,weighting="cox",
              other_markers=as.matrix(nutri5$edad, nutri5$sexo),
              times=c(7,15,30));roc1
m1<-coxph(Surv(FU_time, critical)~score_crit, data=nutri5)
cbind(m1$concordance[6],m1$concordance[6]-(m1$concordance[7]*1.96),
      m1$concordance[6]+(m1$concordance[7]*1.96))

roc2<-timeROC(T=nutri5$FU_time,
              delta=nutri5$critical,marker=nutri5$msl,
              cause=1,weighting="cox",
              other_markers=as.matrix(nutri5$edad, nutri5$sexo),
              times=c(7,15,30));roc2
m1<-coxph(Surv(FU_time, critical)~msl, data=nutri5)
cbind(m1$concordance[6],m1$concordance[6]-(m1$concordance[7]*1.96),
      m1$concordance[6]+(m1$concordance[7]*1.96))

roc3<-timeROC(T=nutri5$FU_time,
              delta=nutri5$critical,marker=nutri5$abc_goals,
              other_markers=as.matrix(nutri5$edad, nutri5$sexo),
              cause=1,weighting="cox",
              times=c(7,15,30));roc3
m1<-coxph(Surv(FU_time, critical)~abc_goals, data=nutri5)
cbind(m1$concordance[6],m1$concordance[6]-(m1$concordance[7]*1.96),
      m1$concordance[6]+(m1$concordance[7]*1.96))

roc4<-timeROC(T=nutri5$FU_time,
              delta=nutri5$critical,
              marker=nutri5$Charlson,
              cause=1,weighting="cox",
              other_markers=as.matrix(nutri5$edad, nutri5$sexo),
              times=c(7,15,30));roc4
m1<-coxph(Surv(FU_time, critical)~Charlson, data=nutri5)
cbind(m1$concordance[6],m1$concordance[6]-(m1$concordance[7]*1.96),
      m1$concordance[6]+(m1$concordance[7]*1.96))

roc5<-timeROC(T=nutri5$FU_time,
              delta=nutri5$critical,
              marker=nutri5$qsofa,
              cause=1,weighting="cox",
              other_markers=as.matrix(nutri5$edad, nutri5$sexo),
              times=c(7,15,30));roc5
m1<-coxph(Surv(FU_time, critical)~qsofa, data=nutri5)
cbind(m1$concordance[6],m1$concordance[6]-(m1$concordance[7]*1.96),
      m1$concordance[6]+(m1$concordance[7]*1.96))

roc6<-timeROC(T=nutri5$FU_time,
              delta=nutri5$critical,
              marker=nutri5$news,
              cause=1,weighting="cox",
              other_markers=as.matrix(nutri5$sexo),
              times=c(7,15,30));roc6
m1<-coxph(Surv(FU_time, critical)~news, data=nutri5)
cbind(m1$concordance[6],m1$concordance[6]-(m1$concordance[7]*1.96),
      m1$concordance[6]+(m1$concordance[7]*1.96))

roc7<-timeROC(T=nutri5$FU_time,
              delta=nutri5$critical,
              marker=nutri5$news2,
              cause=1,weighting="cox",
              other_markers=as.matrix(nutri5$sexo),
              times=c(7,15,30));roc7
m1<-coxph(Surv(FU_time, critical)~news2, data=nutri5)
cbind(m1$concordance[6],m1$concordance[6]-(m1$concordance[7]*1.96),
      m1$concordance[6]+(m1$concordance[7]*1.96))

#### Decision curves ####
library(rmda)
m1 <- decision_curve(critical~score_crit,data = nutri5,
                     thresholds = seq(0, 1, by = .01),
                     confidence.intervals = 'none')
m2 <- decision_curve(critical~msl,data = nutri5,
                     thresholds = seq(0, 1, by = .01),
                     confidence.intervals = 'none')
m3 <- decision_curve(critical~qsofa,data = nutri5,
                     thresholds = seq(0, 1, by = .01),
                     confidence.intervals = 'none')
m4 <- decision_curve(critical~abc_goals,data = nutri5,
                     thresholds = seq(0, 1, by = .01),
                     confidence.intervals = 'none')
m5 <- decision_curve(critical~Charlson,data = nutri5,
                     thresholds = seq(0, 1, by = .01),
                     confidence.intervals = 'none')
m6 <- decision_curve(critical~news2,data = nutri5,
                     thresholds = seq(0, 1, by = .01),
                     confidence.intervals = 'none')
m7 <- decision_curve(critical~news,data = nutri5,
                     thresholds = seq(0, 1, by = .01),
                     confidence.intervals = 'none')


colfunc <- colorRampPalette(c("yellow", "black"))
fig1d<-as.ggplot(~plot_decision_curve( list(m1,m2,m3,m4,m5,m7,m6),
                     curve.names = c('Nutri-CoV', 'MSL-COVID-19', 'q-SOFA', 
                                     'ABC-GOALS', 'Chalson CI', 'NEWS', 'NEWS-2'),
                     col = c(colfunc(8)),
                     lty = c(2,1,2,1,2,1,2),
                     lwd = c(3,2, 3,2,3,2,3),
                     legend.position = 'topright',
                     standardize = T, xlim = c(0,0.9), ylim=c(-0.2,1)))+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

f1<-ggarrange(fig1a,fig1b, labels = LETTERS[1:2])

ggsave(file = "Figure2.jpg", 
       f1,
       bg = "transparent",
       width = 50, 
       height = 20,
       units=c("cm"),
       dpi = 500,
       limitsize = FALSE)

f3<-ggarrange(fig1c,fig1d, labels = LETTERS[1:2])

ggsave(file = "Figure3.jpg", 
       f3,
       bg = "transparent",
       width = 50, 
       height = 20,
       units=c("cm"),
       dpi = 500,
       limitsize = FALSE)


#### Sensitivity analyses ####
### Overall###
m0<-coxph(Surv(FU_time, critical)~score_crit, data=nutri4_1)
cbind(m0$concordance[6],m0$concordance[6]-(m0$concordance[7]*1.96),
      m0$concordance[6]+(m0$concordance[7]*1.96))

m0<-cph(Surv(FU_time, critical)~score_crit, data=nutri4_1, x=T, y=T, surv=T,time.inc=15)
set.seed(123);rms::validate(fit = m0, method="boot", B=1000)

### Non-imputed ###
nutri6<-nutri2[nutri$fechavaloracioninicio>"2020-06-04",]
nutri6$score_msl[nutri6$msl<=0]<-0
nutri6$score_msl[nutri6$msl>=1 & nutri6$msl<=3]<-1;nutri6$score_msl[nutri6$msl>=4 & nutri6$msl<=7]<-2
nutri6$score_msl[nutri6$msl>=8 & nutri6$msl<=11]<-3;nutri6$score_msl[nutri6$msl>=12]<-4
nutri6$fr_cat[nutri6$fr>24 & nutri6$fr<30]<-1;nutri6$fr_cat[nutri6$fr<=24]<-0; nutri6$fr_cat[nutri6$fr>=30]<-2
nutri6$sat[nutri6$sao2aa<85]<-2;nutri6$sat[nutri6$sao2aa>=85 & nutri6$sao2aa<92]<-1;nutri6$sat[nutri6$sao2aa>=92]<-0
nutri6$score_crit<-(nutri6$fr_cat==2)*2+nutri6$score_msl+(nutri6$sat==1)*5+(nutri6$sat==2)*9

### Non-imputed ###
m0<-coxph(Surv(FU_time, critical)~score_crit, data=nutri6)
summary(m0)
cbind(m0$concordance[6],m0$concordance[6]-(m0$concordance[7]*1.96),
      m0$concordance[6]+(m0$concordance[7]*1.96))

m0<-cph(Surv(FU_time, critical)~score_crit, data=nutri6, x=T, y=T, surv=T,time.inc=15)
set.seed(123);rms::validate(fit = m0, method="boot", B=1000)

### Older than 60 years ###
m0<-coxph(Surv(FU_time, critical)~score_crit, data=nutri5 %>% filter(edad>60) )
cbind(m0$concordance[6],m0$concordance[6]-(m0$concordance[7]*1.96),
      m0$concordance[6]+(m0$concordance[7]*1.96))

m0<-cph(Surv(FU_time, critical)~score_crit, data=nutri5 %>% filter(edad>60), x=T, y=T, surv=T,time.inc=15)
set.seed(123);rms::validate(fit = m0, method="boot", B=1000)

### Younger than 60 years ###
m0<-coxph(Surv(FU_time, critical)~score_crit, data=nutri5 %>% filter(edad<=60) )
cbind(m0$concordance[6],m0$concordance[6]-(m0$concordance[7]*1.96),
      m0$concordance[6]+(m0$concordance[7]*1.96))

m0<-cph(Surv(FU_time, critical)~score_crit, data=nutri5 %>% filter(edad<=60), x=T, y=T, surv=T,time.inc=15)
set.seed(123);rms::validate(fit = m0, method="boot", B=1000)

### No comorbidities ###
m0<-coxph(Surv(FU_time, critical)~score_crit, data=nutri5 %>% filter(comorb==0))
cbind(m0$concordance[6],m0$concordance[6]-(m0$concordance[7]*1.96),
      m0$concordance[6]+(m0$concordance[7]*1.96))

m0<-cph(Surv(FU_time, critical)~score_crit, data=nutri5 %>% filter(comorb==0), x=T, y=T, surv=T,time.inc=15)
set.seed(123);rms::validate(fit = m0, method="boot", B=1000)

### Comorbidities ###
m0<-coxph(Surv(FU_time, critical)~score_crit, data=nutri5 %>% filter(comorb==1))
cbind(m0$concordance[6],m0$concordance[6]-(m0$concordance[7]*1.96),
      m0$concordance[6]+(m0$concordance[7]*1.96))

m0<-cph(Surv(FU_time, critical)~score_crit, data=nutri5 %>% filter(comorb==1), x=T, y=T, surv=T,time.inc=15)
set.seed(123);rms::validate(fit = m0, method="boot", B=1000)

### Outcome at admission ###
m0<-coxph(Surv(FU_time, critical)~score_crit, data=nutri5 %>% filter(hospitalizacion==1))
cbind(m0$concordance[6],m0$concordance[6]-(m0$concordance[7]*1.96),
      m0$concordance[6]+(m0$concordance[7]*1.96))

m0<-cph(Surv(FU_time, critical)~score_crit, data=nutri5 %>% filter(hospitalizacion==1), x=T, y=T, surv=T,time.inc=15)
set.seed(123);rms::validate(fit = m0, method="boot", B=1000)

### Cases who were not discharged ###
nutri$score_msl[nutri$msl<=0]<-0
nutri$score_msl[nutri$msl>=1 & nutri$msl<=3]<-1;nutri$score_msl[nutri$msl>=4 & nutri$msl<=7]<-2
nutri$score_msl[nutri$msl>=8 & nutri$msl<=11]<-3;nutri$score_msl[nutri$msl>=12]<-4
nutri$fr_cat[nutri$fr>=20 & nutri$fr<30]<-1;nutri$fr_cat[nutri$fr<20]<-0; nutri$fr_cat[nutri$fr>=30]<-2
nutri$sat[nutri$sao2aa<85]<-2;nutri$sat[nutri$sao2aa>=85 & nutri$sao2aa<92]<-1;nutri$sat[nutri$sao2aa>=92]<-0
nutri$score_crit<-(nutri$fr_cat==2)*2+nutri$score_msl+(nutri$sat==1)*5+(nutri$sat==2)*9
sum(!is.na(nutri$fechaegreso))
m0<-coxph(Surv(FU_time, critical)~score_crit, data=nutri %>% filter(fechavaloracioninicio>"2020-06-04" & (!is.na(fechaegreso) | hospitalizacion==0)))
cbind(m0$concordance[6],m0$concordance[6]-(m0$concordance[7]*1.96),
      m0$concordance[6]+(m0$concordance[7]*1.96))
summary(m0)
m0<-cph(Surv(FU_time, critical)~score_crit, data=nutri %>% filter(fechavaloracioninicio>"2020-06-04" & (!is.na(fechaegreso) | hospitalizacion==0)), x=T, y=T, surv=T,time.inc=15)
set.seed(123);rms::validate(fit = m0, method="boot", B=1000)
