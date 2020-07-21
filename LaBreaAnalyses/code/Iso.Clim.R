options(scipen=999)
library(ggplot2)
library(tidyverse)
library(scales)
library(RColorBrewer)

iso_dat<-read.csv("data/processed/Interpolation/IsoClim_interpolated_rev.csv")# All rabbits and Squirrels with mean calibrated ages < 50 ka BP

####Intraspecific diet ###

## Compare desert cottontails to brush rabbits

iso_rabbits<- iso_dat[grep("S. audubonii|S. bachmani", iso_dat$Species), c(3:6)]

#data normally distributed?
C_aud <- iso_rabbits$d13C[which(iso_rabbits$Species=="S. audubonii")] #yes
C_bach <- iso_rabbits$d13C[which(iso_rabbits$Species=="S. bachmani")] #No
N_aud <- iso_rabbits$d15N[which(iso_rabbits$Species=="S. audubonii")] #sort of
N_bach <- iso_rabbits$d15N[which(iso_rabbits$Species=="S. bachmani")] #sort of
#Use non-parametric tests

library(e1071)
library(ggpubr)

ggqqplot(N_bach)
shapiro.test(N_bach) # if > 0.05, then yes  
kurtosis(N_bach) # if within -2 and 2, then yes 
skewness(C_bach) # if within -2 and 2, then yes
mean(C_bach)
median(C_bach) #if similar to mean, than yes 


#carbon
pdf("output/Isoplots/Rabbit_spp_d13C.pdf", width=5, height=4)
boxplot(d13C~Species, data=iso_rabbits, ylab= expression({delta}^13*C~'\u2030'))
dev.off()
wilcox.test(d13C~Species, data=iso_rabbits, correct = TRUE, alternative = "two.sided")
#t.test(d13C~Taxon, mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F, data=iso_rabbits)

#nitrogen
pdf("output/Isoplots/Rabbit_spp_d15N.pdf", width=5, height=4)
boxplot(d15N~Species,data=iso_rabbits, ylab=expression({delta}^15*N~'\u2030'))
dev.off()
wilcox.test(d15N~Species, data=iso_rabbits, correct = TRUE)


## compare all rabbits to squirrels
C_syl <- iso_dat$d13C[which(iso_dat$Taxon=="Sylvilagus ")] #No
C_oto <- iso_dat$d13C[which(iso_dat$Taxon=="Otospermophilus")] #No
N_syl <- iso_dat$d15N[which(iso_dat$Taxon=="Sylvilagus ")] #Yes
N_oto <- iso_dat$d15N[which(iso_dat$Taxon=="Otospermophilus")] #sort of
#Use non-parametric test

ggqqplot(N_oto)
shapiro.test(N_oto) # if > 0.05, then yes  
kurtosis(N_oto) # if within -2 and 2, then yes 
skewness(N_oto) # if within -2 and 2, then yes
mean(N_oto)
median(N_oto) #if similar to mean, than yes 


#carbon
pdf("output/Isoplots/Rabbit_squirrel_d13C.pdf", width=5, height=4)
boxplot(d13C~Taxon,data=iso_dat, ylab=expression({delta}^13*C~'\u2030'))
dev.off()
wilcox.test(d13C~Taxon, data=iso_dat, correct = TRUE)

#nitrogen
pdf("output/Isoplots/Rabbit_squirrel_d15N.pdf", width=5, height=4)
boxplot(d15N~Taxon,data=iso_dat, ylab=expression({delta}^15*N~'\u2030'))
dev.off()
wilcox.test(d15N~Taxon, data=iso_dat, correct = TRUE)


#### Isotopes and Climate ###

# Methods potential text:
# For carbon and nitrogen, I fit a linear model that originally included oxygen, taxon, and the interaction between the two variables as independent variables. I then performed stepwise regression to determine a final model. 
# Results potential text
# Stepwise regression indicated that there was no significant interaction between 18O and taxon for either carbon or nitrogen stable isotope values. For carbon, variation in 13C was significantly associated with both 18O and taxon (stats from summary(carbon.lm.final)). For nitrogen, neither 18O nor taxon explained significant variation in 15N, though taxon as a single variable was marginally significant (stats from summary(nitrogen.lm.taxon)).


# models - Carbon

# model with interaction term included
carbon.lm.all.interaction<-lm(d13C~pach.d18O_mean * Taxon, data=iso_dat)
summary(carbon.lm.all.interaction)

# model with no interaction term 
# NOTE: this is what the final model is in the end, so THIS  IS WHAT YOU SHOULD REPORT IN THE PAPER.
carbon.lm.all<-lm(d13C~pach.d18O_mean + Taxon, data=iso_dat)
summary(carbon.lm.all)

# climate-only model
carbon.lm.clim<-lm(d13C~pach.d18O_mean, data=iso_dat)
summary(carbon.lm.clim)

# taxon-only model
carbon.lm.taxon<-lm(d13C~Taxon, data=iso_dat)
summary(carbon.lm.taxon)

#stepwise regression 
carbon.lm.final <- step(lm(d13C~pach.d18O_mean*Taxon, data=iso_dat, direction="both"))
summary(carbon.lm.final)
# --> this shows that the d180 + Taxon model (no interaction) is the best final model. 

# t-test of residuals from climate only model
c.t.clim <- t.test(carbon.lm.clim$residuals ~ iso_dat$Taxon[which(!is.na(iso_dat$pach.d18O_mean))])


# models - Nitrogen

# model with interaction term included
nitrogen.lm.all.interaction<-lm(d15N~pach.d18O_mean * Taxon, data=iso_dat)
summary(nitrogen.lm.all.interaction)

# model with no interaction term 
nitrogen.lm.all.additive<-lm(d15N~pach.d18O_mean + Taxon, data=iso_dat)
summary(nitrogen.lm.all.additive)

# climate-only model
nitrogen.lm.clim<-lm(d15N~pach.d18O_mean, data=iso_dat)
summary(nitrogen.lm.clim)

# taxon-only model
nitrogen.lm.taxon<-lm(d15N~Taxon, data=iso_dat)
summary(nitrogen.lm.taxon)

#stepwise regression --> this shows that there is not a good final model
nitrogen.lm.final <- step(lm(d15N~pach.d18O_mean*Taxon, data=iso_dat, direction="both"))
# --> there is not a good final model, so I am just going to treat nitrogen the same as carbon for plotting.
nitrogen.lm.final <- nitrogen.lm.all.additive

# t-test of residuals from climate-only model
n.t.clim <- t.test(nitrogen.lm.clim$residuals ~ iso_dat$Taxon[which(!is.na(iso_dat$pach.d18O_mean))])


# JESSICA's NEW final plot ####

pdf("output/Isoplots/lm_carbon_nitrogen_all.pdf", width=8, height=8)

par(mfcol=c(2,2), mar=c(4,4,1,1), cex.axis=0.8)

# carbon
plot(d13C~pach.d18O_mean, data=iso_dat, pch=16, type="n", xlab="", ylab="")
points(d13C~pach.d18O_mean, 
       data=iso_dat[which(iso_dat$Taxon == "Sylvilagus "),], 
       pch=16, col="orange")
points(d13C~pach.d18O_mean, 
       data=iso_dat[which(iso_dat$Taxon == "Otospermophilus"),], 
       pch=16, col="royalblue2")
abline(carbon.lm.final, lty=2)
abline(lm(d13C~pach.d18O_mean, data=iso_dat), lty=1)
mtext(expression({delta}^18*O~'\u2030'), side=1, line=2.25) 
mtext(expression({delta}^13*C~'\u2030'), side=2, line=2) 
legend("bottomleft", legend = c("Otospermophilus", "Sylvilagus"),
       col = c("royalblue2","orange"), pch = 16, 
       bty = "n", cex = 0.8)

boxplot(carbon.lm.clim$residuals ~ iso_dat$Taxon[which(!is.na(iso_dat$pach.d18O_mean))], 
        xlab="", ylab="")
stripchart(carbon.lm.clim$residuals ~ iso_dat$Taxon[which(!is.na(iso_dat$pach.d18O_mean))], vertical=TRUE, add=TRUE, method="stack", col=c("royalblue2","orange"), pch=16)
mtext("Taxon", side=1, line=2.25)
mtext(expression('Residuals ('~{delta}^13*C~'\u2030'~' ~ '~{delta}^18*O~'\u2030'~')'), side=2, line=2.25, cex=0.8)
legend("topright", legend=paste0("t=", round(c.t.clim$statistic,2), "; df=", round(c.t.clim$parameter,2), "; p=", round(c.t.clim$p.value,2)), bty = "n", cex = 0.8)


# nitrogen
plot(d15N~pach.d18O_mean, data=iso_dat, pch=16, 
     xlab = "", ylab = "", type="n")
points(d15N~pach.d18O_mean, 
       data=iso_dat[which(iso_dat$Taxon == "Sylvilagus "),], 
       pch=16, col="orange")
points(d15N~pach.d18O_mean, 
       data=iso_dat[which(iso_dat$Taxon == "Otospermophilus"),], 
       pch=16, col="royalblue2")
abline(nitrogen.lm.final, lty=2)
abline(lm(d15N~pach.d18O_mean, data=iso_dat), lty=1)
legend("topleft", legend = c("Otospermophilus", "Sylvilagus"),
       col = c("royalblue2","orange"), pch = 16, 
       bty = "n", cex = 0.8)
mtext(expression({delta}^18*O~'\u2030'), side=1, line=2.25)
mtext(expression({delta}^15*N~'\u2030'), side=2, line=2.25)

boxplot(nitrogen.lm.clim$residuals ~ iso_dat$Taxon[which(!is.na(iso_dat$pach.d18O_mean))],
        xlab="", ylab="")
stripchart(nitrogen.lm.clim$residuals ~ iso_dat$Taxon[which(!is.na(iso_dat$pach.d18O_mean))], vertical=TRUE, add=TRUE, method="stack", col=c("royalblue2","orange"), pch=16)
mtext("Taxon", side=1, line=2.25)
mtext(expression('Residuals ('~{delta}^15*N~'\u2030'~' ~ '~{delta}^18*O~'\u2030'~')'), side=2, line=2.25, cex=0.8)
legend("topright", legend=paste0("t=", round(n.t.clim$statistic,2), "; df=", round(n.t.clim$parameter,2), "; p=", round(n.t.clim$p.value,2)), bty = "n", cex = 0.8)

dev.off()

# I haven't put A-D labels on them yet.
# Figure caption text.
# Figure 4. The relationship between isotope niche and climate. A) & C) show the relationship between 13C or 15N, respectively, and 18O.  In both panels, the dashed line indicates the fitted relationship between 13C or 15N and 18O from the final model which  includes Taxon as an independent variable. The solid line indicates the fitted relationship between 13C or 15N and 18O from a linear model that just includes 18O as the independent variable. B) & D) indicate the residuals from the climate-only model, plotted by taxon.



### all taxa ###

#carbon and climate
pdf("output/Isoplots/lm_carbon_all.pdf", width=5, height=4)
plot(d13C~pach.d18O_mean, data=iso_dat, pch=16, 
     xlab = expression({delta}^18*O~'\u2030', ), ylab = expression({delta}^13*C~'\u2030'), type="n")
points(d13C~pach.d18O_mean, data=iso_dat[which(iso_dat$Taxon == "Sylvilagus "),], pch=16, col="orange")
points(d13C~pach.d18O_mean, data=iso_dat[which(iso_dat$Taxon == "Otospermophilus"),], pch=16, col="royalblue2")
abline(carbon.lm.final)
abline(lm(d13C~pach.d18O_mean, data=iso_dat), lty=2, col="red")
legend("bottomleft", legend = c("Otospermophilus", "Sylvilagus"),
       col = c("royalblue2","orange"), pch = 16, 
       bty = "n", cex = 0.8)
dev.off()


# plot the residuals from the linear model by taxon
boxplot(carbon.lm$residuals ~ iso_dat$Taxon[which(!is.na(iso_dat$pach.d18O_mean))],
        xlab="Taxon", 
        ylab=expression('Residuals from the '~{delta}^13*C~'\u2030'~' linear model'))
stripchart(carbon.lm$residuals ~ iso_dat$Taxon[which(!is.na(iso_dat$pach.d18O_mean))], vertical=TRUE, add=TRUE, method="stack", col=c("royalblue2","orange"), pch=16)
t.test(carbon.lm$residuals ~ iso_dat$Taxon[which(!is.na(iso_dat$pach.d18O_mean))])

ccf(iso_dat$pach.d18O_mean, iso_dat$d13C)

#nitrogen and climate
pdf("output/Isoplots/lm_nitrogen_all.pdf", width=5, height=4)
plot(d15N~pach.d18O_mean, data=iso_dat, pch=16,
     xlab = expression({delta}^18*O~'\u2030', ), ylab = expression({delta}^15*N~'\u2030'))
nitrogen.lm<-lm(d15N~pach.d18O_mean, data=iso_dat)
summary(nitrogen.lm)
abline(lm(d15N~pach.d18O_mean, data=iso_dat))
dev.off()

ccf(data0$pach.d18O_mean, data0$d15N)

#plot carbon and climate thru time

ylim.prim <- c(0, 3)   
ylim.sec <- c(-23, -17)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

pdf("output/Isoplots/carbon_climate_time.pdf", width=7, height=4)
C<-ggplot(iso_dat, aes(Calibrated_mean_age, pach.d18O_mean)) +
  geom_line(aes(y = a + d13C*b), color = "brown") +
  geom_line(aes(y = pach.d18O_mean), color = "blue") +
  scale_y_reverse(name = expression({delta}^18*O~'\u2030'), sec.axis = sec_axis(~ (. - a)/b, name = expression({delta}^13*C~'\u2030'))) +
  scale_x_reverse("Calibrated Years Before Present")+
  theme_light() 
C + theme(axis.text.y.left = element_text(color="blue"),
          axis.text.y.right = element_text(color="brown"))
 dev.off()
 
#plot nitrogen thru time
ylim.prim <- c(0, 3)   
ylim.sec <- c(3, 14)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

p<-ggplot(data0, aes(Calibrated_mean_age, pach.d18O_mean)) +
  geom_line(aes(y = a + d15N*b), color = "red") +
  geom_line(aes(y = pach.d18O_mean), color = "blue") +
  scale_y_reverse("d18O Neogloboquadrina pachyderma", sec.axis = sec_axis(~ (. - a)/b, name = expression({delta}^15*N~'\u2030'))) +
  scale_x_reverse("Calibrated Years Before Present")+
  ggtitle("RLB Small Mammal Nitrogen and Climate")+
  theme_light() 
p + theme(axis.text.y.left = element_text(color="blue"),
          axis.text.y.right = element_text(color="red"))


### Just squirrles ###

#carbon and climate
pdf("output/Isoplots/lm_carbon_Oto.pdf", width=5, height=4)
plot(d13C[Taxon=='Otospermophilus']~pach.d18O_mean[Taxon=='Otospermophilus'], 
     xlab = expression({delta}^18*O~'\u2030', ), ylab = expression({delta}^13*C~'\u2030'),
     data=iso_dat, pch=16)
carbon.lm<-lm(d13C[Taxon=='Otospermophilus']~pach.d18O_mean[Taxon=='Otospermophilus'], data=iso_dat)
summary(carbon.lm)
abline(lm(d13C[Taxon=='Otospermophilus']~pach.d18O_mean[Taxon=='Otospermophilus'], data=iso_dat))
dev.off()

#ccf(iso_dat$pach.d18O[Taxon=='O. beecheyi'], iso_dat$d13C[Taxon=='O. beecheyi'])

#nitrogen and climate
pdf("output/Isoplots/lm_nitrogen_Oto.pdf", width=5, height=4)
plot(d15N[Taxon=='Otospermophilus']~pach.d18O_mean[Taxon=='Otospermophilus'], 
     xlab = expression({delta}^18*O~'\u2030', ), ylab = expression({delta}^15*N~'\u2030'), 
     data=iso_dat, pch=16)
carbon.lm<-lm(d15N[Taxon=='Otospermophilus']~pach.d18O_mean[Taxon=='Otospermophilus'], data=iso_dat)
summary(carbon.lm)
abline(lm(d15N[Taxon=='Otospermophilus']~pach.d18O_mean[Taxon=='Otospermophilus'], data=iso_dat))
dev.off()

#ccf(iso_dat$pach.d18O[Taxon=='O. beecheyi'], iso_dat$d13C[Taxon=='O. beecheyi'])

#plot squirrel carbon and climate thru time
ylim.prim <- c(0, 3)   
ylim.sec <- c(-23, -17)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

C<-ggplot(data0, aes(Cal_age, pach.d18O)) +
  geom_line(aes(y = a + d13C*b), color = "brown", data=subset(data0,Taxon=="O. beecheyi")) +
  geom_line(aes(y = pach.d18O), color = "blue", data=subset(data0,Taxon=="O. beecheyi")) +
  scale_y_reverse("d18O Neogloboquadrina pachyderma", sec.axis = sec_axis(~ (. - a)/b, name = expression({delta}^13*C~'\u2030'))) +
  scale_x_reverse("Calibrated Years Before Present")+
  ggtitle("RLB Squirrel Carbon and Climate")+
  theme_light() 
C + theme(axis.text.y.left = element_text(color="blue"),
          axis.text.y.right = element_text(color="brown"))

#plot squirrel nitrogen thru time
ylim.prim <- c(0, 3)   
ylim.sec <- c(3, 14)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

n<-ggplot(data0, aes(Cal_age, pach.d18O)) +
  geom_line(aes(y = a + d15N*b), color = "red", data=subset(data0,Taxon=="O. beecheyi")) +
  geom_line(aes(y = pach.d18O), color = "blue", data=subset(data0,Taxon=="O. beecheyi")) +
  scale_y_reverse("d18O Neogloboquadrina pachyderma", sec.axis = sec_axis(~ (. - a)/b, name = expression({delta}^15*N~'\u2030'))) +
  scale_x_reverse("Calibrated Years Before Present")+
  ggtitle("RLB Squirrel Nitrogen and Climate")+
  theme_light() 
n + theme(axis.text.y.left = element_text(color="blue"),
          axis.text.y.right = element_text(color="red"))

### Rabbits (all) ###

#carbon and climate
pdf("output/Isoplots/lm_carbon_Syl.pdf", width=5, height=4)
plot(d13C[Taxon=='Sylvilagus ']~pach.d18O_mean[Taxon=='Sylvilagus '], 
     xlab = expression({delta}^18*O~'\u2030', ), ylab = expression({delta}^13*C~'\u2030'),
     data=iso_dat, pch=16)
carbon.lm<-lm(d13C[Taxon=='Sylvilagus ']~pach.d18O_mean[Taxon=='Sylvilagus '], data=iso_dat)
summary(carbon.lm)
abline(lm(d13C[Taxon=='Sylvilagus ']~pach.d18O_mean[Taxon=='Sylvilagus '], data=iso_dat))
dev.off()

#nitrogen and climate
pdf("output/Isoplots/lm_nitrogen_Syl.pdf", width=5, height=4)
plot(d15N[Taxon=='Sylvilagus ']~pach.d18O_mean[Taxon=='Sylvilagus '], 
     xlab = expression({delta}^18*O~'\u2030', ), ylab = expression({delta}^15*N~'\u2030'),
     data=iso_dat, pch=16)
nitrogen.lm<-lm(d15N[Taxon=='Sylvilagus ']~pach.d18O_mean[Taxon=='Sylvilagus '], data=iso_dat)
summary(nitrogen.lm)
abline(lm(d15N[Taxon=='Sylvilagus ']~pach.d18O_mean[Taxon=='Sylvilagus '], data=iso_dat))
dev.off()

#plot rabbit carbon and climate thru time
ylim.prim <- c(0, 3)   
ylim.sec <- c(-23, -17)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

C<-ggplot(data0, aes(Cal_age, pach.d18O)) +
  geom_line(aes(y = a + d13C*b), color = "brown", data=subset(data0,Taxon=="Sylvilagus sp")) +
  geom_line(aes(y = pach.d18O), color = "blue", data=subset(data0,Taxon=="Sylvilagus sp")) +
  scale_y_reverse("d18O Neogloboquadrina pachyderma", sec.axis = sec_axis(~ (. - a)/b, name = expression({delta}^13*C~'\u2030'))) +
  scale_x_reverse("Calibrated Years Before Present")+
  ggtitle("RLB Rabbit Carbon and Climate")+
  theme_light() 
C + theme(axis.text.y.left = element_text(color="blue"),
          axis.text.y.right = element_text(color="brown"))

#plot rabbit nitrogen thru time
ylim.prim <- c(0, 3)   
ylim.sec <- c(3, 14)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

n<-ggplot(data0, aes(Cal_age, pach.d18O)) +
  geom_line(aes(y = a + d15N*b), color = "red", data=subset(data0,Taxon=="Sylvilagus sp")) +
  geom_line(aes(y = pach.d18O), color = "blue", data=subset(data0,Taxon=="Sylvilagus sp")) +
  scale_y_reverse("d18O Neogloboquadrina pachyderma", sec.axis = sec_axis(~ (. - a)/b, name = expression({delta}^15*N~'\u2030'))) +
  scale_x_reverse("Calibrated Years Before Present")+
  ggtitle("RLB Rabbit Nitrogen and Climate")+
  theme_light() 
n + theme(axis.text.y.left = element_text(color="blue"),
          axis.text.y.right = element_text(color="red"))


### Rabbits (S. audubonii only) ###

#carbon and climate
plot(d13C[Taxon=='S. audubonii']~pach.d18O[Taxon=='S. audubonii'], data=data0, pch=16)
carbon.lm<-lm(d13C[Taxon=='S. audubonii']~pach.d18O[Taxon=='S. audubonii'], data=data0)
summary(carbon.lm)
abline(lm(d13C[Taxon=='S. audubonii']~pach.d18O[Taxon=='S. audubonii'], data=data0))

#nitrogen and climate
plot(d15N[Taxon=='S. audubonii']~pach.d18O[Taxon=='S. audubonii'], data=data0, pch=16)
carbon.lm<-lm(d15N[Taxon=='S. audubonii']~pach.d18O[Taxon=='S. audubonii'], data=data0)
summary(carbon.lm)
abline(lm(d15N[Taxon=='S. audubonii']~pach.d18O[Taxon=='S. audubonii'], data=data0))

#plot S. audubonii carbon and climate thru time
ylim.prim <- c(0, 3)   
ylim.sec <- c(-23, -17)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

C<-ggplot(data0, aes(Cal_age, pach.d18O)) +
  geom_line(aes(y = a + d13C*b), color = "brown", data=subset(data0,Taxon=="S. audubonii")) +
  geom_line(aes(y = pach.d18O), color = "blue", data=subset(data0,Taxon=="S. audubonii")) +
  scale_y_reverse("d18O Neogloboquadrina pachyderma", sec.axis = sec_axis(~ (. - a)/b, name = expression({delta}^13*C~'\u2030'))) +
  scale_x_reverse("Calibrated Years Before Present")+
  ggtitle("RLB Desert Cottontail Carbon and Climate")+
  theme_light() 
C + theme(axis.text.y.left = element_text(color="blue"),
          axis.text.y.right = element_text(color="brown"))

#plot S. audubonii nitrogen thru time
ylim.prim <- c(0, 3)   
ylim.sec <- c(3, 14)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

n<-ggplot(data0, aes(Cal_age, pach.d18O)) +
  geom_line(aes(y = a + d15N*b), color = "red", data=subset(data0,Taxon=="S. audubonii")) +
  geom_line(aes(y = pach.d18O), color = "blue", data=subset(data0,Taxon=="S. audubonii")) +
  scale_y_reverse("d18O Neogloboquadrina pachyderma", sec.axis = sec_axis(~ (. - a)/b, name = expression({delta}^15*N~'\u2030'))) +
  scale_x_reverse("Calibrated Years Before Present")+
  ggtitle("RLB Desert Cottontail Nitrogen and Climate")+
  theme_light() 
n + theme(axis.text.y.left = element_text(color="blue"),
          axis.text.y.right = element_text(color="red"))

### Rabbits (bachmani only) ###
#carbon and climate
plot(d13C[Taxon=='S. bachmani']~pach.d18O[Taxon=='S. bachmani'], data=data0, pch=16)
carbon.lm<-lm(d13C[Taxon=='S. bachmani']~pach.d18O[Taxon=='S. bachmani'], data=data0)
summary(carbon.lm)
abline(lm(d13C[Taxon=='S. bachmani']~pach.d18O[Taxon=='S. bachmani'], data=data0))

#nitrogen and climate
plot(d15N[Taxon=='S. bachmani']~pach.d18O[Taxon=='S. bachmani'], data=data0, pch=16)
carbon.lm<-lm(d15N[Taxon=='S. bachmani']~pach.d18O[Taxon=='S. bachmani'], data=data0)
summary(carbon.lm)
abline(lm(d15N[Taxon=='S. bachmani']~pach.d18O[Taxon=='S. bachmani'], data=data0))

#plot brush rabbit carbon and climate thru time
ylim.prim <- c(0, 3)   
ylim.sec <- c(-23, -17)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

C<-ggplot(data0, aes(Cal_age, pach.d18O)) +
  geom_line(aes(y = a + d13C*b), color = "brown", data=subset(data0,Taxon=="S. bachmani")) +
  geom_line(aes(y = pach.d18O), color = "blue", data=subset(data0,Taxon=="S. bachmani")) +
  scale_y_reverse("d18O Neogloboquadrina pachyderma", sec.axis = sec_axis(~ (. - a)/b, name = expression({delta}^13*C~'\u2030'))) +
  scale_x_reverse("Calibrated Years Before Present")+
  ggtitle("RLB Brush Rabbit Carbon and Climate")+
  theme_light() 
C + theme(axis.text.y.left = element_text(color="blue"),
          axis.text.y.right = element_text(color="brown"))

#plot S. bachmani nitrogen thru time
ylim.prim <- c(0, 3)   
ylim.sec <- c(3, 14)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

n<-ggplot(data0, aes(Cal_age, pach.d18O)) +
  geom_line(aes(y = a + d15N*b), color = "red", data=subset(data0,Taxon=="S. bachmani")) +
  geom_line(aes(y = pach.d18O), color = "blue", data=subset(data0,Taxon=="S. bachmani")) +
  scale_y_reverse("d18O Neogloboquadrina pachyderma", sec.axis = sec_axis(~ (. - a)/b, name = expression({delta}^15*N~'\u2030'))) +
  scale_x_reverse("Calibrated Years Before Present")+
  ggtitle("Brush Rabbit Nitrogen and Climate")+
  theme_light() 
n + theme(axis.text.y.left = element_text(color="blue"),
          axis.text.y.right = element_text(color="red"))


### Rabbit body size ###

#p3 length and climate
data1<-na.omit(data0)

plot(P3_L~pach.d18O_mean, data=data1, xlim=rev(c(0, 3)), pch=16)
size.lm<-lm(P3_L~pach.d18O_mean, data=data1)
summary(size.lm)
abline(lm(P3_L~pach.d18O_mean, data=data1))
ccf(data1$pach.d18O_mean, data1$P3_L)

pach.d18O_lag = lag(data1$pach.d18O,-6)

size.lm.lag<-lm(P3_L~pach.d18O_lag, data=data1)
summary(size.lm.lag)

#remove bachmani
data1_filtered<-subset(data1,Taxon=="S. audubonii")
plot(P3_L~pach.d18O, data=data1_filtered, pch=16)
size.lm<-lm(P3_L~pach.d18O, data=data1_filtered)
summary(size.lm)
abline(lm(P3_L~pach.d18O, data=data1_filtered))
ccf(data1_filtered$pach.d18O, data1_filtered$P3_L)


#plot rabbit p3 length and climate thru time
ylim.prim <- c(0, 3)   
ylim.sec <- c(2, 4)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

C<-ggplot(data1, aes(Cal_age, pach.d18O)) +
  geom_line(aes(y = a + P3_L*b), color = "brown") +
  geom_line(aes(y = pach.d18O), color = "blue") +
  scale_y_continuous("d18O Neogloboquadrina pachyderma", sec.axis = sec_axis(~ (. - a)/b, name = "p3 length")) +
  scale_x_reverse("Calibrated Years Before Present")+
  ggtitle("RLB Rabbit p3 Size and Climate")+
  theme_light() 
C + theme(axis.text.y.left = element_text(color="blue"),
          axis.text.y.right = element_text(color="brown"))

