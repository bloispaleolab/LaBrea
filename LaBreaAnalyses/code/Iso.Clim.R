options(scipen=999)
library(ggplot2)
library(tidyverse)
library(scales)
library(RColorBrewer)

iso_dat<-read.csv("data/processed/Interpolation/IsoClim_interpolated.csv")# All rabbits and Squirrels with mean calibrated ages < 50 ka BP

####Intraspecific diet ###

## Compare desert cottontails to brush rabbits

iso_rabbits<- iso_dat[grep("S. audubonii|S. bachmani", iso_dat$Taxon), c(3:5)]
rownames(iso_rabbits) <- NULL
levels(iso_rabbits$Taxon)
iso_rabbits$Taxon<-droplevels(iso_rabbits$Taxon)# iso_rabbits$Taxon<-factor(iso_rabbits$Taxon)
levels(iso_rabbits$Taxon)

#carbon
boxplot(d13C~Taxon, data=iso_rabbits, main="Rabbit diet niche", 
        xlab="Taxon", ylab="d13C", col = "purple")
wilcox.test(d13C~Taxon, data=iso_rabbits, correct = TRUE)

#nitrogen
boxplot(d15N~Taxon,data=iso_rabbits, main="Rabbit diet niche", 
        xlab="Species", ylab="d15N", col = "purple")
wilcox.test(d15N~Taxon, data=iso_rabbits, correct = TRUE)


hull_species <- iso_rabbits %>%
  group_by(Taxon) %>%
  slice(chull(d13C, d15N))

p<-ggplot(iso_rabbits, aes(x = d13C, y = d15N)) +
  geom_point()+
  labs(x= expression({delta}^13*C~'\u2030'), y= expression({delta}^15*N~'\u2030'))

p + aes(fill = factor(Taxon)) + 
  geom_polygon(data = hull_species, alpha = 0.5)+
  scale_fill_discrete(name = "Taxon", 
                      labels = c("Desert Cottontail", "Brush Rabbit"))


## compare all rabbits to squirrels
iso_dat_2<-read.csv("data/processed/Interpolation/IsoClim_interpolated.csv")
levels(iso_dat_2$Taxon) 
levels(iso_dat_2$Taxon)[levels(iso_dat_2$Taxon)=="S. audubonii"] <- "Sylvilagus"
levels(iso_dat_2$Taxon)[levels(iso_dat_2$Taxon)=="S. bachmani"] <- "Sylvilagus"
levels(iso_dat_2$Taxon)[levels(iso_dat_2$Taxon)=="Leporidae"] <- "Sylvilagus"
levels(iso_dat_2$Taxon)[levels(iso_dat_2$Taxon)=="Sylvilagus sp"] <- "Sylvilagus"
levels(iso_dat_2$Taxon)

#carbon
boxplot(d13C~Taxon,data=iso_dat_2, main="Rabbit and Squirrel diet niche", 
        xlab="Species", ylab="d13C", col = "purple")

wilcox.test(d13C~Taxon, data=iso_dat_2, correct = TRUE)

#nitrogen
boxplot(d15N~Taxon,data=iso_dat_2, main="Rabbit and Squirrel diet niche", 
        xlab="Species", ylab="d15N", col = "purple")

wilcox.test(d15N~Taxon, data=iso_dat_2, correct = TRUE)


hull_species <- iso_dat_2 %>%
  group_by(Taxon) %>%
  slice(chull(d13C, d15N))

p<-ggplot(iso_dat_2, aes(x = d13C, y = d15N)) +
  geom_point()+
  labs(x= expression({delta}^13*C~'\u2030'), y= expression({delta}^15*N~'\u2030'))

p + aes(fill = factor(Taxon)) + 
  geom_polygon(data = hull_species, alpha = 0.5)+
  scale_fill_discrete(name = "Taxon", 
                      labels = c("Rabbit", "Squirrel"))



#### Isotopes and Climate ###

data0<-read.csv("data/processed/Interpolation/IsoClim_interpolated.csv")# All rabbits and Squirrels with mean calibrated ages < 50 ka BP

### all taxa ###

#carbon and climate
plot(d13C~pach.d18O, data=data0, pch=16)
carbon.lm<-lm(d13C~pach.d18O, data=data0)
summary(carbon.lm)
abline(lm(d13C~pach.d18O, data=data0))
ccf(data0$pach.d18O, data0$d13C)

#nitrogen and climate
plot(d15N~pach.d18O, data=data0, pch=16)
nitrogen.lm<-lm(d15N~pach.d18O, data=data0)
summary(nitrogen.lm)
abline(lm(d15N~pach.d18O, data=data0))
ccf(data0$pach.d18O, data0$d15N)

#plot carbon and climate thru time
data0<-read.csv("P23_IsoClim_interpolated.csv")
ylim.prim <- c(0, 3)   
ylim.sec <- c(-23, -17)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

C<-ggplot(data0, aes(Cal_age, pach.d18O)) +
  geom_line(aes(y = a + d13C*b), color = "brown") +
  geom_line(aes(y = pach.d18O), color = "blue") +
  scale_y_reverse("d18O Neogloboquadrina pachyderma", sec.axis = sec_axis(~ (. - a)/b, name = expression({delta}^13*C~'\u2030'))) +
  scale_x_reverse("Calibrated Years Before Present")+
  ggtitle("RLB Small Mammal Carbon and Climate")+
  theme_light() 
C + theme(axis.text.y.left = element_text(color="blue"),
          axis.text.y.right = element_text(color="brown"))
 
#plot nitrogen thru time
ylim.prim <- c(0, 3)   
ylim.sec <- c(3, 14)

b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

p<-ggplot(data0, aes(Cal_age, pach.d18O)) +
  geom_line(aes(y = a + d15N*b), color = "red") +
  geom_line(aes(y = pach.d18O), color = "blue") +
  scale_y_reverse("d18O Neogloboquadrina pachyderma", sec.axis = sec_axis(~ (. - a)/b, name = expression({delta}^15*N~'\u2030'))) +
  scale_x_reverse("Calibrated Years Before Present")+
  ggtitle("RLB Small Mammal Nitrogen and Climate")+
  theme_light() 
p + theme(axis.text.y.left = element_text(color="blue"),
          axis.text.y.right = element_text(color="red"))


### Just squirrles ###

#carbon and climate
plot(d13C[Taxon=='O. beecheyi']~pach.d18O[Taxon=='O. beecheyi'], data=data0, pch=16)
carbon.lm<-lm(d13C[Taxon=='O. beecheyi']~pach.d18O[Taxon=='O. beecheyi'], data=data0)
summary(carbon.lm)
abline(lm(d13C[Taxon=='O. beecheyi']~pach.d18O[Taxon=='O. beecheyi'], data=data0))
#ccf(data0$pach.d18O[Taxon=='O. beecheyi'], data0$d13C[Taxon=='O. beecheyi'])

#nitrogen and climate
plot(d15N[Taxon=='O. beecheyi']~pach.d18O[Taxon=='O. beecheyi'], data=data0, pch=16)
carbon.lm<-lm(d15N[Taxon=='O. beecheyi']~pach.d18O[Taxon=='O. beecheyi'], data=data0)
summary(carbon.lm)
abline(lm(d15N[Taxon=='O. beecheyi']~pach.d18O[Taxon=='O. beecheyi'], data=data0))
#ccf(data0$pach.d18O[Taxon=='O. beecheyi'], data0$d13C[Taxon=='O. beecheyi'])

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
plot(d13C[Taxon=='Sylvilagus']~pach.d18O[Taxon=='Sylvilagus'], data=iso_dat_2, pch=16)
carbon.lm<-lm(d13C[Taxon=='Sylvilagus']~pach.d18O[Taxon=='Sylvilagus'], data=iso_dat_2)
summary(carbon.lm)
abline(lm(d13C[Taxon=='Sylvilagus']~pach.d18O[Taxon=='Sylvilagus'], data=iso_dat_2))

#nitrogen and climate
plot(d15N[Taxon=='Sylvilagus']~pach.d18O[Taxon=='Sylvilagus'], data=iso_dat_2, pch=16)
carbon.lm<-lm(d15N[Taxon=='Sylvilagus']~pach.d18O[Taxon=='Sylvilagus'], data=iso_dat_2)
summary(carbon.lm)
abline(lm(d15N[Taxon=='Sylvilagus']~pach.d18O[Taxon=='Sylvilagus'], data=iso_dat_2))

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

plot(P3_L~pach.d18O, data=data1, xlim=rev(c(0, 3)), pch=16)
size.lm<-lm(P3_L~pach.d18O, data=data1)
summary(size.lm)
abline(lm(P3_L~pach.d18O, data=data1))
ccf(data1$pach.d18O, data1$P3_L)

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

