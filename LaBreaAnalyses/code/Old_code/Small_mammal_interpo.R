# Interpolate rabbit isotopes to calibrated radiocarbon ages of small mammals 

rabbit<- read.delim("data/processed/Interpolation/Rabbit_isotopes.txt")
ages <- read.delim("data/processed/Interpolation/Calibrated_dates.txt") # mean, error, 

xout <- ages$Calibrated.age
rabbit_extracted <- approx(x=rabbit$Calibrated.age, y=rabbit$del13C_permil, method="linear", xout=xout)
rabbit_extracted$y # interpolated d18O at each radiocarbon date

rabbitN_extracted <- approx(x=rabbit$Calibrated.age, y=rabbit$del15N_permil, method="linear", xout=xout)
rabbitN_extracted$y 

rabbit_interpolated <-cbind(rabbit_extracted$y, rabbitN_extracted$y)

write.csv(rabbit_interpolated, file = "data/processed/Interpolation/rabbit_interpolated.csv")


# Interpolate squirrel isotopes to calibrated radiocarbon ages of small mammals 

squirrel<- read.delim("data/processed/Interpolation/Squirrel_isotopes.txt")
ages <- read.delim("data/processed/Interpolation/Calibrated_dates.txt") # mean, error, 

xout <- ages$Calibrated.age
squirrel_extracted <- approx(x=squirrel$Calibrated.age, y=squirrel$del13C_permil, method="linear", xout=xout)
squirrel_extracted$y # interpolated d18O at each radiocarbon date

squirrelN_extracted <- approx(x=squirrel$Calibrated.age, y=squirrel$del15N_permil, method="linear", xout=xout)
squirrelN_extracted$y 

squirrel_interpolated <-cbind(squirrel_extracted$y, squirrelN_extracted$y)

write.csv(squirrel_interpolated, file = "data/processed/Interpolation/squirrel_interpolated.csv")

