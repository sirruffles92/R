Mouse.Experiments.Spring.2018 <- read.csv("~/Documents/Utsav Mouse Experiments/Mouse Experiments Spring 2018.csv")

X07.03.18_Measures <- do.call(rbind, strsplit(as.character(Mouse.Experiments.Spring.2018$X07.03.18), 'X'))
X07.03.18_Edge_1 <- as.numeric(X07.03.18_Measures[,1])
X07.03.18_Edge_2 <- as.numeric(X07.03.18_Measures[,2])
X07.03.18_Volume <- c()

for (i in 1:length(X07.03.18_Edge_1)){
  X07.03.18_Volume[i] <- (max(X07.03.18_Edge_1[i], X07.03.18_Edge_2[i]) * 0.5) * (min(X07.03.18_Edge_1[i], X07.03.18_Edge_2[i]) ^2)
}

X07.13.18_Measures <- do.call(rbind, strsplit(as.character(Mouse.Experiments.Spring.2018$X07.13.18), 'X'))
X07.13.18_Edge_1 <- as.numeric(X07.13.18_Measures[,1])
X07.13.18_Edge_2 <- as.numeric(X07.13.18_Measures[,2])
X07.13.18_Volume <- c()

for (i in 1:length(X07.13.18_Edge_1)){
  X07.13.18_Volume[i] <- (max(X07.13.18_Edge_1[i], X07.13.18_Edge_2[i]) * 0.5) * (min(X07.13.18_Edge_1[i], X07.13.18_Edge_2[i]) ^2)
}

X07.20.18_Measures <- do.call(rbind, strsplit(as.character(Mouse.Experiments.Spring.2018$X07.20.18), 'X'))
X07.20.18_Edge_1 <- as.numeric(X07.20.18_Measures[,1])
X07.20.18_Edge_2 <- as.numeric(X07.20.18_Measures[,2])
X07.20.18_Volume <- c()

for (i in 1:length(X07.20.18_Edge_1)){
  X07.20.18_Volume[i] <- (max(X07.20.18_Edge_1[i], X07.20.18_Edge_2[i]) * 0.5) * (min(X07.20.18_Edge_1[i], X07.20.18_Edge_2[i]) ^2)
}

X07.25.18_Measures <- do.call(rbind, strsplit(as.character(Mouse.Experiments.Spring.2018$X07.25.18), 'X'))
X07.25.18_Edge_1 <- as.numeric(X07.25.18_Measures[,1])
X07.25.18_Edge_2 <- as.numeric(X07.25.18_Measures[,2])
X07.25.18_Volume <- c()

for (i in 1:length(X07.25.18_Edge_1)){
  X07.25.18_Volume[i] <- (max(X07.25.18_Edge_1[i], X07.25.18_Edge_2[i]) * 0.5) * (min(X07.25.18_Edge_1[i], X07.25.18_Edge_2[i]) ^2)
}

X08.03.18_Measures <- do.call(rbind, strsplit(as.character(Mouse.Experiments.Spring.2018$X08.03.18), 'X'))
X08.03.18_Edge_1 <- as.numeric(X08.03.18_Measures[,1])
X08.03.18_Edge_2 <- as.numeric(X08.03.18_Measures[,2])
X08.03.18_Volume <- c()

for (i in 1:length(X07.03.18_Edge_1)){
  X08.03.18_Volume[i] <- (max(X08.03.18_Edge_1[i], X08.03.18_Edge_2[i]) * 0.5) * (min(X08.03.18_Edge_1[i], X08.03.18_Edge_2[i]) ^2)
}


X08.06.18_Measures <- do.call(rbind, strsplit(as.character(Mouse.Experiments.Spring.2018$Tumor.at.Sac), 'X'))

X08.06.18_Edge_1 <- as.numeric(X08.06.18_Measures[,1])
X08.06.18_Edge_2 <- as.numeric(X08.06.18_Measures[,2])
X08.06.18_Volume <- c()

for (i in 1:length(X08.06.18_Edge_1)){
  X08.06.18_Volume[i] <- (max(X08.06.18_Edge_1[i], X08.06.18_Edge_2[i]) * 0.5) * (min(X08.06.18_Edge_1[i], X08.06.18_Edge_2[i]) ^2)
}



Tumor_Data <- data.frame(Mouse.Experiments.Spring.2018$Group, X07.03.18_Volume, X07.13.18_Volume, X07.20.18_Volume, X07.25.18_Volume, X08.03.18_Volume)

days <- c(0, 21, 31, 38, 43, 52, 55)

Tumor_Volumes_CAF23_MFP9_Average <- (subset(Modified_Tumor_Data$Tumor.Volumes, Modified_Tumor_Data$Mouse.Experiments.Spring.2018.Group == 'CAF23 MFP9'))
Tumor_Volumes_CAF23_MFP9_SEM <- c(0,0,0,0,0,0,0)
Tumor_Volumes_CAF23_MFP9_Average <- c(0, mean(Tumor_Volumes_CAF23_MFP9_Average[1:5]),mean(Tumor_Volumes_CAF23_MFP9_Average[6:10]),mean(Tumor_Volumes_CAF23_MFP9_Average[11:15]),
                                      mean(Tumor_Volumes_CAF23_MFP9_Average[16:20]),mean(Tumor_Volumes_CAF23_MFP9_Average[21:25]),mean(Tumor_Volumes_CAF23_MFP9_Average[26:30]))

Tumor_Volumes_CAF23_DT28_Average <- (subset(Modified_Tumor_Data$Tumor.Volumes, Modified_Tumor_Data$Mouse.Experiments.Spring.2018.Group == 'CAF23+DT28'))

Tumor_Volumes_CAF23_DT28_SEM <- c(0, sd(Tumor_Volumes_CAF23_DT28_Average[1:10]) / sqrt(10),sd(Tumor_Volumes_CAF23_DT28_Average[c(12:15,17:20)]) / sqrt(8),
                                  sd(Tumor_Volumes_CAF23_DT28_Average[c(22:25,27:30)]) / sqrt(8),sd(Tumor_Volumes_CAF23_DT28_Average[c(32:35,37:40)]) / sqrt(8),
                                  sd(Tumor_Volumes_CAF23_DT28_Average[c(42:45,47:50)]) / sqrt(8), sd(Tumor_Volumes_CAF23_DT28_Average[c(52:55,57:60)]) / sqrt(8))

Tumor_Volumes_CAF23_DT28_Average <- c(0, mean(Tumor_Volumes_CAF23_DT28_Average[1:10]),mean(Tumor_Volumes_CAF23_DT28_Average[c(12:15,17:20)]),
                                      mean(Tumor_Volumes_CAF23_DT28_Average[c(22:25,27:30)]),mean(Tumor_Volumes_CAF23_DT28_Average[c(32:35,37:40)]),
                                      mean(Tumor_Volumes_CAF23_DT28_Average[c(42:45,47:50)]), mean(Tumor_Volumes_CAF23_DT28_Average[c(52:55,57:60)]))

Tumor_Volumes_DT28_MFP9_CAF23_MFP4_Average <- (subset(Modified_Tumor_Data$Tumor.Volumes, Modified_Tumor_Data$Mouse.Experiments.Spring.2018.Group == ' DT28 MFP9 CAF23 MFP4'))

Tumor_Volumes_DT28_MFP9_CAF23_MFP4_SEM <- c(0, sd(Tumor_Volumes_DT28_MFP9_CAF23_MFP4_Average[c(2,4:5,7:10)]) / sqrt(7),sd(Tumor_Volumes_DT28_MFP9_CAF23_MFP4_Average[c(12,14:15,17:20)]) / sqrt(7),
                                  sd(Tumor_Volumes_DT28_MFP9_CAF23_MFP4_Average[c(22,24:25,27:30)]) / sqrt(7),sd(Tumor_Volumes_DT28_MFP9_CAF23_MFP4_Average[c(32,34:35,37:40)]) / sqrt(7),
                                  sd(Tumor_Volumes_DT28_MFP9_CAF23_MFP4_Average[c(42,44:45,47:50)]) / sqrt(7),sd(Tumor_Volumes_DT28_MFP9_CAF23_MFP4_Average[c(52,54:55,57:60)]) / sqrt(7))
Tumor_Volumes_DT28_MFP9_CAF23_MFP4_Average <- c(0, mean(Tumor_Volumes_DT28_MFP9_CAF23_MFP4_Average[c(2,4:5,7:10)]),mean(Tumor_Volumes_DT28_MFP9_CAF23_MFP4_Average[c(12,14:15,17:20)]),
                                                mean(Tumor_Volumes_DT28_MFP9_CAF23_MFP4_Average[c(22,24:25,27:30)]),mean(Tumor_Volumes_DT28_MFP9_CAF23_MFP4_Average[c(32,34:35,37:40)]),
                                                mean(Tumor_Volumes_DT28_MFP9_CAF23_MFP4_Average[c(42,44:45,47:50)]),mean(Tumor_Volumes_DT28_MFP9_CAF23_MFP4_Average[c(52,54:55,57:60)]))

Tumor_Volumes_DT28_MFP9_Average <- (subset(Modified_Tumor_Data$Tumor.Volumes, Modified_Tumor_Data$Mouse.Experiments.Spring.2018.Group == 'DT28 MFP9'))

Tumor_Volumes_DT28_MFP9_SEM <- c(0, sd(Tumor_Volumes_DT28_MFP9_Average[1:10]) / sqrt(10),sd(Tumor_Volumes_DT28_MFP9_Average[c(12:15,17:20)]) / sqrt(8),
                                            sd(Tumor_Volumes_DT28_MFP9_Average[c(22:25,27:30)]) / sqrt(8),sd(Tumor_Volumes_DT28_MFP9_Average[c(32:35,37:40)]) / sqrt(8),
                                            sd(Tumor_Volumes_DT28_MFP9_Average[c(42:45,47:50)]) / sqrt(8),sd(Tumor_Volumes_DT28_MFP9_Average[c(52:55,57:60)]) / sqrt(8))
Tumor_Volumes_DT28_MFP9_Average <- c(0, mean(Tumor_Volumes_DT28_MFP9_Average[1:10]),mean(Tumor_Volumes_DT28_MFP9_Average[c(12:15,17:20)]),mean(Tumor_Volumes_DT28_MFP9_Average[c(22:25,27:30)]),
                                      mean(Tumor_Volumes_DT28_MFP9_Average[c(32:35,37:40)]),mean(Tumor_Volumes_DT28_MFP9_Average[c(42:45,47:50)]), mean(Tumor_Volumes_DT28_MFP9_Average[c(52:55,57:60)]))



plot(days, Tumor_Volumes_CAF23_MFP9_Average, ylim=c(0,350),xlim=c(0,60),xlab="Days Post Injection", ylab="Tumor Volume (mm^3)", title("Opposite Flank Injection Experiment"))
lines(days, Tumor_Volumes_CAF23_MFP9_Average)
lines(days, Tumor_Volumes_CAF23_DT28_Average, type='o', col='red',pch=17)
arrows(days, Tumor_Volumes_CAF23_DT28_Average + Tumor_Volumes_CAF23_DT28_SEM, days, Tumor_Volumes_CAF23_DT28_Average - Tumor_Volumes_CAF23_DT28_SEM, 
       length=0.05, angle=90, code=3,col="red")
lines(days, Tumor_Volumes_DT28_MFP9_CAF23_MFP4_Average, type='o', col='green',pch=18)
arrows(days, Tumor_Volumes_DT28_MFP9_CAF23_MFP4_Average + Tumor_Volumes_DT28_MFP9_CAF23_MFP4_SEM, days, Tumor_Volumes_DT28_MFP9_CAF23_MFP4_Average - Tumor_Volumes_DT28_MFP9_CAF23_MFP4_SEM, 
       length=0.05, angle=90, code=3,col="green")
lines(days, Tumor_Volumes_DT28_MFP9_Average, type='o', col='blue',pch=20)
arrows(days, Tumor_Volumes_DT28_MFP9_Average + Tumor_Volumes_DT28_MFP9_SEM, days, Tumor_Volumes_DT28_MFP9_Average - Tumor_Volumes_DT28_MFP9_SEM, 
       length=0.05, angle=90, code=3,col="blue")














