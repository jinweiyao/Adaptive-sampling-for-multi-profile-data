source('Model1_simulation.R')
source('Model2_simulation.R')

#start ploting model 1 result
text = c("Proposed", "FPCA", expression("FPCA "*Delta*"= 0.01"), expression("FPCA "*Delta*"= 0.05"),expression("FPCA "*Delta*"= 0.1"))
# r = 8
par(mfcol = c(3,1))
par(mar=c(2,2,2,1))
par(oma = c(2.5,2.5,1,1))
tick<-c(0.7,1.0,1.4,2.1,2.8)
plot(x = tick,y = proposed[1:5],pch = 14,xlab = "",ylab = "", main = paste("Model ", as.roman(1), "with r = 8"),
     ylim = c(0,200),xaxt='n',cex.main = 1.5,cex = 1.2,cex.lab=1.3, cex.axis=1.3,col = rainbow(12)[8])
axis(side = 1, at = tick,cex.axis = 1.3)
lines(x = tick,y = proposed[1:5],lty = 1,col = rainbow(12)[8])
points(x = tick,y = FPCA[1:5],cex = 1.2, pch = 15,col = rainbow(12)[10])
lines(x = tick,y = FPCA[1:5],lty = 3,col = rainbow(12)[10])
points(x = tick,y = FPCA001[1:5],cex = 1.2,pch = 16)
lines(x = tick,y = FPCA001[1:5],lty = 5)
points(x = tick,y = FPCA005[1:5],cex = 1.2,pch = 17)
lines(x = tick,y = FPCA005[1:5],lty = 5)
points(x = tick,y = FPCA01[1:5],cex = 1.2, pch = 18)
lines(x = tick,y = FPCA01[1:5],lty = 5)
legend(2.0, 185, legend= text,text.width = strwidth(text)[1]*2,
       col=c(rainbow(12)[8],rainbow(12)[10],'black','black','black'), lty = c(1,3,5,5,5), pch = c(14,15,16,17,18), cex=1.2)


# r = 6
plot(x = tick,y = proposed[6:10],pch = 14,xlab = "",ylab = "", main = paste("Model ", as.roman(1), "with r = 6"),
     ylim = c(0,200),xaxt='n',cex.lab=1.3, cex = 1.2, cex.main  = 1.5, cex.axis=1.3,col = rainbow(12)[8])
axis(side = 1, at = tick,cex.axis = 1.3)
lines(x = tick,y = proposed[6:10],lty = 1,col = rainbow(12)[8])
points(x = tick,y = FPCA[6:10],cex = 1.2, pch = 15,col = rainbow(12)[10])
lines(x = tick,y = FPCA[6:10],lty = 3,col = rainbow(12)[10])
points(x = tick,y = FPCA001[6:10],cex = 1.2,pch = 16)
lines(x = tick,y = FPCA001[6:10],lty = 5)
points(x = tick,y = FPCA005[6:10],cex = 1.2,pch = 17)
lines(x = tick,y = FPCA005[6:10],lty = 5)
points(x = tick,y = FPCA01[6:10],cex = 1.2, pch = 18)
lines(x = tick,y = FPCA01[6:10],lty = 5)

# r = 4
plot(x = tick,y = proposed[11:15], cex = 1.2,pch = 14,xlab = "",ylab = "", main = paste("Model ", as.roman(1), "with r = 4"),
     ylim = c(0,200),xaxt='n',cex.lab=1.3, cex.main = 1.5,cex.axis=1.3,col = rainbow(12)[8])
axis(side = 1, at = tick,cex.axis = 1.3)
lines(x = tick,y = proposed[11:15],lty = 1,col = rainbow(12)[8])
points(x = tick,y = FPCA[11:15],cex = 1.2,pch = 15,col = rainbow(12)[10])
lines(x = tick,y = FPCA[11:15],lty = 3,col = rainbow(12)[10])
points(x = tick,y = FPCA001[11:15],cex = 1.2,pch = 16)
lines(x = tick,y = FPCA001[11:15],lty = 5)
points(x = tick,y = FPCA005[11:15], cex = 1.2,pch = 17)
lines(x = tick,y = FPCA005[11:15],lty = 5)
points(x = tick,y = FPCA01[11:15],cex = 1.2,pch = 18)
lines(x = tick,y = FPCA01[11:15],lty = 5)
mtext(expression("Shift magnitude "*delta), side = 1, line = 1, outer = TRUE, cex = 1.1 )
mtext("ARL",side = 2,line = 1,outer = TRUE,cex = 1.1,las = 0)


#start ploting model 2 result
plot(x = tick,y = poposed[1:5],pch = 14,xlab = "",ylab = "", main = paste("Model ", as.roman(2), "with r = 8"),
     ylim = c(0,200),xaxt='n',cex.main = 1.5,cex = 1.2,cex.lab=1.3, cex.axis=1.3,col = rainbow(12)[8])
axis(side = 1, at = tick,cex.axis = 1.3)
lines(x = tick,y = poposed[1:5],lty = 1,col = rainbow(12)[8])
points(x = tick,y = FCA[1:5],cex = 1.2, pch = 15,col = rainbow(12)[10])
lines(x = tick,y = FCA[1:5],lty = 3,col = rainbow(12)[10])
points(x = tick,y = FCA001[1:5],cex = 1.2,pch = 16)
lines(x = tick,y = FCA001[1:5],lty = 5)
points(x = tick,y = FCA005[1:5],cex = 1.2,pch = 17)
lines(x = tick,y = FCA005[1:5],lty = 5)
points(x = tick,y = FCA01[1:5],cex = 1.2, pch = 18)
lines(x = tick,y = FCA01[1:5],lty = 5)
legend(2.0, 190, legend= text,text.width = strwidth(text)[1]*2,
       col=c(rainbow(12)[8],rainbow(12)[10],'black','black','black'), lty = c(1,3,5,5,5), pch = c(14,15,16,17,18), cex=1.2)


# r = 6
plot(x = tick,y = poposed[6:10],pch = 14,xlab = "",ylab = "", main = paste("Model ", as.roman(2), "with r = 6"),
     ylim = c(0,200),xaxt='n',cex.lab=1.3, cex = 1.2, cex.main  = 1.5, cex.axis=1.3,col = rainbow(12)[8])
axis(side = 1, at = tick,cex.axis = 1.3)
lines(x = tick,y = poposed[6:10],lty = 1,col = rainbow(12)[8])
points(x = tick,y = FCA[6:10],cex = 1.2, pch = 15,col = rainbow(12)[10])
lines(x = tick,y = FCA[6:10],lty = 3,col = rainbow(12)[10])
points(x = tick,y = FCA001[6:10],cex = 1.2,pch = 16)
lines(x = tick,y = FCA001[6:10],lty = 5)
points(x = tick,y = FCA005[6:10],cex = 1.2,pch = 17)
lines(x = tick,y = FCA005[6:10],lty = 5)
points(x = tick,y = FCA01[6:10],cex = 1.2, pch = 18)
lines(x = tick,y = FCA01[6:10],lty = 5)

# r = 4
plot(x = tick,y = poposed[11:15], cex = 1.2,pch = 14,xlab = "",ylab = "", main = paste("Model ", as.roman(2), "with r = 4"),
     ylim = c(0,200),xaxt='n',cex.lab=1.3, cex.main = 1.5,cex.axis=1.3,col = rainbow(12)[8])
axis(side = 1, at = tick,cex.axis = 1.3)
lines(x = tick,y = poposed[11:15],lty = 1,col = rainbow(12)[8])
points(x = tick,y = FCA[11:15],cex = 1.2,pch = 15,col = rainbow(12)[10])
lines(x = tick,y = FCA[11:15],lty = 3,col = rainbow(12)[10])
points(x = tick,y = FCA001[11:15],cex = 1.2,pch = 16)
lines(x = tick,y = FCA001[11:15],lty = 5)
points(x = tick,y = FCA005[11:15], cex = 1.2,pch = 17)
lines(x = tick,y = FCA005[11:15],lty = 5)
points(x = tick,y = FCA01[11:15],cex = 1.2,pch = 18)
lines(x = tick,y = FCA01[11:15],lty = 5)
mtext(expression("Shift magnitude "*delta), side = 1, line = 1, outer = TRUE, cex = 1.1 )
mtext("ARL",side = 2,line = 1,outer = TRUE,cex = 1.1,las = 0)

out_se #model1 ARL
out_se1 #model2 ARL



