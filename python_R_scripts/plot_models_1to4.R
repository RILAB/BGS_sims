#Note that demographic model names were changed during the course of manuscript writing so that the original input file names and model names therein
#no longer correspond to the demographic model names used in the manuscript. Namely, demographic models 2-9 were re-labeled as models 5-12
#and demographic models 10-12 were re-labeled as models 2-4. This script makes the necessary changes to fit model names to those used for the
#final manuscript, but it still assumes that the input file names and model names therein correspond to the original model names.

args = commandArgs(TRUE)

print(paste('data directory: ', args[1], sep=' '))

#set working directory to directory where files for observed pi/pi0 and xi/xi0 live
setwd(args[1])

#observed pi/pi0 from simulations and bootstrapped 95 percent CIs (script assumes you are in this working directory)
mean_pi_file = 'ratio_mean_pi_generic_0.00082cM.txt' #mean pi/pi0
lower_pi_file = 'ratio_95CI_lower_pi_generic_0.00082cM.txt' #lower 95CI pi/pi0
upper_pi_file = 'ratio_95CI_upper_pi_generic_0.00082cM.txt' #upper 95CI pi/pi0
mean_xi_file ='ratio_mean_xi_generic_0.00082cM.txt' #mean xi/xi0
lower_xi_file = 'ratio_95CI_lower_xi_generic_0.00082cM.txt' #lower 95CI xi/xi0
upper_xi_file = 'ratio_95CI_upper_xi_generic_0.00082cM.txt' #upper 95CI xi/xi0

#Nordborg 1996 pi/pi0 expectations for demographic models
nord96_models10to12 = 'gen_models1_and_10to12_nord96_out.txt'

###########################
##start plot demographies##
###########################
xpos=c(2, seq(12, 102, 10))
tck=-0.075
ymin=90
ymax=1000000
yseq=c(100, 1000, 10000, 100000, 1000000)

#demographic model information
demmodels = read.csv('new_generic_models_0.01gen_models1to5_and_models10to12.txt', header=TRUE)
colnames(demmodels) = c('y1', 'y2', 'y3', 'y4', 'y5', 'y10', 'y11', 'y12', 'gen', 'X')

pdf('models1to4.pdf', height=6, width=5.825)
par(mfrow=c(3,2), oma=c(7,10.2,0,4) + 0.1, mar=c(0,0,0.9,0.9), las=1)

#set an lwd variable for line widths
w = 3

legendcex = 1.15

#make plot for demographic models 1, 2, and 3
plot(demmodels$X, demmodels$y1, ylim=c(ymin, ymax), type='n', xaxt='n', cex=2, log='y', yaxt='n', ylab='')

#make grid lines
abline(h=c(yseq), v=c(xpos), col = "darkgray", lty = "dotted")

axis(2, at=yseq, labels=c('100', '1000', '10000', '100000', '1000000'), cex.axis=1.25)
axis(1, at=xpos, tck=tck, labels=FALSE)

#lines for demographic model population sizes through time
lines(demmodels$X, demmodels$y1, lty=1, lwd=w) #model 1
lines(demmodels$X, demmodels$y10, lty=1, lwd=w, col='red') #model 2
lines(demmodels$X, demmodels$y11, lty=1, lwd=w, col='blue') #model 3

#make legend boxes to identify models and their colors and to identify the 'Expected' pi/pi0 dotted line from Nordborg 1996 pi/pi0 expectations
legend('bottomright', legend=c('Expectedd'), col='white', lty=-1, lwd=w, box.lwd=2, cex=0.9, bg='white', text.col='white') #this is just to make a blank box that's a tad smaller than the default, it will not print
legend('bottomright', legend=c('Expected'), col='black', lty=3, lwd=w, box.lwd=NA, cex=0.95, text.width=27) #print 'Expected' in the blank legend box
legend('topleft', legend=c('model 1', 'model 2', 'model 3'), col=c('black', 'red', 'blue'), lty=c(1, 1, 1), ncol=1, lwd=w, box.lwd=2, cex=legendcex, bg='white')

text(-65, 10000, "Effective Population Size", srt=90, xpd=NA, cex=1.35)
text(-54, 10000, expression(paste("(", italic("N"[e]), "; log-scale)", sep='')), srt=90, xpd=NA, cex=1.35)
box(lwd=2)

#make plot for demographic model 1 and 4
plot(demmodels$X, demmodels$y1, ylim=c(9800, 100000), type='n', xaxt='n', yaxt='n', cex=2, log='y')

#make grid lines
abline(h=c(10000, 20000, 40000, 100000), v=c(xpos), col = "darkgray", lty = "dotted")

axis(1, at=xpos, tck=tck, labels=FALSE)
axis(4, at=c(10000, 20000, 40000, 100000), labels=c('10000', '20000', '40000', '100000'), cex.axis=1.25)

#lines for demographic model population sizes through time
lines(demmodels$X, demmodels$y1, lty=1, lwd=w) #model 1
lines(demmodels$X, demmodels$y12, lty=1, lwd=w, col='purple') #model 4

#make legend boxes to identify models and their colors and to identify the 'Expected' pi/pi0 dotted line from Nordborg 1996 pi/pi0 expectations
legend('bottomright', legend=c('Expectedd'), col='white', lty=-1, lwd=w, box.lwd=2, cex=0.9, bg='white', text.col='white') #this is just to make a blank box that's a tad smaller than the default, it will not print
legend('bottomright', legend=c('Expected'), col='black', lty=3, lwd=w, box.lwd=NA, cex=0.95, text.width=27) #print 'Expected' in the blank legend box
legend('topleft', legend=c('model 1', 'model 4'), col=c('black', 'purple'), lty=c(1, 1), ncol=1, lwd=w, box.lwd=2, cex=legendcex, bg='white')
box(lwd=2)

#########################
##end plot demographies##
#########################

#make function to plot pi/pi0 and xi/xi0 for each generation in the simulation model 
plot_func = function(meanfile, CIlowerfile, CIupperfile, YMIN1, YMAX1, BY1, YMIN2, YMAX2, BY2, LWD, lowerlabel=FALSE, PIS=FALSE){

	#make vector for positions along x-axis for plotting pi/pi0 and xi/xi0 for each generation (404 total)
	x = c(0.8, 2, 3.2, 4, 4.8, 6, 7.2, 8, 8.8, 10, 11.2, 12, 12.8, 14, 15.2, 16, 16.8, 18, 19.2, 20, 20.8, 22, 23.2, 24, 24.8, 26, 27.2, 28, 28.8, 30, 31.2, 32, 32.8, 34, 35.2, 36, 36.8, 38, 39.2, 40, 40.8, 42, 43.2, 44, 44.8, 46, 47.2, 48, 48.8, 50, 51.2, 52, 52.8, 54, 55.2, 56, 56.8, 58, 59.2, 60, 60.8, 62, 63.2, 64, 64.8, 66, 67.2, 68, 68.8, 70, 71.2, 72, 72.8, 74, 75.2, 76, 76.8, 78, 79.2, 80, 80.8, 82, 83.2, 84, 84.8, 86, 87.2, 88, 88.8, 90, 91.2, 92, 92.8, 94, 95.2, 96, 96.8, 98, 99.2, 100, 100.8, 102, 103.2, 104, 104.8, 106, 107.2, 108, 108.8, 110, 111.2, 112, 112.8, 114, 115.2, 116, 116.8, 118, 119.2, 120, 120.8, 122, 123.2, 124, 124.8, 126, 127.2, 128, 128.8, 130, 131.2, 132, 132.8, 134, 135.2, 136, 136.8, 138, 139.2, 140, 140.8, 142, 143.2, 144, 144.8, 146, 147.2, 148, 148.8, 150, 151.2, 152, 152.8, 154, 155.2, 156, 156.8, 158, 159.2, 160, 160.8, 162, 163.2, 164, 164.8, 166, 167.2, 168, 168.8, 170, 171.2, 172, 172.8, 174, 175.2, 176, 176.8, 178, 179.2, 180, 180.8, 182, 183.2, 184, 184.8, 186, 187.2, 188, 188.8, 190, 191.2, 192, 192.8, 194, 195.2, 196, 196.8, 198, 199.2, 200, 200.8, 202, 203.2, 204, 204.8, 206, 206.8, 208, 208.8, 210, 211.2, 212, 212.8, 214, 215.2, 216, 216.8, 218, 219.2, 220, 220.8, 222, 223.2, 224, 224.8, 226, 227.2, 228, 228.8, 230, 231.2, 232, 232.8, 234, 235.2, 236, 236.8, 238, 239.2, 240, 240.8, 242, 243.2, 244, 244.8, 246, 247.2, 248, 248.8, 250, 251.2, 252, 252.8, 254, 255.2, 256, 256.8, 258, 259.2, 260, 260.8, 262, 263.2, 264, 264.8, 266, 267.2, 268, 268.8, 270, 271.2, 272, 272.8, 274, 275.2, 276, 276.8, 278, 279.2, 280, 280.8, 282, 283.2, 284, 284.8, 286, 287.2, 288, 288.8, 290, 291.2, 292, 292.8, 294, 295.2, 296, 296.8, 298, 299.2, 300, 300.8, 302, 303.2, 304, 304.8, 306, 307.2, 308, 308.8, 310, 311.2, 312, 312.8, 314, 315.2, 316, 316.8, 318, 319.2, 320, 320.8, 322, 323.2, 324, 324.8, 326, 327.2, 328, 328.8, 330, 331.2, 332, 332.8, 334, 335.2, 336, 336.8, 338, 339.2, 340, 340.8, 342, 343.2, 344, 344.8, 346, 347.2, 348, 348.8, 350, 351.2, 352, 352.8, 354, 355.2, 356, 356.8, 358, 359.2, 360, 360.8, 362, 363.2, 364, 364.8, 366, 367.2, 368, 368.8, 370, 371.2, 372, 372.8, 374, 375.2, 376, 376.8, 378, 379.2, 380, 380.8, 382, 383.2, 384, 384.8, 386, 387.2, 388, 388.8, 390, 391.2, 392, 392.8, 394, 395.2, 396, 396.8, 398, 399.2, 400, 400.8, 402, 403.2, 404)
	
	tck=-0.075
	
	#y-axis limits and grid points/tick marks for plotting models 1, 2, and 3 
	yseq1=seq(YMIN1, YMAX1, BY1)
	ylab1=as.character(yseq1)
	
	#y-axis limits and grid points/tick marks for plotting models 1 and 4 
	yseq2=seq(YMIN2, YMAX2, BY2)
	ylab2=as.character(yseq2)
	
	#sub function for plotting actual pi/pi0 lines and the 95 percent CI envelopes
	subplot_CI = function(YCONST, X, YDEM1, YDEM2, YCONST_CIlower, YCONST_CIupper, YDEM1_CIlower, YDEM1_CIupper, YDEM2_CIlower, YDEM2_CIupper, XLABEL=lowerlabel, YLABEL=FALSE, XLAB, XPOS, YTYPE=1){
		
		#if YTYPE == 1, use ylimits for plotting models 1, 2, and 3
		#if YTYPE == 2, use ylimits for plotting models 1 and 4
		if(YTYPE == 1){
			YMIN = YMIN1
			YMAX = YMAX1
			yseq = yseq1
			ylab = ylab1
		} else{
			YMIN = YMIN2
			YMAX = YMAX2
			yseq = yseq2
			ylab = ylab2		
		}
		
		#start empty plot
		plot(X, YCONST, ylim=c(YMIN, YMAX), type='n', xaxt='n', yaxt='n', cex=2)
		
		#make grid points
		abline(h=c(yseq), v=c(XPOS), col = "darkgray", lty = "dotted")
		
		#if YLABEL==TRUE, then we are plotting first column of plots so we need tick marks and a y-axis title
		if(YLABEL==TRUE){
			axis(2, at=yseq, labels=ylab, cex.axis=1.25)
			FRACTITLEY1 = 0.525
			FRACTITLEY2 = 0.425
			if(PIS==TRUE){
				text(par("usr")[1]-(FRACTITLEY1*(par("usr")[2]-par("usr")[1])), (YMIN+YMAX)/2, "Relative Diversity", srt=90, xpd=NA, cex=1.5)
				text(par("usr")[1]-(FRACTITLEY2*(par("usr")[2]-par("usr")[1])), (YMIN+YMAX)/2, expression(paste('(', pi, '/', pi[0], ')', sep='')), srt=90, xpd=NA, cex=1.5)
			} else{
				text(par("usr")[1]-(FRACTITLEY1*(par("usr")[2]-par("usr")[1])), (YMIN+YMAX)/2, "Relative Singleton", srt=90, xpd=NA, cex=1.5)
				text(par("usr")[1]-(FRACTITLEY2*(par("usr")[2]-par("usr")[1])), (YMIN+YMAX)/2, expression(paste('Density (', xi, '/', xi[0], ')', sep='')), srt=90, xpd=NA, cex=1.5)				
			}
		} else{ #YLABEL==FALSE, so we are plotting second column of plots, so only need tick marks
			axis(4, at=yseq, cex.axis=1.25)
		}
		
		#if XLABEL==TRUE, then we are plotting bottom row of plots so we need tick marks and an x-axis title
		if(XLABEL == TRUE){
			axis(1, at=XPOS, labels=FALSE)
			FRACX = 0.1*(par("usr")[4]-par("usr")[3])
			FRACTITLEX = 0.4*(par("usr")[4]-par("usr")[3])
			text(XPOS, par("usr")[3] - FRACX, labels=XLAB, srt=45, adj=1, pos=NULL, xpd=NA, cex=1.25)
			text(202, par("usr")[3] - FRACTITLEX, expression(paste('Time (', italic('N'[anc]), ' generations)', sep='')), cex=1.35, xpd=NA)
		} else{ #XLABEL==FALSE, so we are plotting second row of plots, so only need tick marks (but need longer tick marks to connect to bottom row of plots)
			axis(1, at=XPOS, tck=tck, labels=FALSE)
		}
		
		#add envelopes representing 95% CIs for model 1
		polygon(c(X, rev(X)), c(YCONST_CIlower, rev(YCONST_CIupper)), col=rgb(0, 0, 0, 0.4), border=NA)
		#add lines for pi/pi0 (or xi/xi0) for model 1
		lines(X, YCONST, lty=1, lwd=LWD)
		
		#add envelopes representing 95% CIs for other models and lines for pi/pi0 (or xi/xi0)
		if(YDEM1 == YDEM2){# model 4
			polygon(c(X, rev(X)), c(YDEM2_CIlower, rev(YDEM2_CIupper)), col=rgb(160/255, 32/255, 240/255, 0.4), border=NA)
			lines(X, YDEM2, lty=1, lwd=LWD, col='purple')			
		} else{# models 2 and 3
			polygon(c(X, rev(X)), c(YDEM2_CIlower, rev(YDEM2_CIupper)), col=rgb(0, 0, 255/255, 0.4), border=NA)
			lines(X, YDEM2, lty=1, lwd=LWD, col='blue')
			polygon(c(X, rev(X)), c(YDEM1_CIlower, rev(YDEM1_CIupper)), col=rgb(255/255, 0, 0, 0.4), border=NA)
			lines(X, YDEM1, lty=1, lwd=LWD, col='red')
		}
		box(lwd=2)
	}
	
	df_mean = read.table(meanfile, header=FALSE)[,seq(2,405, 1)]
	df_lower = read.table(CIlowerfile, header=FALSE)[,seq(2,405, 1)]
	df_upper = read.table(CIupperfile, header=FALSE)[,seq(2,405, 1)]

	#make vectors for pi/pi0 (or xi/xi0) values and their lower and upper CIs
	y1 = as.numeric(df_mean[1,])
	y10 = as.numeric(df_mean[10,])
	y11 = as.numeric(df_mean[11,])
	y12 = as.numeric(df_mean[12,])
	
	y1lower = as.numeric(df_lower[1,])
	y10lower = as.numeric(df_lower[10,])
	y11lower = as.numeric(df_lower[11,])
	y12lower = as.numeric(df_lower[12,])
			
	y1upper = as.numeric(df_upper[1,])
	y10upper = as.numeric(df_upper[10,])
	y11upper = as.numeric(df_upper[11,])
	y12upper = as.numeric(df_upper[12,])
	
	#make table for Norborg 1996 expectations for pi/pi0
	BGS_est_models10to12 = read.table(nord96_models10to12, stringsAsFactors=F, header=T) 
	
	#make x-axis labels and positions for those labels
	xlab = c('-1.0', '', '-0.8', '', '-0.6', '', '-0.4', '', '-0.2', '', '0.0')
	xpos = c(4, seq(44, 404, 40))
	
	#plot models 1, 2, and 3
	subplot_CI(y1, x, y10, y11, y1lower, y1upper, y10lower, y10upper, y11lower, y11upper, lowerlabel, TRUE, xlab, xpos, YTYPE=1)
	#add Norborg 1996 expectations for pi/pi0
	if(PIS == TRUE){
		lines(BGS_est_models10to12$pos, BGS_est_models10to12$model10, lty=3, lwd=LWD, col='red')
		lines(BGS_est_models10to12$pos, BGS_est_models10to12$model11, lty=3, lwd=LWD, col='blue')
	}
	
	#plot models 1 and 4
	subplot_CI(y1, x, y12, y12, y1lower, y1upper, y12lower, y12upper, y12lower, y12upper, lowerlabel, FALSE, xlab, xpos, YTYPE=2)
	#add Norborg 1996 expectations for pi/pi0
	if(PIS == TRUE){
		lines(BGS_est_models10to12$pos, BGS_est_models10to12$model12, lty=3, lwd=LWD, col='purple')
	}
	axis(2, at=yseq, tck=tck, labels=FALSE)
}


############
##plot xis##
############
plot_func(mean_xi_file, lower_xi_file, upper_xi_file, 0.64, 1.175, 0.1, 0.89, 0.97, 0.015, 3, FALSE, FALSE)

############
##plot pis##
############
plot_func(mean_pi_file, lower_pi_file, upper_pi_file, 0.45, 1, 0.1, 0.54, 0.63, 0.015, 3, TRUE, TRUE)

dev.off()
