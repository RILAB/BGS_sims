#Note that demographic model names were changed during the course of manuscript writing so that the original input file names and model names therein
#no longer correspond to the demographic model names used in the manuscript. Namely, demographic models 2-9 were re-labeled as models 5-12
#and demographic models 10-12 were re-labeled as models 2-4. This script makes the necessary changes to fit model names to those used for the
#final manuscript, but it still assumes that the input file names and model names therein correspond to the original model names.
#Additionally, arguments for this script also use the original model names.

args = commandArgs(TRUE)

print(paste('data directory: ', args[1], sep=' '))
print(paste('demographic models: ', args[2], sep=' '))

#############################
###start plot demographies###
#############################

#set working directory to directory where files for observed pi, pi0, xi, and xi0 live
setwd(args[1])

#demographic models to plot (either '25' for models 5-8 [renamed] or '69' for models 9-12 [renamed] or '1012' for models 2-4 [renamed])
models = args[2]

#get first and last model numbers from our models to plot (this will help to index dfs later on)
if(models == '1012'){
	firstmodel=10
	lastmodel=12
} else{
	firstmodel = as.numeric(strsplit(models, '')[[1]][1])
	lastmodel = as.numeric(strsplit(models, '')[[1]][2])
}

#observed pi, pi0, xi, and xi0 from simulations and bootstrapped 95 percent CIs (for genetic window closest to selected region [0.00082 cM away])
mean_bgs_pi_file = 'bgs_mean_pi_generic_0.00082cM.txt'
lower_bgs_pi_file = 'bgs_95CI_lower_pi_generic_0.00082cM.txt'
upper_bgs_pi_file = 'bgs_95CI_upper_pi_generic_0.00082cM.txt'
mean_bgs_xi_file = 'bgs_mean_xi_generic_0.00082cM.txt'
lower_bgs_xi_file = 'bgs_95CI_lower_xi_generic_0.00082cM.txt'
upper_bgs_xi_file = 'bgs_95CI_upper_xi_generic_0.00082cM.txt'
mean_neut_pi_file = 'neut_mean_pi_generic_0.00082cM.txt'
lower_neut_pi_file = 'neut_95CI_lower_pi_generic_0.00082cM.txt'
upper_neut_pi_file = 'neut_95CI_upper_pi_generic_0.00082cM.txt'
mean_neut_xi_file = 'neut_mean_xi_generic_0.00082cM.txt'
lower_neut_xi_file = 'neut_95CI_lower_xi_generic_0.00082cM.txt'
upper_neut_xi_file =  'neut_95CI_upper_xi_generic_0.00082cM.txt'

#x-axis labels (representing time) and positions for those labels
if(models == '25' | models == '1012'){
	xlab = c('-1.0', '', '-0.8', '', '-0.6', '', '-0.4', '', '-0.2', '', '0.0')
	xpos = c(2, seq(12, 102, 10))
} else {
	xlab = c('-0.1', '', '-0.08', '', '-0.06', '', '-0.04', '', '-0.02', '', '0.0')
	xpos = c(1.000, 1.000, 2.000, 3.000, 4.000, 5.000, 6.000, 7.000, 8.000, 9.000, 10.000, 11.000)
}

tck = -0.064
ymin = 100
ymax = 1000000
yseq = c(100, 1000, 10000, 100000, 1000000)

#demographic model information
demmodels1to5_and_10to12 = read.csv('new_generic_models_0.01gen_models1to5_and_models10to12.txt', header=TRUE)
demmodels6to9 = read.csv('new_generic_models_0.01gen_models6to9.txt', header=TRUE)
colnames(demmodels1to5_and_10to12) = c('y1', 'y2', 'y3', 'y4', 'y5', 'y10', 'y11', 'y12', 'gen', 'X')
colnames(demmodels6to9) = c('y1', 'y6', 'y7', 'y8', 'y9', 'gen', 'X')

if(models == '1012'){
	pdf('model2to4_0.00082cM.pdf', height=6, width=7.5)
	par(mfrow=c(3,3), oma=c(5.5,8.7,0,0) + 0.1, mar=c(0,0,0.9,0.9), las=1)
} else if(models == '25'){
	pdf('model5to8_0.00082cM.pdf', height=6, width=9.5)
	par(mfrow=c(3,4), oma=c(5.5,8.55,0,0) + 0.1, mar=c(0,0,0.9,0.9), las=1)
} else if(models == '69'){
	pdf('model9to12_0.00082cM.pdf', height=6, width=9.5)
	par(mfrow=c(3,4), oma=c(5.5,8.55,0,0) + 0.1, mar=c(0,0,0.9,0.9), las=1)
}

#set an lwd variable for line widths
w=3

#initiate empty plot to plot demographic models 
if(models == '25' | models == '1012'){
	plot(demmodels1to5_and_10to12$X, demmodels1to5_and_10to12$y1, ylim=c(ymin, ymax), type='n', xaxt='n', cex=2, log='y', yaxt='n', ylab='')
} else{
	plot(demmodels6to9$X, demmodels6to9$y1, ylim=c(ymin, ymax), type='n', xaxt='n', cex=2, log='y', yaxt='n', ylab='')
}

#make grid lines
abline(h=c(yseq), v=c(xpos), col="darkgray", lty="dotted")

#this is a plot in the first column of plots, so add y-axis title and labels
text(par("usr")[1]-(0.54*(par("usr")[2]-par("usr")[1])), 10000, 'Effective Population Size', srt=90, xpd=NA, cex=1.4)
text(par("usr")[1]-(0.44*(par("usr")[2]-par("usr")[1])), 10000, expression(paste("(", italic("N"[e]), "; log-scale)", sep='')), srt=90, xpd=NA, cex=1.4)
axis(2, at=yseq, labels=FALSE)
text(par("usr")[1]-(0.076*(par("usr")[2]-par("usr")[1])), yseq, labels=c('100', '1000', '10000', '100000', '1000000'), pos=NULL, xpd=NA, cex=1.25, adj=1)

#add bottom tick marks, these will be longer than the default so that they connect the plots within the figure
axis(1, at=xpos, tck=tck, labels=FALSE)

legendcex = 1.15

#add lines for demographic model population sizes through time and legends indicating model name
if(models == '25'){
	lines(demmodels1to5_and_10to12$X, demmodels1to5_and_10to12$y2, lty=1, lwd=w, col='black')
	legend('topleft', legend=c('model 5'), col=c('black'), lty=1, lwd=w, cex=legendcex, box.lwd=2, bg='white')
} else if(models == '1012'){
	lines(demmodels1to5_and_10to12$X, demmodels1to5_and_10to12$y10, lty=1, lwd=w, col='black')
	legend('topleft', legend=c('model 2'), col=c('black'), lty=1, lwd=w, cex=legendcex, box.lwd=2, bg='white')
} else{
	lines(demmodels6to9$X, demmodels6to9$y6, lty=1, lwd=w, col='black')
	legend('topleft', legend=c('model 9'), col=c('black'), lty=1, lwd=w, cex=legendcex, box.lwd=2, bg='white')
}

box(lwd=2)

#plot other demographic models (similar to above)
if(models == '25'){ #plots for models 5 to 8
	for(i in seq(3, 5)){
		plot(demmodels1to5_and_10to12$X, demmodels1to5_and_10to12$y1, ylim=c(ymin, ymax), type='n', xaxt='n', yaxt='n', cex=2, log='y')
		abline(h=c(yseq), v=c(xpos), col="darkgray", lty="dotted") #make grid lines
		
		#add bottom and left tick marks, these will be longer than the default so that they connect the plots within the figure
		axis(1, at=xpos, tck=tck, labels=FALSE)
		axis(2, at=yseq, tck=tck, labels=FALSE)
		
		lines(demmodels1to5_and_10to12$X, demmodels1to5_and_10to12[,i], lty=1, lwd=w, col='black') #add lines for demographic model population sizes through time
		legend('topleft', legend=c(paste('model ', i+3,  sep='')), col=c('black'), lty=1, lwd=w, cex=legendcex, box.lwd=2, bg='white') #legend indicating model name
		box(lwd=2)
	}
} else if(models == '1012'){ #plots for models 2 to 4
	for(i in seq(7, 8)){
		plot(demmodels1to5_and_10to12$X, demmodels1to5_and_10to12$y1, ylim=c(ymin, ymax), type='n', xaxt='n', yaxt='n', cex=2, log='y')
		abline(h=c(yseq), v=c(xpos), col="darkgray", lty="dotted") #make grid lines
		
		#add bottom and left tick marks, these will be longer than the default so that they connect the plots within the figure
		axis(1, at=xpos, tck=tck, labels=FALSE)
		axis(2, at=yseq, tck=tck, labels=FALSE)
		lines(demmodels1to5_and_10to12$X, demmodels1to5_and_10to12[,i], lty=1, lwd=w, col='black') #add lines for demographic model population sizes through time
		legend('topleft', legend=c(paste('model ', i-4,  sep='')), col=c('black'), lty=1, lwd=w, cex=legendcex, box.lwd=2, bg='white') #legend indicating model name
		box(lwd=2)
	}
} else { #plots for models 9 to 12
	for(i in seq(3, 5)){
		plot(demmodels6to9$X, demmodels6to9$y1, ylim=c(ymin, ymax), type='n', xaxt='n', yaxt='n', cex=2, log='y')
		abline(h=c(yseq), v=c(xpos), col="darkgray", lty="dotted") #make grid lines
		
		#add bottom and left tick marks, these will be longer than the default so that they connect the plots within the figure
		axis(1, at=xpos, tck=tck, labels=FALSE)
		axis(2, at=yseq, tck=tck, labels=FALSE)
		
		lines(demmodels6to9$X, demmodels6to9[,i], lty=1, lwd=w, col='black') #add lines for demographic model population sizes through time
		legend('topleft', legend=c(paste('model ', i+7,  sep='')), col=c('black'), lty=1, cex=legendcex, lwd=w, box.lwd=2, bg='white') #legend indicating model name
		box(lwd=2)
	}
}

###########################
###end plot demographies###
###########################


#if plotting xi and pi for models 2-8, need all generations from data file
#if plotting xi and pi for models 9-12, only need last 42 generations from data file (so start at 364)
if(models == '25' | models == '1012'){
	start = 2
} else{
	start = 364
}

#make function to plot pi, pi0, xi, and xi0 for each generation in the simulation model 
plot_func = function(meanfile_neut, CIlowerfile_neut, CIupperfile_neut, meanfile_bgs, CIlowerfile_bgs, CIupperfile_bgs, YMIN, YMAX, BY, LWD, lowerlabel=FALSE){

	if(models == '25' | models == '1012'){
		#make vector for positions along x-axis for plotting pi/pi0 and xi/xi0 for each generation (404 total for models 1-8)
		x = c(0.8, 2, 3.2, 4, 4.8, 6, 7.2, 8, 8.8, 10, 11.2, 12, 12.8, 14, 15.2, 16, 16.8, 18, 19.2, 20, 20.8, 22, 23.2, 24, 24.8, 26, 27.2, 28, 28.8, 30, 31.2, 32, 32.8, 34, 35.2, 36, 36.8, 38, 39.2, 40, 40.8, 42, 43.2, 44, 44.8, 46, 47.2, 48, 48.8, 50, 51.2, 52, 52.8, 54, 55.2, 56, 56.8, 58, 59.2, 60, 60.8, 62, 63.2, 64, 64.8, 66, 67.2, 68, 68.8, 70, 71.2, 72, 72.8, 74, 75.2, 76, 76.8, 78, 79.2, 80, 80.8, 82, 83.2, 84, 84.8, 86, 87.2, 88, 88.8, 90, 91.2, 92, 92.8, 94, 95.2, 96, 96.8, 98, 99.2, 100, 100.8, 102, 103.2, 104, 104.8, 106, 107.2, 108, 108.8, 110, 111.2, 112, 112.8, 114, 115.2, 116, 116.8, 118, 119.2, 120, 120.8, 122, 123.2, 124, 124.8, 126, 127.2, 128, 128.8, 130, 131.2, 132, 132.8, 134, 135.2, 136, 136.8, 138, 139.2, 140, 140.8, 142, 143.2, 144, 144.8, 146, 147.2, 148, 148.8, 150, 151.2, 152, 152.8, 154, 155.2, 156, 156.8, 158, 159.2, 160, 160.8, 162, 163.2, 164, 164.8, 166, 167.2, 168, 168.8, 170, 171.2, 172, 172.8, 174, 175.2, 176, 176.8, 178, 179.2, 180, 180.8, 182, 183.2, 184, 184.8, 186, 187.2, 188, 188.8, 190, 191.2, 192, 192.8, 194, 195.2, 196, 196.8, 198, 199.2, 200, 200.8, 202, 203.2, 204, 204.8, 206, 206.8, 208, 208.8, 210, 211.2, 212, 212.8, 214, 215.2, 216, 216.8, 218, 219.2, 220, 220.8, 222, 223.2, 224, 224.8, 226, 227.2, 228, 228.8, 230, 231.2, 232, 232.8, 234, 235.2, 236, 236.8, 238, 239.2, 240, 240.8, 242, 243.2, 244, 244.8, 246, 247.2, 248, 248.8, 250, 251.2, 252, 252.8, 254, 255.2, 256, 256.8, 258, 259.2, 260, 260.8, 262, 263.2, 264, 264.8, 266, 267.2, 268, 268.8, 270, 271.2, 272, 272.8, 274, 275.2, 276, 276.8, 278, 279.2, 280, 280.8, 282, 283.2, 284, 284.8, 286, 287.2, 288, 288.8, 290, 291.2, 292, 292.8, 294, 295.2, 296, 296.8, 298, 299.2, 300, 300.8, 302, 303.2, 304, 304.8, 306, 307.2, 308, 308.8, 310, 311.2, 312, 312.8, 314, 315.2, 316, 316.8, 318, 319.2, 320, 320.8, 322, 323.2, 324, 324.8, 326, 327.2, 328, 328.8, 330, 331.2, 332, 332.8, 334, 335.2, 336, 336.8, 338, 339.2, 340, 340.8, 342, 343.2, 344, 344.8, 346, 347.2, 348, 348.8, 350, 351.2, 352, 352.8, 354, 355.2, 356, 356.8, 358, 359.2, 360, 360.8, 362, 363.2, 364, 364.8, 366, 367.2, 368, 368.8, 370, 371.2, 372, 372.8, 374, 375.2, 376, 376.8, 378, 379.2, 380, 380.8, 382, 383.2, 384, 384.8, 386, 387.2, 388, 388.8, 390, 391.2, 392, 392.8, 394, 395.2, 396, 396.8, 398, 399.2, 400, 400.8, 402, 403.2, 404)
		#make vector of positions for placing x-axis labels
		xpos = seq(4, 404, 40)
		#x-axis labels
		xlab = c('-1.0', '', '-0.8', '', '-0.6', '', '-0.4', '', '-0.2', '', '0.0')
	} else{
		#make vector for positions along x-axis for plotting pi/pi0 and xi/xi0 for each generation (42 total for models 9-12)
		x = c(1.5, 2, 2.8, 4, 5.2, 6, 6.8, 8, 9.2, 10, 10.8, 12, 13.2, 14, 14.8, 16, 17.2, 18, 18.8, 20, 21.2, 22, 22.8, 24, 25.2, 26, 26.8, 28, 29.2, 30, 30.8, 32, 33.2, 34, 34.8, 36, 37.2, 38, 38.8, 40, 41.2, 42)
		#make vector of positions for placing x-axis labels
		xpos = seq(2, 42, 4)
		#x-axis labels
		xlab = c('-0.1', '', '-0.08', '', '-0.06', '', '-0.04', '', '-0.02', '', '0.0')
	}

	#make vector for y-axis labels and positions
	yseq=seq(YMIN, YMAX, BY)

	#y-axis labels in E-notation
	ylabs = as.character(formatC(yseq, format = "e", digits = 1))
	
	if(ylabs[1] == "0.0e+00"){
		ylabs[1] = '0.0'
	}
	
	subplot_CI = function(YCONST, X, YDEM1, YDEM2, YCONST_CIlower, YCONST_CIupper, YDEM1_CIlower, YDEM1_CIupper, YDEM2_CIlower, YDEM2_CIupper, XLABEL=lowerlabel, YLABEL=FALSE){
		
		#start empty plot
		if(models == '25' | models == '1012'){ #models 2-8
			plot(X, YCONST, ylim=c(YMIN, YMAX), type='n', xaxt='n', yaxt='n', cex=2)
		} else{ #models 9-12 (note shortened xlim due to lesser number of generations to plot)
			plot(X, YCONST, ylim=c(YMIN, YMAX), xlim=c(1.6, 42), type='n', xaxt='n', yaxt='n', cex=2)
		}
		
		#make grid lines
		abline(h=c(yseq), v=c(xpos), col="darkgray", lty="dotted")
		
		#if YLABEL==TRUE, then we are plotting first column of plots so we need tick marks and a y-axis title		
		if(YLABEL==TRUE){
			axis(2, at=yseq, labels=FALSE)
			text(par("usr")[1]-(0.076*(par("usr")[2]-par("usr")[1])), yseq, labels=ylabs, pos=NULL, xpd=NA, cex=1.25, adj=1)
			if(lowerlabel == FALSE){ #plotting row of singleton (xi) plots so add corresponding yaxis title
			 	text(par("usr")[1]-(0.49*(par("usr")[2]-par("usr")[1])), (YMIN+YMAX)/2, expression(paste('Singleton Density (', xi, ')', sep='')), srt=90, xpd=NA, cex=1.4)
			} else{ #plotting row of diversity (pi) plots so add corresponding yaxis title
			 	text(par("usr")[1]-(0.49*(par("usr")[2]-par("usr")[1])), (YMIN+YMAX)/2, expression(paste('Diversity (', pi, ')', sep='')), srt=90, xpd=NA, cex=1.4)
			 }
		}
		
		#if XLABEL==TRUE, then we are plotting bottom row of plots so we need tick marks and an x-axis title
		if(XLABEL == TRUE){
			axis(1, at=xpos, labels=FALSE)
			text(xpos, par("usr")[3] + 2.5*par("usr")[3], labels=xlab, srt=45, adj=1, pos=NULL, xpd=NA, cex=1.25)
			if(models == '25' | models == '1012'){ #models 2-8
				text(207, par("usr")[3] + 9.5*par("usr")[3],  expression(paste('Time (', italic('N'[anc]), ' generations)', sep='')), xpd=NA, cex=1.4)
			} else{ #models 9-12 
				text(22.25, par("usr")[3] + 9.5*par("usr")[3],  expression(paste('Time (', italic('N'[anc]), ' generations)', sep='')), xpd=NA, cex=1.4)
			}
		} else{ #XLABEL==FALSE, so we are plotting second row of plots, so only need tick marks (but need longer tick marks to connect to bottom row of plots)
			axis(1, at=xpos, tck=tck, labels=FALSE)
		}

		#add envelopes and lines representing 95% CIs and values for pi (or xi) (under BGS)
		polygon(c(X, rev(X)), c(YDEM2_CIlower, rev(YDEM2_CIupper)), col=rgb(134/255, 37/255, 203/255, 0.5), border=NA)
		lines(X, YDEM2, lty=1, lwd=LWD, col='darkviolet')
		
		#add envelopes and lines representing 95% CIs and values for pi0 (or xi0) (under neutrality)
		polygon(c(X, rev(X)), c(YDEM1_CIlower, rev(YDEM1_CIupper)), col=rgb(240/255, 133/255, 50/255, 0.5), border=NA)
		lines(X, YDEM1, lty=1, lwd=LWD, col='darkorange1')
		box(lwd=2)
	}
	
	#read data files with pi and pi0 (or xi and xi0) values
	df_mean_bgs = read.table(meanfile_bgs, header=FALSE)[,seq(start,405, 1)][c(1, seq(firstmodel, lastmodel)), ]
	df_lower_bgs = read.table(CIlowerfile_bgs, header=FALSE)[,seq(start,405, 1)][c(1, seq(firstmodel, lastmodel)), ]
	df_upper_bgs = read.table(CIupperfile_bgs, header=FALSE)[,seq(start,405, 1)][c(1, seq(firstmodel, lastmodel)), ]
	
	df_mean_neut = read.table(meanfile_neut, header=FALSE)[,seq(start,405, 1)][c(1, seq(firstmodel, lastmodel)), ]
	df_lower_neut = read.table(CIlowerfile_neut, header=FALSE)[,seq(start,405, 1)][c(1, seq(firstmodel, lastmodel)), ]
	df_upper_neut = read.table(CIupperfile_neut, header=FALSE)[,seq(start,405, 1)][c(1, seq(firstmodel, lastmodel)), ]

	#divide by 10000 (locus size) to get per-site pi (or xi) statistics
	df_mean_bgs = df_mean_bgs/10000
	df_lower_bgs = df_lower_bgs/10000
	df_upper_bgs = df_upper_bgs/10000
	df_mean_neut = df_mean_neut/10000
	df_lower_neut = df_lower_neut/10000
	df_upper_neut = df_upper_neut/10000	

	#make vectors for pi and pi0 (or xi and xi0) values and their lower and upper CIs
	y1_bgs = as.numeric(df_mean_bgs[1,])
	y2_bgs = as.numeric(df_mean_bgs[2,])
	y3_bgs = as.numeric(df_mean_bgs[3,])
	y4_bgs = as.numeric(df_mean_bgs[4,])
	
	y1lower_bgs = as.numeric(df_lower_bgs[1,])
	y2lower_bgs = as.numeric(df_lower_bgs[2,])
	y3lower_bgs = as.numeric(df_lower_bgs[3,])
	y4lower_bgs = as.numeric(df_lower_bgs[4,])
	
	y1upper_bgs = as.numeric(df_upper_bgs[1,])
	y2upper_bgs = as.numeric(df_upper_bgs[2,])
	y3upper_bgs = as.numeric(df_upper_bgs[3,])
	y4upper_bgs = as.numeric(df_upper_bgs[4,])
	
	y1_neut = as.numeric(df_mean_neut[1,])
	y2_neut = as.numeric(df_mean_neut[2,])
	y3_neut = as.numeric(df_mean_neut[3,])
	y4_neut = as.numeric(df_mean_neut[4,])
	
	y1lower_neut = as.numeric(df_lower_neut[1,])
	y2lower_neut = as.numeric(df_lower_neut[2,])
	y3lower_neut = as.numeric(df_lower_neut[3,])
	y4lower_neut = as.numeric(df_lower_neut[4,])
	
	y1upper_neut = as.numeric(df_upper_neut[1,])
	y2upper_neut = as.numeric(df_upper_neut[2,])
	y3upper_neut = as.numeric(df_upper_neut[3,])
	y4upper_neut = as.numeric(df_upper_neut[4,])
	
	if(models != '1012'){	
		y5_bgs = as.numeric(df_mean_bgs[5,])
		y5lower_bgs = as.numeric(df_lower_bgs[5,])
		y5upper_bgs = as.numeric(df_upper_bgs[5,])
		y5_neut = as.numeric(df_mean_neut[5,])
		y5lower_neut = as.numeric(df_lower_neut[5,])
		y5upper_neut = as.numeric(df_upper_neut[5,])
	}
		
	#plot models pi and pi0 (or xi and xi0)
	subplot_CI (y1_neut, x, y2_neut, y2_bgs, y1lower_neut, y1upper_neut, y2lower_neut, y2upper_neut, y2lower_bgs, y2upper_bgs, lowerlabel, TRUE)
	subplot_CI (y1_neut, x, y3_neut, y3_bgs, y1lower_neut, y1upper_neut, y3lower_neut, y3upper_neut, y3lower_bgs, y3upper_bgs, lowerlabel, FALSE)
	axis(2, at=yseq, tck=tck, labels=FALSE)
	subplot_CI (y1_neut, x, y4_neut, y4_bgs, y1lower_neut, y1upper_neut, y4lower_neut, y4upper_neut, y4lower_bgs, y4upper_bgs, lowerlabel, FALSE)
	axis(2, at=yseq, tck=tck, labels=FALSE)
	
	if(models != '1012'){ #for models 5-12, we have an additional column of plots
		subplot_CI (y1_neut, x, y5_neut, y5_bgs, y1lower_neut, y1upper_neut, y5lower_neut, y5upper_neut, y5lower_bgs, y5upper_bgs, lowerlabel, FALSE)
		axis(2, at=yseq, tck=tck, labels=FALSE)
	}
}

##################
##start plot xis##
##################

df_mean_bgs = read.table(mean_bgs_xi_file, header=FALSE)[,seq(start,405, 1)][seq(firstmodel, lastmodel), ]
df_mean_neut = read.table(mean_neut_xi_file, header=FALSE)[,seq(start,405, 1)][seq(firstmodel, lastmodel), ]

#make a matrix of maximum xi and xi0 values to get upper y-axis limit 
if(models == '1012'){
	my_min_max_mat = matrix(nrow=6, ncol=2)

	for(i in seq(1,3)){
	my_min_max_mat[i, 1] = min(df_mean_neut[i,])
	my_min_max_mat[i, 2] = max(df_mean_neut[i,])
	}

	for(i in seq(4,6)){
	my_min_max_mat[i, 1] = min(df_mean_bgs[i-3,])
	my_min_max_mat[i, 2] = max(df_mean_bgs[i-3,])
	}

} else{
	my_min_max_mat = matrix(nrow=8, ncol=2)
	
	for(i in seq(1,4)){
	my_min_max_mat[i, 1] = min(df_mean_neut[i,])
	my_min_max_mat[i, 2] = max(df_mean_neut[i,])
	}
	
	for(i in seq(5,8)){
	my_min_max_mat[i, 1] = min(df_mean_bgs[i-4,])
	my_min_max_mat[i, 2] = max(df_mean_bgs[i-4,])
	}
}

#get upper y-axis limit for plotting
ymax = max(my_min_max_mat/10000) + 0.08*max(my_min_max_mat/10000)

ymin=0
numdigits=4
by = round((ymax-ymin)/5, digits=numdigits)

plot_func(mean_neut_xi_file, lower_neut_xi_file, upper_neut_xi_file, mean_bgs_xi_file, lower_bgs_xi_file, upper_bgs_xi_file, ymin, ymax, by, w, FALSE)

################
##end plot xis##
################

##################
##start plot pis##
##################

df_mean_bgs = read.table(mean_bgs_pi_file, header=FALSE)[,seq(start,405, 1)][seq(firstmodel, lastmodel), ]
df_mean_neut = read.table(mean_neut_pi_file, header=FALSE)[,seq(start,405, 1)][seq(firstmodel, lastmodel), ]

#make a matrix of maximum pi and pi0 values to get upper y-axis limit 
if(models == '1012'){
	my_min_max_mat = matrix(nrow=6, ncol=2)

	for(i in seq(1,3)){
	my_min_max_mat[i, 1] = min(df_mean_neut[i,])
	my_min_max_mat[i, 2] = max(df_mean_neut[i,])
	}

	for(i in seq(4,6)){
	my_min_max_mat[i, 1] = min(df_mean_bgs[i-3,])
	my_min_max_mat[i, 2] = max(df_mean_bgs[i-3,])
	}

} else{
	my_min_max_mat = matrix(nrow=8, ncol=2)
	
	for(i in seq(1,4)){
	my_min_max_mat[i, 1] = min(df_mean_neut[i,])
	my_min_max_mat[i, 2] = max(df_mean_neut[i,])
	}
	
	for(i in seq(5,8)){
	my_min_max_mat[i, 1] = min(df_mean_bgs[i-4,])
	my_min_max_mat[i, 2] = max(df_mean_bgs[i-4,])
	}
}

#get upper y-axis limit for plotting
ymax = max(my_min_max_mat/10000) + 0.08*max(my_min_max_mat/10000)

by = round((ymax-ymin)/5, digits=numdigits)

plot_func(mean_neut_pi_file, lower_neut_pi_file, upper_neut_pi_file, mean_bgs_pi_file, lower_bgs_pi_file, upper_bgs_pi_file, ymin, ymax, by, w, TRUE)

################
##end plot pis##
################


#make a subplot function to plot insets for pi, pi0, xi, and xi0 with truncated x- and y-axes
subplot_func = function(meanfile_neut, CIlowerfile_neut, CIupperfile_neut, meanfile_bgs, CIlowerfile_bgs, CIupperfile_bgs, XMIN, XMAX, YMIN, YMAX, LWD, cols, LOG){
	
	if(models == '25' | models == '1012'){
		#make vector for positions along x-axis for plotting pi/pi0 and xi/xi0 for each generation (404 total for models 1-8)
		x = c(0.8, 2, 3.2, 4, 4.8, 6, 7.2, 8, 8.8, 10, 11.2, 12, 12.8, 14, 15.2, 16, 16.8, 18, 19.2, 20, 20.8, 22, 23.2, 24, 24.8, 26, 27.2, 28, 28.8, 30, 31.2, 32, 32.8, 34, 35.2, 36, 36.8, 38, 39.2, 40, 40.8, 42, 43.2, 44, 44.8, 46, 47.2, 48, 48.8, 50, 51.2, 52, 52.8, 54, 55.2, 56, 56.8, 58, 59.2, 60, 60.8, 62, 63.2, 64, 64.8, 66, 67.2, 68, 68.8, 70, 71.2, 72, 72.8, 74, 75.2, 76, 76.8, 78, 79.2, 80, 80.8, 82, 83.2, 84, 84.8, 86, 87.2, 88, 88.8, 90, 91.2, 92, 92.8, 94, 95.2, 96, 96.8, 98, 99.2, 100, 100.8, 102, 103.2, 104, 104.8, 106, 107.2, 108, 108.8, 110, 111.2, 112, 112.8, 114, 115.2, 116, 116.8, 118, 119.2, 120, 120.8, 122, 123.2, 124, 124.8, 126, 127.2, 128, 128.8, 130, 131.2, 132, 132.8, 134, 135.2, 136, 136.8, 138, 139.2, 140, 140.8, 142, 143.2, 144, 144.8, 146, 147.2, 148, 148.8, 150, 151.2, 152, 152.8, 154, 155.2, 156, 156.8, 158, 159.2, 160, 160.8, 162, 163.2, 164, 164.8, 166, 167.2, 168, 168.8, 170, 171.2, 172, 172.8, 174, 175.2, 176, 176.8, 178, 179.2, 180, 180.8, 182, 183.2, 184, 184.8, 186, 187.2, 188, 188.8, 190, 191.2, 192, 192.8, 194, 195.2, 196, 196.8, 198, 199.2, 200, 200.8, 202, 203.2, 204, 204.8, 206, 206.8, 208, 208.8, 210, 211.2, 212, 212.8, 214, 215.2, 216, 216.8, 218, 219.2, 220, 220.8, 222, 223.2, 224, 224.8, 226, 227.2, 228, 228.8, 230, 231.2, 232, 232.8, 234, 235.2, 236, 236.8, 238, 239.2, 240, 240.8, 242, 243.2, 244, 244.8, 246, 247.2, 248, 248.8, 250, 251.2, 252, 252.8, 254, 255.2, 256, 256.8, 258, 259.2, 260, 260.8, 262, 263.2, 264, 264.8, 266, 267.2, 268, 268.8, 270, 271.2, 272, 272.8, 274, 275.2, 276, 276.8, 278, 279.2, 280, 280.8, 282, 283.2, 284, 284.8, 286, 287.2, 288, 288.8, 290, 291.2, 292, 292.8, 294, 295.2, 296, 296.8, 298, 299.2, 300, 300.8, 302, 303.2, 304, 304.8, 306, 307.2, 308, 308.8, 310, 311.2, 312, 312.8, 314, 315.2, 316, 316.8, 318, 319.2, 320, 320.8, 322, 323.2, 324, 324.8, 326, 327.2, 328, 328.8, 330, 331.2, 332, 332.8, 334, 335.2, 336, 336.8, 338, 339.2, 340, 340.8, 342, 343.2, 344, 344.8, 346, 347.2, 348, 348.8, 350, 351.2, 352, 352.8, 354, 355.2, 356, 356.8, 358, 359.2, 360, 360.8, 362, 363.2, 364, 364.8, 366, 367.2, 368, 368.8, 370, 371.2, 372, 372.8, 374, 375.2, 376, 376.8, 378, 379.2, 380, 380.8, 382, 383.2, 384, 384.8, 386, 387.2, 388, 388.8, 390, 391.2, 392, 392.8, 394, 395.2, 396, 396.8, 398, 399.2, 400, 400.8, 402, 403.2, 404)
		#make vector of positions for placing x-axis labels
		xpos = seq(4, 404, 20)
		#x-axis labels
		xlab = c('-1.0', '-0.95', '-0.9', '-0.85', '-0.8', '-0.75', '-0.7', '-0.65', '-0.6', '', '', '', '-0.4', '', '', '', '-0.2', '', '', '', '0.0')
		if(XMAX > 100){
			#make vector of positions for placing x-axis labels (more sparsely distributed)
			xpos = seq(4, 404, 80)
			#x-axis labels
			xlab = c('-1.0', '-0.8', '-0.6', '-0.4', '-0.2', '0.0')
		}
	} else{
		#make vector for positions along x-axis for plotting pi/pi0 and xi/xi0 for each generation (42 total for models 9-12)
		x = c(1.5, 2, 2.8, 4, 5.2, 6, 6.8, 8, 9.2, 10, 10.8, 12, 13.2, 14, 14.8, 16, 17.2, 18, 18.8, 20, 21.2, 22, 22.8, 24, 25.2, 26, 26.8, 28, 29.2, 30, 30.8, 32, 33.2, 34, 34.8, 36, 37.2, 38, 38.8, 40, 41.2, 42)
		#make vector of positions for placing x-axis labels
		xpos = seq(2, 42, 4)
		#x-axis labels
		xlab = c('-0.1', '-0.09', '-0.08', '-0.07', '-0.06', '', '-0.04', '', '-0.02', '', '0.0')
	}

	#get distance between locations for grid points/axis tick marks along y-axis
	#if log-scaled
	LOGBY = (log(YMAX) - log(YMIN))/3
	#if linear
	BY = (YMAX-YMIN)/3
	
	#make vector of positions for grid points/axis tick marks along y-axis and make vector of corresponding labels
	if(LOG == TRUE){ #if log-scaled
		yseq = seq(log(YMIN), log(YMAX), by=LOGBY)
		ylabs = as.character(exp(yseq))
	} else{ #if linear
		ylabs = as.character(seq(YMIN, YMAX, by=BY))
		yseq=as.numeric(ylabs)
	}

	#if plot inset is for first column or third column of plots, adjust y-axis limits and rename
	if(cols == '1' | cols == '3'){
		YMINsub = YMIN+((YMAX-YMIN)*0.02525)
		YMAXsub = YMAX-((YMAX-YMIN)*0.04214)
	} else{ #don't adjust y-axis limits but still rename
		YMINsub = YMIN
		YMAXsub = YMAX
	}
	
	#relabel y-axis labels in E-notation
	ylabs = formatC(as.numeric(ylabs), format = "e", digits = 1)
	
	subplot_CI = function(YCONST, X, YDEM1, YDEM2, YCONST_CIlower, YCONST_CIupper, YDEM1_CIlower, YDEM1_CIupper, YDEM2_CIlower, YDEM2_CIupper){
		#start empty plot
		if(LOG == TRUE){ #log-scale y-limits
			plot(X, YCONST, ylim=c(log(YMINsub), log(YMAXsub)), xlim=c(XMIN, XMAX), type='n', xaxt='n', yaxt='n', cex=2)
		} else { #linear y-limits
			plot(X, YCONST, ylim=c(YMINsub, YMAXsub), xlim=c(XMIN, XMAX), type='n', xaxt='n', yaxt='n', cex=2)
		}
		
		#to make a blank canvas for our insets that won't show the background plots and grid lines
		rect(-1000, -1000, 1100, 1100, col='white')
		
		#make grid lines and axis tick marks with labels
		abline(h=c(yseq), v=c(xpos), col="darkgray", lty="dotted")
		axis(2, at=yseq, labels=ylabs, mgp=c(3, .75, 0))
		axis(1, at=xpos, labels=xlab, mgp=c(3, .5, 0))

		#add envelopes and lines representing 95% CIs and values for pi (or xi) under BGS (violet) and neutrality (orange)
		polygon(c(X, rev(X)), c(YDEM2_CIlower, rev(YDEM2_CIupper)), col=rgb(134/255, 37/255, 203/255, 0.5), border=NA)
		lines(X, YDEM2, lty=1, lwd=LWD, col='darkviolet')
		polygon(c(X, rev(X)), c(YDEM1_CIlower, rev(YDEM1_CIupper)), col=rgb(240/255, 133/255, 50/255, 0.5), border=NA)
		lines(X, YDEM1, lty=1, lwd=LWD, col='darkorange1')
		box(lwd=2)
	}
	
	#read data files with pi and pi0 (or xi and xi0) values	
	df_mean_bgs = read.table(meanfile_bgs, header=FALSE)[,seq(start,405, 1)][c(1, seq(firstmodel, lastmodel)), ]
	df_lower_bgs = read.table(CIlowerfile_bgs, header=FALSE)[,seq(start,405, 1)][c(1, seq(firstmodel, lastmodel)), ]
	df_upper_bgs = read.table(CIupperfile_bgs, header=FALSE)[,seq(start,405, 1)][c(1, seq(firstmodel, lastmodel)), ]
	
	df_mean_neut = read.table(meanfile_neut, header=FALSE)[,seq(start,405, 1)][c(1, seq(firstmodel, lastmodel)), ]
	df_lower_neut = read.table(CIlowerfile_neut, header=FALSE)[,seq(start,405, 1)][c(1, seq(firstmodel, lastmodel)), ]
	df_upper_neut = read.table(CIupperfile_neut, header=FALSE)[,seq(start,405, 1)][c(1, seq(firstmodel, lastmodel)), ]
	
	#divide by 10000 (locus size) to get per-site diversity
	df_mean_bgs = df_mean_bgs/10000
	df_lower_bgs = df_lower_bgs/10000
	df_upper_bgs = df_upper_bgs/10000
	df_mean_neut = df_mean_neut/10000
	df_lower_neut = df_lower_neut/10000
	df_upper_neut = df_upper_neut/10000	

	#log-transform dfs if LOG==TRUE
	if(LOG == TRUE){
		df_mean_bgs = log(df_mean_bgs)
		df_lower_bgs = log(df_lower_bgs)
		df_upper_bgs = log(df_upper_bgs)

		df_mean_neut = log(df_mean_neut)
		df_lower_neut = log(df_lower_neut)
		df_upper_neut = log(df_upper_neut)
	}

	#make vectors for pi and pi0 (or xi and xi0) values and their lower and upper CIs
	y1_bgs = as.numeric(df_mean_bgs[1,])
	y2_bgs = as.numeric(df_mean_bgs[2,])
	y3_bgs = as.numeric(df_mean_bgs[3,])
	y4_bgs = as.numeric(df_mean_bgs[4,])
	
	y1lower_bgs = as.numeric(df_lower_bgs[1,])
	y2lower_bgs = as.numeric(df_lower_bgs[2,])
	y3lower_bgs = as.numeric(df_lower_bgs[3,])
	y4lower_bgs = as.numeric(df_lower_bgs[4,])
	
	y1upper_bgs = as.numeric(df_upper_bgs[1,])
	y2upper_bgs = as.numeric(df_upper_bgs[2,])
	y3upper_bgs = as.numeric(df_upper_bgs[3,])
	y4upper_bgs = as.numeric(df_upper_bgs[4,])
	
	y1_neut = as.numeric(df_mean_neut[1,])
	y2_neut = as.numeric(df_mean_neut[2,])
	y3_neut = as.numeric(df_mean_neut[3,])
	y4_neut = as.numeric(df_mean_neut[4,])
	
	y1lower_neut = as.numeric(df_lower_neut[1,])
	y2lower_neut = as.numeric(df_lower_neut[2,])
	y3lower_neut = as.numeric(df_lower_neut[3,])
	y4lower_neut = as.numeric(df_lower_neut[4,])
	
	y1upper_neut = as.numeric(df_upper_neut[1,])
	y2upper_neut = as.numeric(df_upper_neut[2,])
	y3upper_neut = as.numeric(df_upper_neut[3,])
	y4upper_neut = as.numeric(df_upper_neut[4,])
	
	if(models != '1012'){	
		y5_bgs = as.numeric(df_mean_bgs[5,])
		y5lower_bgs = as.numeric(df_lower_bgs[5,])
		y5upper_bgs = as.numeric(df_upper_bgs[5,])
		y5_neut = as.numeric(df_mean_neut[5,])
		y5lower_neut = as.numeric(df_lower_neut[5,])
		y5upper_neut = as.numeric(df_upper_neut[5,])
	}
	
	#make inset plots for models' pi and pi0 (or xi and xi0) if they are in the plot column where cols == 'column number'
	if(cols=='1'){
		subplot_CI(y1_neut, x, y2_neut, y2_bgs, y1lower_neut, y1upper_neut, y2lower_neut, y2upper_neut, y2lower_bgs, y2upper_bgs)
	} else if(cols=='2'){
		subplot_CI(y1_neut, x, y3_neut, y3_bgs, y1lower_neut, y1upper_neut, y3lower_neut, y3upper_neut, y3lower_bgs, y3upper_bgs)
	} else if(cols=='3'){
		subplot_CI(y1_neut, x, y4_neut, y4_bgs, y1lower_neut, y1upper_neut, y4lower_neut, y4upper_neut, y4lower_bgs, y4upper_bgs)
	} else if(cols=='4'){
		subplot_CI(y1_neut, x, y5_neut, y5_bgs, y1lower_neut, y1upper_neut, y5lower_neut, y5upper_neut, y5lower_bgs, y5upper_bgs)
	}
}

#set x-limits for inset plots
if(models == '25' | models == '1012'){
	xmin = 5
	xmax = 50
} else {
	xmin = 2
	xmax = 20
}

#update variable for lwd (linewidth)
w=2

#set y-limits for inset plots
if(models == '69'){
	ymin1 = 0.000245
	ymax1 = 0.0006
	ymin2 = 6e-05
	ymax2 = 0.000275
	ymin3 = 0.00013
	ymax3 = 0.0003
	ymin4 = 2e-05
	ymax4 = 0.0001
} else if(models == '1012'){
	ymin1 = 0.00012
	ymax1 = 0.0002
	ymin2 = 0.0000175
	ymax2 = 0.00005
	ymin3 = 0.0011
	ymax3 = 0.0028
} else{
	ymin1 = 0.00015
	ymax1 = 0.00035
	ymin2 = 0.000025
	ymax2 = 0.00009
	ymin3 = 0.00013
	ymax3 = 0.0003
	ymin4 = 2e-05
	ymax4 = 0.0001
}


###########################
##start plot inset plots###
###########################

if(models != '1012'){
	if(models == '25'){
		#inset plot for row 2, column 1 (xi for generic model 5)
		par(fig=c(0.06+.01, 0.215-.01, 0.49, 0.652), new=T) #coordinates for inset
		subplot_func(mean_neut_xi_file, lower_neut_xi_file, upper_neut_xi_file, mean_bgs_xi_file, lower_bgs_xi_file, upper_bgs_xi_file, xmin, xmax, ymin1, ymax1, w, '1', TRUE)
	} else{
		#inset plot for row 2, column 1 (xi for generic model 9)
		par(fig=c(0.06+.01, 0.215-.03, 0.49, 0.652), new=T) #coordinates for inset
		subplot_func(mean_neut_xi_file, lower_neut_xi_file, upper_neut_xi_file, mean_bgs_xi_file, lower_bgs_xi_file, upper_bgs_xi_file, xmin, 16.7, ymin1, ymax1, w, '1', TRUE)
	}
	#inset plot for row 2, column 2 (xi for generic model 6 and 10)
	par(fig=c(0.31+.01, 0.465-.01, 0.49, 0.652), new=T) #coordinates for inset
	subplot_func(mean_neut_xi_file, lower_neut_xi_file, upper_neut_xi_file, mean_bgs_xi_file, lower_bgs_xi_file, upper_bgs_xi_file, xmin, xmax, ymin2, ymax2, w, '2', TRUE)
	
	#inset plot for row 2, column 3 (xi for generic model 7 and 11)
	par(fig=c(0.56+.01, 0.715-.01, 0.49, 0.652), new=T) #coordinates for inset
	subplot_func(mean_neut_xi_file, lower_neut_xi_file, upper_neut_xi_file, mean_bgs_xi_file, lower_bgs_xi_file, upper_bgs_xi_file, xmin, xmax, ymin3, ymax3, w, '3', TRUE)
	
	#inset plot for row 2, column 4 (xi for generic model 8 and 12)
	par(fig=c(0.81+.01, 0.965-.01, 0.49, 0.652), new=T) #coordinates for inset
	subplot_func(mean_neut_xi_file, lower_neut_xi_file, upper_neut_xi_file, mean_bgs_xi_file, lower_bgs_xi_file, upper_bgs_xi_file, xmin, xmax, ymin4, ymax4, w, '4', TRUE)
} else {
	#inset plot for row 2, column 1 (xi and xi0 for model 2)
	par(fig=c(0.117, 0.292, 0.49, 0.652), new=T) #coordinates for inset
	subplot_func(mean_neut_xi_file, lower_neut_xi_file, upper_neut_xi_file, mean_bgs_xi_file, lower_bgs_xi_file, upper_bgs_xi_file, xmin, xmax, ymin1, ymax1, w, '1', TRUE)
	
	#inset plot for row 2, column 2 (xi and xi0 for model 3)
	par(fig=c(0.45, 0.625, 0.49, 0.652), new=T) #coordinates for inset
	subplot_func(mean_neut_xi_file, lower_neut_xi_file, upper_neut_xi_file, mean_bgs_xi_file, lower_bgs_xi_file, upper_bgs_xi_file, xmin, xmax, ymin2, ymax2, w, '2', TRUE)
	
	#inset plot for row 2, column 3 (xi and xi0 for model 4)
	par(fig=c(0.784, 0.959, 0.49-0.0875, 0.652-0.0875), new=T) #coordinates for inset
	subplot_func(mean_neut_xi_file, lower_neut_xi_file, upper_neut_xi_file, mean_bgs_xi_file, lower_bgs_xi_file, upper_bgs_xi_file, xmin-xmin*0.5, xmax, ymin3, ymax3, w, '3', TRUE)
	
	#inset plot for row 3, column 2 (pi and pi0 for model 3)
	par(fig=c(0.45, 0.625, 0.49-.3325, 0.652-.3325), new=T) #coordinates for inset
	subplot_func(mean_neut_pi_file, lower_neut_pi_file, upper_neut_pi_file, mean_bgs_pi_file, lower_bgs_pi_file, upper_bgs_pi_file, xmin, 400, 0.0000215, 0.00003, w, '2', TRUE)
}

#########################
##end plot inset plots###
#########################

dev.off()
