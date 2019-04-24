library(plotrix)
library(viridis)
library(RColorBrewer)

#Note that demographic model names were changed during the course of manuscript writing so that the original input file names and model names therein
#no longer correspond to the demographic model names used in the manuscript. Namely, demographic models 2-9 were re-labeled as models 5-12
#and demographic models 10-12 were re-labeled as models 2-4. This script makes the necessary changes to fit model names to those used for the
#final manuscript, but it still assumes that the input file names and model names therein correspond to the original model names.
#Additionally, arguments for this script also use the original model names.

args = commandArgs(TRUE)

print(paste('data directory: ', args[1], sep=' '))
print(paste('diversity data file: ', args[2], sep=' '))
print(paste('singelton data file: ', args[3], sep=' '))
print(paste('demographic model: ', args[4], sep=' '))
print(paste('y-axis lower limit for pi/pi0: ', args[5], sep=' '))
print(paste('y-axis upper limit for pi/pi0: ', args[6], sep=' '))
print(paste('y-axis grid spacing for pi/pi0: ', args[7], sep=' '))
print(paste('y-axis lower limit for xi/xi0: ', args[8], sep=' '))
print(paste('y-axis upper limit for xi/xi0: ', args[9], sep=' '))
print(paste('y-axis grid spacing for xi/xi0: ', args[10], sep=' '))
print(paste('total number of generations in simulation model: ', args[11], sep=' '))

#set working directory to directory where data files live
setwd(args[1])

#total number of generations in the simulation model
total=as.numeric(args[11])

#set axis tick mark size
tck=-0.04

#set colors for each generation in the simulation model, set colors for demographic model, set colors for legend
#(all colors in simulation results, demographic model, and legend will correspond to each other)
mycolors = rev(viridis_pal(option = "D")(total-3)[seq(1, total-3, 1)])
mycolors_dem = rev(viridis_pal(option = "D")(102)[seq(1, 102, 1)])
legend_mycolors = rev(viridis_pal(option = "D")(1500)[seq(1, 1500, 1)])

#make function to plot each generation in the simulation model across all genetic distance windows (20 windows total)
plot_func = function(RECOMBFILE_PI, RECOMBFILE_XI, YSEQ_PI, YSEQ_XI, YMIN_PI, YMAX_PI, YMIN_XI, YMAX_XI, XPOS, XLAB, XLABPOS, XTITLEPOS, TOTAL, PI_BREAK, XI_BREAK1, XI_BREAK2, MODELTYPE){
	#set different sizes for text
	empty_legendcex=1.15
	gradient_legendcex=1.15
	plotlabel_cex=1.2
	xlab_ylab_values_cex=1.25
	axis_title_cex=1.5
	if(modeltype == '12'){
		empty_legendcex=empty_legendcex*0.8
		gradient_legendcex=gradient_legendcex*0.8
		plotlabel_cex=plotlabel_cex*0.8
		xlab_ylab_values_cex=xlab_ylab_values_cex*0.8
		axis_title_cex=axis_title_cex*0.8
	}
	
	#set y-axis grid line spacing
	YSEQ_PI_minor = seq(YSEQ_PI[1]+((YSEQ_PI[2]-YSEQ_PI[1])/2), YSEQ_PI[length(YSEQ_PI)]+((YSEQ_PI[2]-YSEQ_PI[1])/2), YSEQ_PI[2]-YSEQ_PI[1])
	YSEQ_XI_minor = seq(YSEQ_XI[1]+((YSEQ_XI[2]-YSEQ_XI[1])/2), YSEQ_XI[length(YSEQ_XI)]+((YSEQ_XI[2]-YSEQ_XI[1])/2), YSEQ_XI[2]-YSEQ_XI[1])
	
	#initiate blank plot for first plot of pi/pi0 across all genetic distance windows
	plot(seq(1, 20, 1), RECOMBFILE_PI[1,], ylim=c(YMIN_PI, YMAX_PI), type='n', xaxt='n', yaxt='n')
	axis(2, at=YSEQ_PI, tck=-0.025, labels=as.character(YSEQ_PI), cex.axis=xlab_ylab_values_cex)
	axis(2, at=YSEQ_PI_minor, tck=-0.025, labels=FALSE)
	axis(1, at=XPOS, tck=tck, labels=FALSE)
	if(!(modeltype %in% c('6', '7', '8', '9', '12'))){
		axis(4, at=YSEQ_PI, tck=tck, labels=FALSE)
		axis(4, at=YSEQ_PI_minor, tck=tck, labels=FALSE)
	}
	abline(h=c(YSEQ_PI, YSEQ_PI_minor), v=c(XPOS), col='darkgray', lty="dotted") #make plot grids
	
	#below is for making an empty legend box, then placing in a plot label (i.e., 'A', 'B', 'C', 'D', etc.) that corresponds to the gradient legend and shows which generations are being plotted
	#get legend box coordinates
	my.legend.size = legend('bottomright', legend=expression(italic('plot')), col='', pch=NA, cex=empty_legendcex, box.lwd=2, bg='white') #note label 'plot' is just filler, it will not be printed
	#make an empty legend box
	legend('bottomright', legend=expression(italic('plot')), text.col='white', col='white', pch=NA, cex=empty_legendcex, box.lwd=2, bg='white') #note label 'plot' is just filler, it will not be printed
	#fill in legend box with a plot label
	text(my.legend.size$text$x+0.08, my.legend.size$text$y, 'A', cex=plotlabel_cex)
	
	#add y-axis title
	text(-7.7, (YMIN_PI+YMAX_PI)/2, "Relative Diversity", srt=90, xpd=NA, cex=axis_title_cex)
	text(-6, (YMIN_PI+YMAX_PI)/2, expression(paste('(', pi, '/', pi[0], ')', sep='')), srt=90, xpd=NA, cex=axis_title_cex)
	box(lwd=2)
	
	#now go through each generation and plot pi/pi0 across all genetic distance windows using lines() function
	#start with generation 1
	lines(seq(1, 20, 1), RECOMBFILE_PI[1,], lwd=2, col=mycolors[1]) #note color of each line will change each generation so that it matches the gradient legend 
	
	#we will use this to keep track of any change in direction in pi/pi0 (see below)
	minimumx=min(RECOMBFILE_PI[,1])
	
	if(modeltype %in% c('6', '7', '8', '9', '12')){
		#these models will only have a single plot!!
		
		#go through each generation from 2 to TOTAL-3
		#TOTAL = number of generations; we use TOTAL-3 because the first 3 generations occur before demography begins (i.e., they are part of the burn-in)
		for(j in 2:(TOTAL-3)){
			lines(seq(1, 20, 1), RECOMBFILE_PI[j,], lwd=2, col=mycolors[j]) #note color of each line will change each generation so that it matches the gradient legend
		}
		#plotted all generations for these models
		#now redraw in dashed gray and red lines the first and last generation to make it easier to ID where the first generation begins and where the last generation ends on the plot
		lines(seq(1, 20, 1), RECOMBFILE_PI[1,], lwd=2, col='red', lty=2) #first generation
		lines(seq(1, 20, 1), RECOMBFILE_PI[(TOTAL-3),], lwd=3, col='red', lty=3) #last generation
		
		#empty plot to make room for gradient legend
		plot(1:10, 1:10, pch='', bty='l', type='n', xaxt='n', yaxt='n', axes = FALSE)
		#additional empty plot to make room for gradient legend
		if(modeltype != '12'){
			plot(1:10, 1:10, pch='', bty='l', type='n', xaxt='n', yaxt='n', axes = FALSE)
		}
		
	} else{
		#these models will have more than one plot!!
		
		#go through each generation
		#TOTAL = number of generations; we use TOTAL-3 below because the first 3 generations occur before demography begins (i.e., they are part of the burn-in)
		for(j in 2:(TOTAL-3)){
			#make a variable, 'newstart', to account for change in direction of pi/pi0
			#if a change in direction is encountered, make a new plot (this will prevent overlap of generations on the y-axis)
			newstart=RECOMBFILE_PI[j,1]
			
			if(newstart == minimumx){ #we have reached a nadir in pi/pi0, make a new plot!
				FIRSTJ = j #this is the new generation to start from for plotting our new plot
				
				#redraw in dashed gray and red lines the first and last generation to make it easier to ID where the first generation begins and where the last generation ends on the previous plot
				lines(seq(1, 20, 1), RECOMBFILE_PI[1,], lwd=2, col='red', lty=2) #first generation
				lines(seq(1, 20, 1), RECOMBFILE_PI[j-1,], lwd=3, col='red', lty=3) #last generation
				
				#start making our next plot, beginning where we left off!
				plot(seq(1, 20, 1), RECOMBFILE_PI[j,], ylim=c(YMIN_PI, YMAX_PI), type='n', xaxt='n', yaxt='n')
				axis(2, at=YSEQ_PI, tck=tck, labels=FALSE)
				axis(2, at=YSEQ_PI_minor, tck=tck, labels=FALSE)
				axis(1, at=XPOS, tck=tck, labels=FALSE)
				abline(h=c(YSEQ_PI, YSEQ_PI_minor), v=c(XPOS), col='darkgray', lty="dotted") #make plot grids
				
				#below is for making an empty legend box, then placing in a plot label (i.e., 'A', 'B', 'C', 'D', etc.) that corresponds to the gradient legend and shows which generations are being plotted
				#get legend box coordinates
				my.legend.size = legend('bottomright', legend=expression(italic('plot')), col='', pch=NA, cex=empty_legendcex, box.lwd=2, bg='white') #note label 'plot' is just filler, it will not be printed
				#make an empty legend box
				legend('bottomright', legend=expression(italic('plot')), text.col='white', col='white', pch=NA, cex=empty_legendcex, box.lwd=2, bg='white') #note label 'plot' is just filler, it will not be printed
				#fill in legend box with a plot label
				text(my.legend.size$text$x+0.08, my.legend.size$text$y, 'B', cex=plotlabel_cex)
				box(lwd=2)
				
				#start plotting lines for the new plot
				lines(seq(1, 20, 1), RECOMBFILE_PI[j,], lwd=2, col=mycolors[j])
				lines(seq(1, 20, 1), RECOMBFILE_PI[j,], lwd=2, col='red', lty=2)
			} else{
				#plot each generation's pi/pi0 across the 20 windows
				lines(seq(1, 20, 1), RECOMBFILE_PI[j,], lwd=2, col=mycolors[j])
			}
		}
		
		#we have plottoed all generations
		#we now now need to redraw in dashed gray and red lines the first and last generation to make it easier to ID where the first generation begins and where the last generation ends on the last plot
		lines(seq(1, 20, 1), RECOMBFILE_PI[FIRSTJ,], lwd=2, col='red', lty=2) #first generation
		lines(seq(1, 20, 1), RECOMBFILE_PI[(TOTAL-3),], lwd=3, col='red', lty=3) #last generation
		
		#make an empty plot for placing the gradient legend
		plot(1:10, 1:10, pch='', bty='l', type='n', xaxt='n', yaxt='n', axes = FALSE)
	}
	
	##START: make gradient legend (very complicated!)
	if(total > 100){
		START = '-1.0'
		N2 = '-0.5'
	} else{
		START = '-0.1'
		N2 = '-0.05'
	}
	
	END = '0.0'
	NANC = expression(italic('N'[anc]))
	
	#reference points for placing gradient legedn
	ENDBOUND1 = 3.5
	ENDBOUND2 = 8
	BOT = 0.9
	TOP = 8.4
	MID = (BOT+TOP)/2
	
	text(5.5, 10.1, expression(paste('Time (', italic('N'[anc]), ' generations)', sep='')), cex=axis_title_cex, xpd=NA)
	text(ENDBOUND1-0.25, 9.1, expression(paste('(', pi, '/', pi[0], ')', sep='')), cex=axis_title_cex)
	text(ENDBOUND2-0.25, 9.1, expression(paste('(', xi, '/', xi[0], ')', sep='')), cex=axis_title_cex)
	if(MODELTYPE %in% c('6', '7', '8', '9', '12')){
		ADD=1.1
	} else{
		ADD=0.9
	}
	
	#place labels for generation time
	if(MODELTYPE %in% c('4', '5')){
		text(ENDBOUND1+1.06, (0.05*(TOP-BOT))+BOT+0.0925, '-0.95', cex=plotlabel_cex)
		text(ENDBOUND2+1.06, (0.05*(TOP-BOT))+BOT+0.0925, '-0.95', cex=plotlabel_cex)
		text(ENDBOUND1+0.9, BOT+0.01, START, cex=plotlabel_cex)
		text(ENDBOUND2+0.9, BOT+0.01, START, cex=plotlabel_cex)
	} else {
		text(ENDBOUND1+0.9, BOT+0.05, START, cex=plotlabel_cex)
		text(ENDBOUND2+0.9, BOT+0.05, START, cex=plotlabel_cex)
	}
	text(ENDBOUND1+ADD, MID+0.05, N2, cex=plotlabel_cex)
	text(ENDBOUND1+0.8, TOP+0.01, END, cex=plotlabel_cex)
	text(ENDBOUND2+ADD, MID+0.05, N2, cex=plotlabel_cex)
	text(ENDBOUND2+0.8, TOP+0.01, END, cex=plotlabel_cex)
	
	#use color.legend function to make actual color gradient
	color.legend(ENDBOUND1-0.5, BOT, ENDBOUND1, TOP, c(''), legend_mycolors, align="rb", gradient="y")
	color.legend(ENDBOUND2-0.5, BOT, ENDBOUND2, TOP, c(''), legend_mycolors, align="rb", gradient="y")
	
	#draw segments and make labels showing where each plot for pi/pi0 is being referenced on the color gradient (every color represents a distinct generation)
	#note variable PI_BREAK tells us where one plot begins and the other ends on the color gradient (i.e., where a nadir exists in pi/pi0)
	segments(ENDBOUND1, BOT, ENDBOUND1+0.25, BOT)
	segments(ENDBOUND1, MID, ENDBOUND1+0.25, MID)
	segments(ENDBOUND1, TOP, ENDBOUND1+0.25, TOP)
	segments(ENDBOUND2, BOT, ENDBOUND2+0.25, BOT)
	segments(ENDBOUND2, MID, ENDBOUND2+0.25, MID)
	segments(ENDBOUND2, TOP, ENDBOUND2+0.25, TOP)
	if(MODELTYPE %in% c('4', '5')){
		segments(ENDBOUND1, (0.05*(TOP-BOT))+BOT, ENDBOUND1+0.25, (0.05*(TOP-BOT))+BOT)
		segments(ENDBOUND2, (0.05*(TOP-BOT))+BOT, ENDBOUND2+0.25, (0.05*(TOP-BOT))+BOT)
	}
	
	SEGSTART_PI = ENDBOUND1-0.95
	
	if(MODELTYPE %in% c('6', '7', '8', '9', '12')){
		CENTERPOINT1_PI = (BOT+TOP)/2
		segments(SEGSTART_PI,BOT,SEGSTART_PI,TOP)
	} else{
		CENTERPOINT1_PI = (BOT+((PI_BREAK/(total-3))*7.5)+BOT-0.05)/2
		segments(SEGSTART_PI,BOT,SEGSTART_PI,((PI_BREAK/(total-3))*7.5)+BOT-0.05)
		segments(SEGSTART_PI,((PI_BREAK/(total-3))*7.5)+BOT+0.05,SEGSTART_PI,TOP)
		CENTERPOINT2_PI = (((PI_BREAK/(total-3))*7.5)+BOT+0.05+TOP)/2
		segments(SEGSTART_PI, ((PI_BREAK/(total-3))*7.5)+BOT-0.05, SEGSTART_PI+0.25, ((PI_BREAK/(total-3))*7.5)+BOT-0.05)
		segments(SEGSTART_PI, ((PI_BREAK/(total-3))*7.5)+BOT+0.05, SEGSTART_PI+0.25, ((PI_BREAK/(total-3))*7.5)+BOT+0.05)
		text(SEGSTART_PI-0.4, CENTERPOINT2_PI, 'B', cex=gradient_legendcex)
	}
	
	segments(SEGSTART_PI, BOT, SEGSTART_PI+0.25, BOT)
	segments(SEGSTART_PI, TOP, SEGSTART_PI+0.25, TOP)
	text(SEGSTART_PI-0.4, CENTERPOINT1_PI, 'A', cex=gradient_legendcex)
	
	#draw segments and make labels showing where each plot for xi/xi0 is being referenced on the color gradient (every color represents a distinct generation)
	#note variables XI_BREAK1 and XI_BREAK2 tell us where one plot begins and the other ends on the color gradient (i.e., where a nadir or peak exists in xi/xi0)
	SEGSTART_XI = ENDBOUND2-0.95
	
	segments(SEGSTART_XI,BOT,SEGSTART_XI+0.25,BOT)
		if(((XI_BREAK1/(total-3))*7.5)+BOT-0.05 < 1){
			segments(SEGSTART_XI,BOT,SEGSTART_XI,BOT+0.1)
			segments(SEGSTART_XI,BOT+0.1,SEGSTART_XI+0.25,BOT+0.1)
			segments(SEGSTART_XI,BOT+0.1+0.1,SEGSTART_XI+0.25,BOT+0.1+0.1)
			segments(SEGSTART_XI,((XI_BREAK2/(total-3))*7.5)+BOT-0.05,SEGSTART_XI,BOT+0.1+0.1)
			CENTERPOINT1_XI = (BOT+(BOT+0.05))/2
			CENTERPOINT2_XI = (BOT+0.05+0.05+((XI_BREAK2/(total-3))*7.5)+BOT-0.05)/2
		} else{
			segments(SEGSTART_XI,BOT,SEGSTART_XI,((XI_BREAK1/(total-3))*7.5)+BOT-0.05)
			segments(SEGSTART_XI,((XI_BREAK1/(total-3))*7.5)+BOT-0.05,SEGSTART_XI+0.25,((XI_BREAK1/(total-3))*7.5)+BOT-0.05)
			segments(SEGSTART_XI,((XI_BREAK1/(total-3))*7.5)+BOT+0.05,SEGSTART_XI+0.25,((XI_BREAK1/(total-3))*7.5)+BOT+0.05)
			segments(SEGSTART_XI,((XI_BREAK1/(total-3))*7.5)+BOT+0.05,SEGSTART_XI,((XI_BREAK2/(total-3))*7.5)+BOT-0.05)
			CENTERPOINT1_XI = (BOT+((XI_BREAK1/(total-3))*7.5)+BOT-0.05)/2
			CENTERPOINT2_XI = (((XI_BREAK1/(total-3))*7.5)+BOT+0.05 + ((XI_BREAK2/(total-3))*7.5)+BOT-0.05)/2
		}
		if(MODELTYPE != '12'){
			segments(SEGSTART_XI,((XI_BREAK2/(total-3))*7.5)+BOT-0.05,SEGSTART_XI+0.25,((XI_BREAK2/(total-3))*7.5)+BOT-0.05)
			segments(SEGSTART_XI,((XI_BREAK2/(total-3))*7.5)+BOT+0.05,SEGSTART_XI,TOP)
			segments(SEGSTART_XI,((XI_BREAK2/(total-3))*7.5)+BOT+0.05,SEGSTART_XI+0.25,((XI_BREAK2/(total-3))*7.5)+BOT+0.05)
			CENTERPOINT3_XI = (((XI_BREAK2/(total-3))*7.5)+BOT+0.05 + TOP)/2
			if(MODELTYPE %in% c('6', '7', '8', '9')){
				text(SEGSTART_XI-0.4, CENTERPOINT2_XI, 'C', cex=gradient_legendcex)
				text(SEGSTART_XI-0.4, CENTERPOINT3_XI, 'D', cex=gradient_legendcex)	
			} else{
				text(SEGSTART_XI-0.4, CENTERPOINT2_XI, 'D', cex=gradient_legendcex)
				text(SEGSTART_XI-0.4, CENTERPOINT3_XI, 'E', cex=gradient_legendcex)
			}
		} else{
			text(SEGSTART_XI-0.4, CENTERPOINT2_XI, 'C', cex=gradient_legendcex)
		}
	segments(SEGSTART_XI, TOP, SEGSTART_XI+0.25, TOP)
	if(MODELTYPE %in% c('6', '7', '8', '9', '12')){
		text(SEGSTART_XI-0.4, CENTERPOINT1_XI, 'B', cex=gradient_legendcex)		
	} else{
		text(SEGSTART_XI-0.4, CENTERPOINT1_XI, 'C', cex=gradient_legendcex)
	}
	##END: make gradient legend
	
	
	#initiate blank plot for first plot of xi/xi0 across all genetic distance windows
	plot(seq(1, 20, 1), RECOMBFILE_XI[1,], ylim=c(YMIN_XI, YMAX_XI), type='n', xaxt='n', yaxt='n')
	axis(2, at=YSEQ_XI, tck=-0.025, labels=as.character(YSEQ_XI), cex.axis=xlab_ylab_values_cex)
	axis(2, at=YSEQ_XI_minor, tck=-0.025, labels=FALSE)
	axis(1, at=XPOS, tck=-0.025, labels=FALSE)
	text(XPOS, XLABPOS, labels=XLAB, srt=45, adj=1, pos=NULL, xpd=NA, cex=xlab_ylab_values_cex)
	axis(3, at=XPOS, tck=tck, labels=FALSE)
	axis(4, at=YSEQ_XI, tck=tck, labels=FALSE)
	axis(4, at=YSEQ_XI_minor, tck=tck, labels=FALSE)
	abline(h=c(YSEQ_XI, YSEQ_XI_minor), v=c(XPOS), col='darkgray', lty="dotted") #make plot grids

	#below is for making an empty legend box, then placing in a plot label (i.e., 'A', 'B', 'C', 'D', etc.) that corresponds to the gradient legend and shows which generations are being plotted
	#get legend box coordinates
	my.legend.size = legend('bottomright', legend=expression(italic('plot')), col='', pch=NA, cex=empty_legendcex, box.lwd=2, bg='white') #note label 'plot' is just filler, it will not be printed
	#make an empty legend box
	legend('bottomright', legend=expression(italic('plot')), text.col='white', col='white', pch=NA, cex=empty_legendcex, box.lwd=2, bg='white') #note label 'plot' is just filler, it will not be printed
	if(MODELTYPE %in% c('6', '7', '8', '9', '12')){
		#fill in legend box with a plot label
		text(my.legend.size$text$x+0.08, my.legend.size$text$y, 'B', cex=plotlabel_cex)
	} else{	
		#fill in legend box with a plot label
		text(my.legend.size$text$x+0.08, my.legend.size$text$y, 'C', cex=plotlabel_cex)
	}
	
	#add y-axis title
	text(-7.7, (YMIN_XI+YMAX_XI)/2, "Relative Singleton Density", srt=90, xpd=NA, cex=axis_title_cex)
	text(-6, (YMIN_XI+YMAX_XI)/2, expression(paste('(', xi, '/', xi[0], ')', sep='')), srt=90, xpd=NA, cex=axis_title_cex)
	
	#add x-axis title	
	text(10.5, XTITLEPOS,  "Genetic Distance (cM)", xpd=NA, cex=axis_title_cex)
	box(lwd=2)
	
	#now go through each generation and plot xi/xi0 across all genetic distance windows using lines() function
	#start with generation 1
	lines(seq(1, 20, 1), RECOMBFILE_XI[1,], lwd=2, col=mycolors[1]) #note color of each line will change each generation so that it matches the gradient legend
	
	#we will use these to keep track of any change in direction in xi/xi0 (see below)
	minimumx=min(RECOMBFILE_XI[,1])
	maximumx=max(RECOMBFILE_XI[,1])
	
	#go through each generation
	#TOTAL = number of generations; we use TOTAL-3 below because the first 3 generations occur before demography begins (i.e., they are part of the burn-in)
	for(j in 2:(TOTAL-3)){
		if(MODELTYPE != '12'){
			#there will be three plots for these models!
			#make a variable, 'newstart', to account for change in direction of xi/xi0
			#if a change in direction is encountered, make a new plot (this will prevent overlap of generations on the y-axis)
			newstart=RECOMBFILE_XI[j,1]
			if(newstart == minimumx){ #we have reached a nadir in xi/xi0, make a new plot!
				FIRSTJ = j #this is the new generation to start from for plotting our new plot
				
				#redraw in dashed gray and red lines the first and last generation to make it easier to ID where the first generation begins and where the last generation ends on the previous plot
				if(modeltype != '7'){
					lines(seq(1, 20, 1), RECOMBFILE_XI[1,], lwd=2, col='red', lty=2) #first generation
					lines(seq(1, 20, 1), RECOMBFILE_XI[j-1,], lwd=3, col='red', lty=3) #last generation
				}
				
				#start making our next plot, beginning where we left off!
				plot(seq(1, 20, 1), RECOMBFILE_XI[j,], ylim=c(YMIN_XI, YMAX_XI), type='n', xaxt='n', yaxt='n')
				abline(h=c(YSEQ_XI, YSEQ_XI_minor), v=c(XPOS), col='darkgray', lty="dotted") #make plot grids
				
				#below is for making an empty legend box, then placing in a plot label (i.e., 'A', 'B', 'C', 'D', etc.) that corresponds to the gradient legend and shows which generations are being plotted
				#get legend box coordinates
				my.legend.size = legend('bottomright', legend=expression(italic('plot')), col='', pch=NA, cex=empty_legendcex, box.lwd=2, bg='white') #note label 'plot' is just filler, it will not be printed
				#make an empty legend box
				legend('bottomright', legend=expression(italic('plot')), text.col='white', col='white', pch=NA, cex=empty_legendcex, box.lwd=2, bg='white') #note label 'plot' is just filler, it will not be printed
				if(MODELTYPE %in% c('6', '7', '8', '9')){
					#fill in legend box with a plot label
					text(my.legend.size$text$x+0.08, my.legend.size$text$y, 'C', cex=plotlabel_cex)
				}else {
					#fill in legend box with a plot label
					text(my.legend.size$text$x+0.08, my.legend.size$text$y, 'D', cex=plotlabel_cex)				
					}
					
				#add x-axis title				
				text(10.5, XTITLEPOS,  "Genetic Distance (cM)", xpd=NA, cex=axis_title_cex)
				box(lwd=2)
				
				axis(2, at=c(YSEQ_XI, YSEQ_XI_minor), tck=tck, labels=FALSE)
				if(!(modeltype %in% c('6', '7', '8', '9'))){
					axis(3, at=XPOS, tck=tck, labels=FALSE)
				}
				axis(4, at=c(YSEQ_XI, YSEQ_XI_minor), tck=tck, labels=FALSE)
				axis(1, at=XPOS, tck=-0.025, labels=FALSE)
				text(XPOS, XLABPOS, labels=XLAB, srt=45, adj=1, pos=NULL, xpd=NA, cex=xlab_ylab_values_cex)
				
				#plot xi/xi0 across the 20 windows for new plot
				lines(seq(1, 20, 1), RECOMBFILE_XI[j,], lwd=2, col=mycolors[j])

			} else if(newstart == maximumx){ #we have reached a peak in xi/xi0, make a new plot!
				#redraw in dashed gray and red lines the first and last generation to make it easier to ID where the first generation begins and where the last generation ends on the previous plot
				lines(seq(1, 20, 1), RECOMBFILE_XI[FIRSTJ,], lwd=2, col='red', lty=2) #first generation
				lines(seq(1, 20, 1), RECOMBFILE_XI[j-1,], lwd=3, col='red', lty=3) #last generation
				
				#start making our next plot, beginning where we left off!
				plot(seq(1, 20, 1), RECOMBFILE_XI[j,], ylim=c(YMIN_XI, YMAX_XI), type='n', xaxt='n', yaxt='n')
				abline(h=c(YSEQ_XI, YSEQ_XI_minor), v=c(XPOS), col='darkgray', lty="dotted") #make plot grids
				
				#below is for making an empty legend box, then placing in a plot label (i.e., 'A', 'B', 'C', 'D', etc.) that corresponds to the gradient legend and shows which generations are being plotted
				#get legend box coordinates
				my.legend.size = legend('bottomright', legend=expression(italic('plot')), col='', pch=NA, cex=empty_legendcex, box.lwd=2, bg='white') #note label 'plot' is just filler, it will not be printed
				#make an empty legend box
				legend('bottomright', legend=expression(italic('plot')), text.col='white', col='white', pch=NA, cex=empty_legendcex, box.lwd=2, bg='white') #note label 'plot' is just filler, it will not be printed
				if(MODELTYPE %in% c('6', '7', '8', '9')){
					#fill in legend box with a plot label
					text(my.legend.size$text$x+0.08, my.legend.size$text$y, 'D', cex=plotlabel_cex)				
				} else{
					#fill in legend box with a plot label
					text(my.legend.size$text$x+0.08, my.legend.size$text$y, 'E', cex=plotlabel_cex)
				}
				
				#add x-axis title						
				text(10.5, XTITLEPOS,  "Genetic Distance (cM)", xpd=NA, cex=axis_title_cex)
				box(lwd=2)
				
				axis(2, at=c(YSEQ_XI, YSEQ_XI_minor), tck=tck, labels=FALSE)
				axis(1, at=XPOS, tck=-0.025, labels=FALSE)
				text(XPOS, XLABPOS, labels=XLAB, srt=45, adj=1, pos=NULL, xpd=NA, cex=xlab_ylab_values_cex)
				
				#plot xi/xi0 across the 20 windows for new plot
				lines(seq(1, 20, 1), RECOMBFILE_XI[j,], lwd=2, col=mycolors[j])
				
				#keep track of first generation for this plot
				FIRSTJ = j
			} else{
				
				#plot each generation's xi/xi0 across the 20 windows
				lines(seq(1, 20, 1), RECOMBFILE_XI[j,], lwd=2, col=mycolors[j])
			}
		} else{
			
			#make a variable, 'newstart', to account for change in direction of xi/xi0
			#if a change in direction is encountered, make a new plot (this will prevent overlap of generations on the y-axis)
			newstart=RECOMBFILE_XI[j,1]
			if(newstart == maximumx){ #we have reached a nadir in xi/xi0, make a new plot!
				
				#redraw in dashed gray and red lines the first and last generation to make it easier to ID where the first generation begins and where the last generation ends on the previous plot
				lines(seq(1, 20, 1), RECOMBFILE_XI[1,], lwd=2, col='red', lty=2) #first generation
				lines(seq(1, 20, 1), RECOMBFILE_XI[j-1,], lwd=3, col='red', lty=3) #last generation
				FIRSTJ = j #this is the new generation to start from for plotting our new plot
				
				#start making our next plot, beginning where we left off!
				plot(seq(1, 20, 1), RECOMBFILE_XI[j,], ylim=c(YMIN_XI, YMAX_XI), type='n', xaxt='n', yaxt='n')
				abline(h=c(YSEQ_XI, YSEQ_XI_minor), v=c(XPOS), col='darkgray', lty="dotted") #make plot grids
				
				#below is for making an empty legend box, then placing in a plot label (i.e., 'A', 'B', 'C', 'D', etc.) that corresponds to the gradient legend and shows which generations are being plotted
				#get legend box coordinates
				my.legend.size = legend('bottomright', legend=expression(italic('plot')), col='', pch=NA, cex=empty_legendcex, box.lwd=2, bg='white') #note label 'plot' is just filler, it will not be printed
				#make an empty legend box
				legend('bottomright', legend=expression(italic('plot')), text.col='white', col='white', pch=NA, cex=empty_legendcex, box.lwd=2, bg='white') #note label 'plot' is just filler, it will not be printed
				#fill in legend box with a plot label
				text(my.legend.size$text$x+0.08, my.legend.size$text$y, 'C', cex=plotlabel_cex)
				
				#add x-axis title	
				text(10.5, XTITLEPOS,  "Genetic Distance (cM)", xpd=NA, cex=axis_title_cex)
				box(lwd=2)
				
				axis(2, at=c(YSEQ_XI, YSEQ_XI_minor), tck=tck, labels=FALSE)
				axis(1, at=XPOS, tck=-0.025, labels=FALSE)
				text(XPOS, XLABPOS, labels=XLAB, srt=45, adj=1, pos=NULL, xpd=NA, cex=xlab_ylab_values_cex)
				
				#plot xi/xi0 across the 20 windows for new plot
				lines(seq(1, 20, 1), RECOMBFILE_XI[j,], lwd=2, col=mycolors[j])
			} else{
				#plot each generation's xi/xi0 across the 20 windows
				lines(seq(1, 20, 1), RECOMBFILE_XI[j,], lwd=2, col=mycolors[j])
			}
		}
	}
	#redraw in dashed gray and red lines the first and last generation to make it easier to ID where the first generation begins and where the last generation ends on the previous plot
	lines(seq(1, 20, 1), RECOMBFILE_XI[FIRSTJ,], lwd=2, col='red', lty=2) #first generation
	lines(seq(1, 20, 1), RECOMBFILE_XI[(TOTAL-3),], lwd=3, col='red', lty=3) #last generation
	points(0.00984, 1.04, pch=1)
}

#make function to plot demographic model corresponding to the simulation model
models_func = function(MODELTYPE){
	ymin = 100
	ymax = 1000000
	yseq = c(100, 1000, 10000, 100000, 1000000)
	
	if(MODELTYPE %in% c('2', '3', '4', '5', '10', '11', '12')){
		if(MODELTYPE %in% c('4', '8', '5', '9')){
			xlab = c('-1.0', '-0.95', '-0.9', '-0.8', '-0.7', '-0.6', '-0.5', '-0.4', '-0.3', '-0.2', '-0.1', '0.0')
			xpos = c(2, 7, seq(12, 102, 10))
		} else{
			xlab = c('-1.0', '-0.9', '-0.8', '-0.7', '-0.6', '-0.5', '-0.4', '-0.3', '-0.2', '-0.1', '0.0')
			xpos = c(2, seq(12, 102, 10))			
		}
	} else {
		xlab = c('-0.1', '-0.09', '-0.08', '-0.07', '-0.06', '-0.05', '-0.04', '-0.03', '-0.02', '-0.01', '0.0')
		xpos = c(2, seq(4.5, 102, 2.5))
		xpos = c(1.000, 2.000, 3.000, 4.000, 5.000, 6.000, 7.000, 8.000, 9.000, 10.000, 11.000)
	}

	#read in files with demographic model information (i.e, change in N for each generation)
	demmodels1to5 = read.csv('new_generic_models_0.01gen_models1to5_and_models10to12.txt', header=TRUE)
	demmodels6to9 = read.csv('new_generic_models_0.001gen_models6to9.txt', header=TRUE)
	colnames(demmodels1to5) = c('y1', 'y2', 'y3', 'y4', 'y5', 'y10', 'y11', 'y12', 'gen', 'X')
	colnames(demmodels6to9) = c('y1', 'y6', 'y7', 'y8', 'y9', 'gen', 'X')
	
	#set an lwd variable for line widths
	w=3
	
	#make empty plot first
	if(MODELTYPE %in% c('2', '3', '4', '5', '10', '11', '12')){
		plot(demmodels1to5$X, demmodels1to5$y1, ylim=c(ymin, ymax), type='n', xaxt='n', cex=2, log='y', yaxt='n', ylab='', xlab=expression(paste('Time (', italic('N'[anc]), ' generations)', sep='')))
	} else{
		plot(demmodels6to9$X, demmodels6to9$y1, ylim=c(ymin, ymax), type='n', xaxt='n', cex=2, log='y', yaxt='n', ylab='', xlab=expression(paste('Time (', italic('N'[anc]), ' generations)', sep='')))
	}
	
	abline(h=c(yseq), v=c(xpos), col='darkgray', lty="dotted") #make plot grids
	axis(2, at=yseq, tck=-0.025, labels=FALSE)
	text(par("usr")[1]-(0.02*(par("usr")[2]-par("usr")[1])), yseq, labels=c('100', '1000', '10000', '100000', '1000000'), pos=NULL, xpd=NA, cex=0.8, adj=1)
	
	if(MODELTYPE %in% c('4', '5')){
		text(xpos+(0.019*(xpos[length(xpos)]-xpos[1])), c(25, 37, rep(25, 10)), labels=xlab, pos=NULL, xpd=NA, cex=c(0.8, 0.75, rep(0.8, 10)), adj=1)	
	} else{
		text(xpos+(0.019*(xpos[length(xpos)]-xpos[1])), 25, labels=xlab, pos=NULL, xpd=NA, cex=0.8, adj=1)	
	}
	text((par("usr")[1]+par("usr")[2])/2, 9, expression(paste('Time (', italic('N'[anc]), ' generations)', sep='')), cex=1, xpd=NA)
	if(MODELTYPE == '12'){
		text(-25, 10000, "Effective Population Size", srt=90, xpd=NA, cex=1)
		text(-21, 10000, expression(paste("(", italic("N"[e]), "; log-scale)", sep='')), srt=90, xpd=NA, cex=1)
	} else if(MODELTYPE %in% c('6', '7', '8', '9')){
		text(-1.1, 10000, "Effective Population Size", srt=90, xpd=NA, cex=1)
		text(-0.8, 10000, expression(paste("(", italic("N"[e]), "; log-scale)", sep='')), srt=90, xpd=NA, cex=1)		
	} else{
		text(-19, 10000, "Effective Population Size", srt=90, xpd=NA, cex=1)
		text(-16, 10000, expression(paste("(", italic("N"[e]), "; log-scale)", sep='')), srt=90, xpd=NA, cex=1)
	}
	axis(1, at=xpos, tck=-0.025, labels=FALSE)
	
	#cex size for plot label (i.e., model name)
	legendcex = 0.8
	
	#below is for making an empty legend box, then placing in a plot label (i.e., 'model 2', 'model 3', etc.) that corresponds to the demographic model being plotted
	if(MODELTYPE %in% c('7', '8', '9')){
		#need more room for labels for these MODELTYPEs because will have 2 digits after changing demographic model name
		#get legend box coordinates
		my.legend.size = legend('topleft', legend=paste('mid ', MODELTYPE, sep=''), text.col='white', col='white', pch=NA, cex=legendcex, box.lwd=2, bg='white') #note legend label here is just filler, it will not be printed
		#make an empty legend box
		legend('topleft', legend=paste('mod ', MODELTYPE, sep=''), text.col='white', col='white', pch=NA, cex=legendcex, box.lwd=2, bg='white') #note legend label here is just filler, it will not be printed
	} else{
		#get legend box coordinates
		my.legend.size = legend('topleft', legend=paste('mid', MODELTYPE, sep=''), text.col='white', col='white', pch=NA, cex=legendcex, box.lwd=2, bg='white') #note legend label here is just filler, it will not be printed
		#make an empty legend box
		legend('topleft', legend=paste('mod', MODELTYPE, sep=''), text.col='white', col='white', pch=NA, cex=legendcex, box.lwd=2, bg='white') #note legend label here is just filler, it will not be printed
	}
	
	#make line for the population size for the demographic model that changes in color through time so that it corresponds to the gradient legend of the pi/pi0 and xi/xi0 plots
	#and fill in legend box with the demographic model number (i.e., model label)
	if(MODELTYPE %in% c('2', '3', '4', '5')){
		color.scale.lines(demmodels1to5$X,demmodels1to5[,as.numeric(MODELTYPE)],col=mycolors_dem,lwd=w) #make population size line for model
		
		#fill in legend box with model label
		FILLTYPE = as.numeric(MODELTYPE)+3 #need to change because demographic model labels were changed during the course of manuscript writing
		text(my.legend.size$text$x+1.2, 760000, paste('model ', FILLTYPE, sep=''), cex=legendcex)
	} else if(MODELTYPE %in% c('10', '11', '12')){
		color.scale.lines(demmodels1to5$X,demmodels1to5[,as.numeric(MODELTYPE)-4],col=mycolors_dem,lwd=w) #make population size line for model
		
		#fill in legend box with model label
		FILLTYPE = as.numeric(MODELTYPE)-8 #need to change because demographic model labels were changed during the course of manuscript writing
		if(MODELTYPE == '12'){
			text(my.legend.size$text$x+2.8, 760000, paste('model ', FILLTYPE, sep=''), cex=legendcex)
		} else{
			text(my.legend.size$text$x+1.9, 760000, paste('model ', FILLTYPE, sep=''), cex=legendcex)	
		}
	} else{
		color.scale.lines(demmodels6to9$X,demmodels6to9[,as.numeric(MODELTYPE)-4],col=mycolors_dem,lwd=w) #make population size line for model
		
		#fill in legend box with model label
		FILLTYPE = as.numeric(MODELTYPE)+3 #need to change because demographic model labels were changed during the course of manuscript writing
		if(MODELTYPE %in% c('7', '8', '9')){
			text(my.legend.size$text$x+0.17, 760000, paste('model ', FILLTYPE, sep=''), cex=legendcex)
		} else{
			text(my.legend.size$text$x+0.14, 760000, paste('model ', FILLTYPE, sep=''), cex=legendcex)
		}
	}
	box(lwd=2)
}

#read in files
#only read in lines 4-total because lines 1-3 represent the first 3 generations that occur before demography begins (i.e., they are part of the burn-in)
#only read in columns 2-21 (representing 20 windows of genetic distance) because the first column is a label column
recombfile_pi = read.table(args[2], stringsAsFactors=F)[seq(4, total, 1),seq(2, 21, 1)]
recombfile_xi = read.table(args[3], stringsAsFactors=F)[seq(4, total, 1),seq(2, 21, 1)]

#y-axis limits and spacing for grid for pi/pi0 and xi/xi0 plots
ymin_pi=as.numeric(args[5])
ymax_pi=as.numeric(args[6])
yby_pi=as.numeric(args[7])
ymin_xi=as.numeric(args[8])
ymax_xi=as.numeric(args[9])
yby_xi=as.numeric(args[10])
yseq_pi=seq(ymin_pi, ymax_pi, yby_pi)
yseq_xi=seq(ymin_xi, ymax_xi, yby_xi)

modeltype=args[4]

#make positions on x-axis and labels for genetic distance markers (in cM)
xpos=seq(1, 20, 1)
xlabpos=ymin_xi-(0.1*(ymax_xi-ymin_xi))
xtitlepos=ymin_xi-(0.4*(ymax_xi-ymin_xi))
xlab=c('0.00082', '', '0.00246', '', '0.0041', '', '0.00574', '', '0.00738', '', '', '0.00984', '', '0.01148', '', '0.01312', '', '0.01476', '', '0.0164')


#below we need to assign variables (i.e., pi_break,  xi_break1, xi_break2) to make the segment markers telling us where each plot lies with respect to the color gradient legend
for(j in 1:(total-3)){
	newstart=recombfile_pi[j,1]
	if(newstart == min(recombfile_pi[,1])){
		pi_break = j
	}
}

if(modeltype %in% c('6', '7', '8', '9', '12')){
	pi_break = total-3
}

for(j in 1:(total-3)){
	newstart=recombfile_xi[j,1]
	if(modeltype == '12'){
		if(newstart == max(recombfile_xi[,1])){
			xi_break1 = j
		}
	} else{
		if(newstart == min(recombfile_xi[,1])){
			xi_break1 = j
		}
	}
	if(newstart == max(recombfile_xi[,1])){
		xi_break2 = j
	}
}

if(modeltype == '12'){
	xi_break2 = total
}

if(modeltype %in% c('10', '11', '12')){
	newmodelname = as.numeric(modeltype)-8 #need to change because demographic model labels were changed during the course of manuscript writing
} else{
	newmodelname = as.numeric(modeltype)+3 #need to change because demographic model labels were changed during the course of manuscript writing
}

#make pdf of demographic model plot
if(modeltype == '12'){
	pdf(paste('model', newmodelname, 'I.pdf', sep=''), height=3.025, width=6.25)
} else{
	pdf(paste('model', newmodelname, 'I.pdf', sep=''), height=3.04, width=8.68)
	}
par(oma=c(2.75,6,0,0) + 0.1, mar=c(0.1,0,0.7,0.7), las=1)
models_func(modeltype)
dev.off()

if(modeltype == '12'){
	WIDTH = 6.36
	HEIGHT = 6.035
} else{
	WIDTH = 8.5
	HEIGHT = 5.75
}

#make pdf of pi/pi0 and xi/xi0 plot
pdf(paste('model', newmodelname, 'II.pdf', sep=''), height=HEIGHT, width=WIDTH)
if(modeltype != '12'){
	par(mfrow=c(2,3), oma=c(7,8,0,0) + 0.1, mar=c(0,0,0.9,0.9), las=1)
} else{
	par(mfrow=c(2,2), oma=c(7,8,0,0) + 0.1, mar=c(0,0,0.9,0.9), las=1)
}
plot_func(recombfile_pi, recombfile_xi, yseq_pi, yseq_xi, ymin_pi, ymax_pi, ymin_xi, ymax_xi, xpos, xlab, xlabpos, xtitlepos, total, pi_break, xi_break1, xi_break2, modeltype)
dev.off()

#now combine demographic model plot and pi/pi0/xi/xi0 plot using your favorite pdf editor (or just powerpoint)!
