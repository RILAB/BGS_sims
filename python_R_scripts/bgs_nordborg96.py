from scipy.stats import gamma
from scipy.integrate import quad, dblquad
from math import exp, sqrt
import sys

#This script takes three arguments:
#1) Ne: the effective population size
#2) THRESH: a threshold for selection such that 2*Ne*s = THRESH, where s is the selection coefficient
#3) R: the length of the selected flanking region (in basepairs) causing background selection

#The output generated is an expectation of pi/pi_0 (i.e., B) for background selection derived from Eq. 14 of Nordborg et al. 1996, with distribution of fitness effects (i.e., distribution of s) based on those of Boyko et al. 2008 and Torgerson et al. 2009

#USAGE: python ./bgs_nordborg96.py <Ne> <THRESH> <R>
#EXAMPLE USAGE: python ./bgs_nordborg96.py 10000 4 2000000

#arguments
Ne = int(sys.argv[1]) #effective size of population
THRESH = float(sys.argv[2]) #threshold for selection
R = float(sys.argv[3]) #size of region where selected loci fall (number of basepairs)
	
def mapDist(z): #r(z) in Nordborg et al., 1996
	return z*(8.2e-10) #(rec. rate for Torres et al. 2019 simulations)

#DFE PARAMS (mean=a_SHAPE*b_SHAPE)
a_SHAPE_TORGERSON = 0.0415
b_SCALE_TORGERSON = 0.012481927710843374 
a_SHAPE_BOYKO = 0.184
b_SCALE_BOYKO = 0.15992391304347828
BOYKO_stddv = sqrt(a_SHAPE_BOYKO*b_SCALE_BOYKO**2)
BOYKO_avrg = a_SHAPE_BOYKO*b_SCALE_BOYKO
TORGERSON_stddv = sqrt(a_SHAPE_TORGERSON*b_SCALE_TORGERSON**2)
TORGERSON_avrg = a_SHAPE_TORGERSON*b_SCALE_TORGERSON

#equation 14 in Nordborg et al., 1996 (note, only integration to left of neutral locus necessary for our specific case)
#BOYKO et al. DFE
def func_integ_BOYKO_leftside(z, t):
	return 1/t*(1.0/(pow(1+mapDist(z)*(1.0-t)/t, 2.0)))*gamma.pdf(x=t, a=a_SHAPE_BOYKO, scale=b_SCALE_BOYKO)

#TORGERSON et al. DFE
def func_integ_TORGERSON_leftside(z, t):
	return 1/t*(1.0/(pow(1+mapDist(z)*(1.0-t)/t, 2.0)))*gamma.pdf(x=t, a=a_SHAPE_TORGERSON, scale=b_SCALE_TORGERSON)

TRUNC = THRESH/(2*Ne) #truncate DFE so all loci <= TRUNC are treated as neutral
U_BOYKO = (1-gamma.cdf(x=TRUNC, a=a_SHAPE_BOYKO, scale=b_SCALE_BOYKO))*2*R*0.07*1.66e-8 #total deleterious mutation rate for coding loci (BOYKO DFE) (mu = 1.66e-8)
U_TORGERSON = (1-gamma.cdf(x=TRUNC, a=a_SHAPE_TORGERSON, scale=b_SCALE_TORGERSON))*2*R*0.13*1.66e-8 #total deleterious mutation rate conserved non-coding loci (TORGERSON DFE) (mu = 1.66e-8)

P = 1.0 #proportion of deleterious loci to left of neutral locus (ALL if P = 1.0)

#get w for coding loci, integrate first over DFE starting at TRUNC and also integrate over entire region from 0 to P*R 
wBOYKO = dblquad(func_integ_BOYKO_leftside, TRUNC, BOYKO_avrg + 100*BOYKO_stddv, lambda x: 0, lambda x: P*R)[0]
wBOYKO = exp(-U_BOYKO/2.0/R*wBOYKO)

#get w for conserved non-coding loci, integrate first over DFE starting at TRUNC and also integrate over entire region from 0 to P*R 
wTORGERSON = dblquad(func_integ_TORGERSON_leftside, TRUNC, TORGERSON_avrg + 100*TORGERSON_stddv, lambda x: 0, lambda x: P*R)[0]
wTORGERSON = exp(-U_TORGERSON/2.0/R*wTORGERSON)

wTOT = wBOYKO*wTORGERSON #wTOT is pi/pi_0 (B)

print(wTOT)
