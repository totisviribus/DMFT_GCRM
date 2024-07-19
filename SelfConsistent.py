import cmath
import numpy as np
from scipy.optimize import minimize
from scipy.integrate import quad
import random as rd
import sys
####

gamma_inv=float(sys.argv[1])
K=float(sys.argv[2])
sigK=0.10
h=float(sys.argv[3])
muC=1.00
m=1.00
sigM=0.10
D=1.00
sigD=0.00

########
def PN(N,g,H,sigN2):
	result = ( H / np.sqrt(2.*np.pi*sigN2) ) * np.exp( - (H*N - g)**2/(2.*sigN2) ) 
	return result

def NormN(x, sigC):
	Phi,N,N2,Nu,R,R2,Chi = x

	g = muC*R - m
	sigC2 = sigC**2
	sigM2 = sigM**2
	sigN2 = sigM**2 + sigC2*R2
	H = h + sigC2*Chi

	Norm, err = quad(lambda X: PN(X,g,H,sigN2) , 0., np.inf )
	return Norm

#####
## Independent variables : Phi, N, N2, Nu, R, R2, Chi
def Consumers(x, sigC):
	Phi,N,N2,Nu,R,R2,Chi = x

	g = muC*R - m
	sigC2 = sigC**2
	sigM2 = sigM**2
	sigN2 = sigM**2 + sigC2*R2
	sigN = np.sqrt(sigN2)
	H = h + sigC2*Chi

	PPhi = NormN(x, sigC)

	NN, errN = quad(lambda X: X * PN(X,g,H,sigN2) , 0., np.inf )

	NN2, errN2 = quad(lambda X: X**2 * PN(X,g,H,sigN2) , 0., np.inf )

	NNu, errNu = quad(lambda X: 1./np.sqrt(2.*np.pi*sigN2) * np.exp( - (H*X - g)**2/(2.*sigN2) ) , 0., np.inf )

	return PPhi, NN, NN2, NNu

######
def PR(R,Deff,A,sigK2,sigR2):
	result=1./np.sqrt(2.*np.pi)*np.abs( ( (A*R**2 + K)*sigR2*R + sigK2*(2.*A*R + Deff) ) / ( sigR2*R**2 + sigK2 )**1.50 ) * np.exp( - ( A*R**2 + Deff*R - K )**2 / ( 2.*(sigR2*R**2 + sigK2) ) )
	return result

def NormR(x, sigC):
	Phi,N,N2,Nu,R,R2,Chi = x

	Deff = D + gamma_inv*muC*N
	sigD2 = sigD**2
	sigC2 = sigC**2
	sigK2 = sigK**2
	sigR2 = sigD2 + gamma_inv*sigC2*N2
	A = gamma_inv*sigC2*Nu

	Norm, err = quad(lambda X: PR(X,Deff,A,sigK2,sigR2), 0., np.inf)

	return Norm
#####
## Independent variables : Phi, N, N2, Nu, R, R2, Chi
def Resources(x, sigC):
	Phi,N,N2,Nu,R,R2,Chi = x

	Deff = D + gamma_inv*muC*N
	sigD2 = sigD**2
	sigC2 = sigC**2
	sigK2 = sigK**2
	sigR2 = sigD2 + gamma_inv*sigC2*N2
	A = gamma_inv*sigC2*Nu

	RR, errR = quad(lambda X: X * PR(X,Deff,A,sigK2,sigR2), 0., np.inf)

	RR2, errR2 = quad(lambda X: X**2 * PR(X,Deff,A,sigK2,sigR2), 0., np.inf)

	CChi, errChi = quad(lambda X: X / np.sqrt( 2.*np.pi*(sigR2*X**2 + sigK2) )  * np.exp( - ( A*X**2 + Deff*X - K )**2 / ( 2.*(sigR2*X**2 + sigK2) ) ), 0., np.inf)

	return RR, RR2, CChi

########
########
def OptFtn(x,sigC):
	Phi,N,N2,Nu,R,R2,Chi = x
	PPhi, NN, NN2, NNu = Consumers(x, sigC)
	RR, RR2, CChi = Resources(x, sigC)

	opt_Phi = ( 1. - Phi/PPhi )**2
	opt_N   = ( 1. - N/NN )**2
	opt_N2  = ( 1. - N2/NN2 )**2
	opt_Nu  = ( 1. - Nu/NNu )**2
	opt_R   = ( 1. - R/RR )**2
	opt_R2  = ( 1. - R2/RR2 )**2
	opt_Chi = ( 1. - Chi/CChi )**2

	out=[opt_Phi, opt_N, opt_N2, opt_Nu, opt_R, opt_R2, opt_Chi]
	return sum(out)

##########
def main():
	sigC_List=list(reversed([i/500 for i in range(0,600)]))
	wfname="./g%.2lf_K%.2lf_h%.2lf.txt" %(gamma_inv,K,h)

## Randomly Choose Initial Start Point: Phi, N, N2, Nu, R, R2, Chi
#	bounds
	bnds=((0.0, 1.0), (0.0, np.inf), (0.0, np.inf), (0.0, np.inf), (0.0, np.inf), (0.0, np.inf), (0.0, np.inf))
#	initial
	x = np.array([rd.uniform(0,1),	rd.uniform(0,1),	rd.uniform(0,1),	rd.uniform(0,1),	rd.uniform(0,1),	rd.uniform(0,1),	rd.uniform(0,1)])

#
	for sigC in sigC_List:
		res = minimize(OptFtn, x, method='Nelder-Mead', args=(sigC), tol=1e-30, bounds=bnds, options={'disp':False, 'adaptive':False})

		x = res.x
		phi, n, n2, nu, rr, r2, chi = res.x
		err = res.fun

		wf=open(wfname,'a')
		wf.write("%.5lf\t\t%.6le\t\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\t%.8le\n" %(sigC,err,phi,n,n2,nu,rr,r2,chi))
		wf.close()
######
if __name__ == "__main__":
	main()

