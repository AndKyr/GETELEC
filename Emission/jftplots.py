import numpy as np
import matplotlib.pyplot as plt
import os

def plotcurve():
	#read data from file
	x,y,reg= np.loadtxt("J-Fplot.csv",delimiter=',',unpack='True')	

	xf=x[np.where(reg==-1)]
	xi=x[np.where(reg==0)]
	xt=x[np.where(reg==1)]

	yf=y[np.where(reg==-1)]
	yi=y[np.where(reg==0)]
	yt=y[np.where(reg==1)]
	plt.plot(xt,yt,'r',label="ft")  # connect points with a blue line
	plt.plot(xf,yf,'b',label="ff")
	plt.plot(xi,yi,'k',label="fi")

def writeparams(F,W,R,T,gamma):
	f=open('paramin.csv','w')
	print >>f, W,R,T,gamma
	print >>f, str(list(F))[1:-1]
	f.close()

def readJ():
	x,y,reg = np.loadtxt("J-Fplot.csv",delimiter=',',unpack='True')
	return y

font=40
beta=1.;
Fmacmin=2.
Fmacmax=8.
Fmac=1./np.linspace(1./Fmacmax,1./Fmacmin,500)


Ti=[300. , 500.]
Wi=[4.5, 4, 3.5]
#betai=[66.5,10.]
R=200.0
T=00.0
gamma=10.0
ij=0
plt.figure(1)
for j in range(len(Wi)):
	W=Wi[j]
	#beta=betai[j]
	F=Fmac*beta
	if j:
		Ti=[300.]
	for i in range(len(Ti)):
		T=Ti[i]
		writeparams(F,W,R,T,gamma)
		os.system("make main")
		logJ=readJ()
		leglab='$\phi=' + str(W)+'eV$, $T='+ str(int(T)) + 'K$'
		plt.semilogy(1/Fmac,np.exp(logJ),label=leglab,linewidth=2)
		ij+=1

plt.xlabel(r'$1/F (m/GV)$',fontsize=font)
plt.ylabel(r'$J (A/nm^2)$',fontsize=font)
plt.xticks(fontsize=font)
plt.yticks(fontsize=font)
plt.legend(fontsize=font,loc=3)

plt.show()#plot results

