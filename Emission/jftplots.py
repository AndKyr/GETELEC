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

color='rgbkm'
beta=10;
Fmacmin=6e-3
Fmacmax=1.8e-1
Fmac=1/np.linspace(1/Fmacmax,1/Fmacmin,500)


Ti=[600. , 700.]
Wi=[4.5, 1.5]
betai=[66.5,10.]
R=5.0
T=500.0
gamma=10.0
ij=0
plt.figure(1)
for j in range(len(Wi)):
	W=Wi[j]
	beta=betai[j]
	F=Fmac*beta
	for i in range(len(Ti)):
		T=Ti[i]
		writeparams(F,W,R,T,gamma)
		os.system("make main")
		logJ=readJ()
		plt.semilogy(1/Fmac,np.exp(logJ),color[ij])
		ij+=1

plt.show()#plot results

