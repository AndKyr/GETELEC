import numpy as np
import matplotlib.pyplot as plt
import os

font=37

def writeparams(F,W,R,T,gamma):
	f=open('paramin.csv','w')
	print >>f, W,R,T,gamma, F
	f.close()

def readspect():
	data= np.genfromtxt("spectra.csv",delimiter=',')
	E=data[:,0]
	integE=data[:,1]
	integ=data[:,2]
	return (E,integE,integ)



Fi=[0.1, 0.2, .5, 3., 8.]
W=4.5
gamma=10.0
T=300.0
R=100.0

plt.figure(1)
for F in Fi:
	writeparams(F,W,R,T,gamma)
	os.system("make spectroscopy")
	(E,integE,integ)=readspect()
	if round(F)==F:
		Fstr=str(int(round(F)))
	else:
		Fstr=str(F)
	plt.plot(E,integ/max(integ),linewidth=2,label=r'$F=$' +Fstr+ r'$GV/m$')

	
plt.text(-0.3,0.03,r'$E_F$',fontsize=font)
plt.plot([0.,0.],[0.,1.],'k',linewidth=2)
plt.xlabel(r'$E (eV)$',fontsize=font)
plt.ylabel(r'$j(E) (a.u.)$',fontsize=font)
plt.xticks(fontsize=font)
plt.yticks(fontsize=font)
plt.axis([-2.5,4.5,0.,1.],linewidth=2)
plt.legend(fontsize=font,loc=2)


plt.show()#plot results
