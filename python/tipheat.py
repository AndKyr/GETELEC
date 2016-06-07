#!/usr/bin/env python

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rcParams["font.family"] = "Serif"
matplotlib.rcParams["font.size"] = 35
matplotlib.rcParams["axes.labelsize"] = 35
matplotlib.rcParams["xtick.labelsize"] = 35
matplotlib.rcParams["ytick.labelsize"] = 35
matplotlib.rcParams["legend.fontsize"] = 35

fig = plt.figure(figsize=(20,15))
ax = fig.gca()
ax.grid()
ax.set_xlabel(r"$z [nm]$")
ax.set_ylabel(r"$heat [W]$")
ax.set_title(r"FN-plot test")

x = np.array([0.0000000000000000E+00,0.1836421299574484E+00,0.3672842599148967E+00,0.5509263898723451E+00,0.7345685198297934E+00,0.9182106497872418E+00,0.1101852779744690E+01,0.1285494909702138E+01,0.1469137039659587E+01,0.1652779169617035E+01,0.1836421299574484E+01,0.2020063429531932E+01,0.2203705559489380E+01,0.2387347689446829E+01,0.2570989819404277E+01,0.2754631949361726E+01,0.2938274079319174E+01,0.3121916209276622E+01,0.3305558339234071E+01,0.3489200469191519E+01,0.3672842599148967E+01,0.3856484729106416E+01,0.4040126859063863E+01,0.4223768989021313E+01,0.4407411118978761E+01,0.4591053248936209E+01,0.4774695378893657E+01,0.4958337508851105E+01,0.5141979638808554E+01,0.5325621768766003E+01,0.5509263898723451E+01,0.5692906028680899E+01,0.5876548158638347E+01,0.6060190288595796E+01,0.6243832418553244E+01,0.6427474548510692E+01,0.6611116678468141E+01,0.6794758808425589E+01,0.6978400938383038E+01,0.7162043068340486E+01,0.7345685198297934E+01,0.7529327328255382E+01,0.7712969458212831E+01,0.7896611588170280E+01,0.8080253718127727E+01,0.8263895848085177E+01,0.8447537978042625E+01,0.8631180108000073E+01,0.8814822237957522E+01,0.8998464367914970E+01,0.9182106497872418E+01,0.9365748627829866E+01,0.9549390757787314E+01,0.9733032887744763E+01,0.9916675017702211E+01,0.1010031714765966E+02,0.1028395927761711E+02,0.1046760140757456E+02,0.1065124353753201E+02,0.1083488566748945E+02,0.1101852779744690E+02,0.1120216992740435E+02,0.1138581205736180E+02,0.1156945418731925E+02])
y = np.array([0.2337265383759011-309,0.0000000000000000E+00,0.2738412957739241E-13,0.3105416756454558E-13,0.3105416756454511E-13,0.3105416756452439E-13,0.3105416756415570E-13,0.3105416754253219E-13,0.3105416363359640E-13,0.3105416752795989E-13,0.3105416590636309E-13,0.3105413111636850E-13,0.3105416756301921E-13,0.3105416755638241E-13,0.3105413778084647E-13,0.3105016194233809E-13,0.3096118057230434E-13,0.3005030120791851E-13,0.2215526878218883E-13,-0.6254685062071281E-14,-0.1487756429110761E-12,-0.1250061963362855E-13,0.4230822779014739E-10,0.1316964608507596E-10,0.8823518161044656E-09,0.2337895613260235E-08,0.1502443528263776E-07,0.1859073910101161E-07,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.0000000000000000E+00,0.1157654531250000E+02])

ax.plot(x,y,"b-",linewidth=2,markersize=3,label=r"$current$")
ax.set_yscale("log")

ax.legend(loc="best")

plt.savefig("png/tipheat.png")
plt.show()
