
import math
import sys
print(sys.version_info)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.font_manager as fm
import matplotlib.axis
from matplotlib.pyplot import gca
from matplotlib.ticker import FormatStrFormatter
from matplotlib import cm
import matplotlib.colors as colors

def info(tuparr):
    for arr in tuparr:
        print(arr.min(),arr.mean(),arr.max())
#
#------------------
#  STAR FILE   #
#-----------------
A2 = open('./stars_ran_20.dat', 'r')
lines = A2.readlines()
A2.close()

data0 = [line.split() for line in lines] 
datos  = np.asfarray(data0)
x = datos[:,0]
y = datos[:,1]
z = datos[:,2]

#~~~~~~~~
# CIRCULO
#~~~~~~~~
e=0.                  
rmax=3.e17           
rmin=(1.-e)*rmax
pos_angle=0. #90.-65.   
theta = np.deg2rad(np.arange(0.0, 360.0, 1.0))
x_circ= rmax * np.cos(theta)                      #Parameterized equation of ellipse
y_circ = rmin * np.sin(theta)
rtheta = np.radians(pos_angle)
xprime = x_circ*np.cos(rtheta) - y_circ*np.sin(rtheta)  #Rotate to desired position angle
yprime = x_circ*np.sin(rtheta) + y_circ*np.cos(rtheta)

#-----------
#   Figura
#-----------
plt.figure(figsize=(11,11.))
plt.plot(x,z,'o',color='blueviolet',markersize=7,alpha=0.7)


plt.plot(xprime,yprime,'-',color='gray',linewidth=1.5)
#
plt.axis([-3.e17,3e17, -3e17,3e17])
#
plt.show()



































