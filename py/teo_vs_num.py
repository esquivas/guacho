import numpy as np
import matplotlib.pyplot as plt
import struct
from guacho_utils import *

#path = '../prueba/grav/split_all/mhd/10G/BIN/'
path = '../prueba/grav/split_all/mhd/2temp/BIN/'
#path = '../prueba/grav/non_split/mhd/2temp/BIN/'

rhoi = readbin3d_all(nout=0,neq=0,path=path,verbose=False, mhd=True)
vxi  = readbin3d_all(nout=0,neq=1,path=path,verbose=False, mhd=True)
vyi  = readbin3d_all(nout=0,neq=2,path=path,verbose=False, mhd=True)
vzi  = readbin3d_all(nout=0,neq=3,path=path,verbose=False, mhd=True)
Pi   = readbin3d_all(nout=0,neq=4,path=path,verbose=False, mhd=True)

rhof = readbin3d_all(nout=100,neq=0,path=path,verbose=False, mhd=True)
vxf  = readbin3d_all(nout=100,neq=1,path=path,verbose=False, mhd=True)
vyf  = readbin3d_all(nout=100,neq=2,path=path,verbose=False, mhd=True)
vzf  = readbin3d_all(nout=100,neq=3,path=path,verbose=False, mhd=True)
Pf   = readbin3d_all(nout=100,neq=4,path=path,verbose=False, mhd=True)

gamma     = 5./3.
amh       = 1.66e-24    
Kb        = 1.38e-16     
Rg        = 8.3145e7     
Ggrav     = 6.67259e-8
Msun      = 1.99E33 
Rsun      = 6.955e10
GM        = Ggrav*Msun
g         = GM/Rsun/Rsun
Temp      = 1.e4
n_0       = 1.e19
mu_number = 0.6
mu_1      = 1. # Totalmente neutro
mu_2      = .5 # Totalmente ionizado

P_0 = n_0*amh*Rg*Temp

ny = 100
y_phys = 1.e9 #100 Mm
dy = y_phys/ny

#y = np.linspace(0,1,ny)
y = np.arange(5000000,1015000000,dy)

#y = y*y_phys

#Perfil de T
T   = np.zeros(y.shape)
P   = np.zeros(y.shape)
rho = np.zeros(y.shape)
mu = np.zeros(y.shape)
smut = np.zeros(y.shape)


Pt = np.zeros(y.shape)
Tempi = np.zeros(y.shape)
Tempf = np.zeros(y.shape)

for i in range(len(T)):
  mu[i] = 0.6
  ycoord = (i + 0.5)*dy
  if(ycoord<=2.e8):
    T[i]  = 1.E4
    #Tempi[i]=Pi[1,i,ny/2]*mu/rhoi[1,i,ny/2]/Rg 
    #Tempf[i]=Pf[1,i,ny/2]*mu/rhof[1,i,ny/2]/Rg
  else:
    T[i]  = 1.E6
    #Tempi[i]=Pi[1,i,ny/2]*mu/rhoi[1,i,ny/2]/Rg 
    #Tempf[i]=Pf[1,i,ny/2]*mu/rhof[1,i,ny/2]/Rg 
	

mut = 0.  
for i in range(len(T)):
#    t_invs = t_invs+(1/T[i])
    mut = mut+(mu[i]/T[i])
#    rTs[i] = t_invs	
    smut[i] = mut

Pt  = P_0*np.exp((-g*dy/Rg)*smut)
rho = mu*Pt/Rg/T
Tempi=Pi*mu_number/rhoi/Rg 
Tempf=Pf*mu_number/rhof/Rg



plt.ion()

plt.figure(1)
plt.clf()
plt.plot(y[1:],rhoi[1,0:,ny/2],'o-',label='Inicial')
plt.plot(y[1:],rhof[1,0:,ny/2],'.-',label='Final')
plt.plot(y[1:],rho[1:],'--',linewidth=2,label='Analitica')
plt.yscale('log')
plt.title('Densidad')
plt.legend()

plt.figure(2)
plt.clf()
plt.plot(y[1:],Pi[1,0:,ny/2],'o-',label='Inicial')
plt.plot(y[1:],Pf[1,0:,ny/2],'.-',label='Final')
plt.plot(y[1:],Pt[1:],'--',linewidth=2,label='Analitica')
plt.yscale('log')
plt.title('Presion')
plt.legend()

plt.figure(3)
plt.clf()
plt.plot(y[1:],Tempi[1,::,ny/2],'o-',label='Inicial')
plt.plot(y[1:],Tempf[1,::,ny/2],'.-',label='Final')
plt.plot(y[1:],T[1:],'--',linewidth=2,label='Analitica')
#plt.yscale('log')
plt.title('Temperatura')
plt.legend(loc='lower right')
