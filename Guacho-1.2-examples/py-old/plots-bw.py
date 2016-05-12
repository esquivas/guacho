from guachoUtils import *
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.colors import LogNorm


mpl.rc('xtick',labelsize=8)
mpl.rc('ytick',labelsize=8)
mpl.rc('text',usetex=True)
#mpl.rc('font',**{'family':'serif','serif':['Computer Modern']})

def add_inner_title(ax, title, loc, size=None, **kwargs):
    from matplotlib.offsetbox import AnchoredText
    from matplotlib.patheffects import withStroke
    if size is None:
        size = dict(size=plt.rcParams['ytick.labelsize'])
    at = AnchoredText(title, loc=loc, prop=size,
                      pad=0., borderpad=0.5,
                      frameon=False, **kwargs)
    ax.add_artist(at)
    at.txt._text.set_path_effects([withStroke(foreground="w", linewidth=2)])
    return at


nxtot=400
nytot=4
nztot=4
mpiX=16
mpiY=1
mpiZ=1
neqs=8


cv=1.5
gam=(cv+1.)/cv
vsc=1.
psc=1.

firstrun=True
plt.ion()


path='/Users/esquivel/Desktop/Storage-Diable/MHD-EXO/BW-HLLE/BIN/'
path1='/Users/esquivel/Desktop/Storage-Diable/MHD-EXO/BW-HLLD/BIN/'

x_ax=np.linspace(0,1,nxtot)

for nout in range(5,6):
    #  this is the ionization rate
    if (firstrun):
        denT= coplot3d(2,nytot/2,0,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path)
        vx=   coplot3d(2,nytot/2,1,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path)/denT
        vy=   coplot3d(2,nytot/2,2,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path)/denT
        vz=   coplot3d(2,nytot/2,3,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path)/denT        
        Etot= coplot3d(2,nytot/2,4,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path)
        bx=   coplot3d(2,nytot/2,5,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path)
        by=   coplot3d(2,nytot/2,6,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path)
        bz=   coplot3d(2,nytot/2,7,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path)
        Pgas=  (Etot-0.5*(denT*(vx*vx+vy*vy+vz*vz)+bx*bx+by*by+bz*bz))/cv

        denT1= coplot3d(2,nytot/2,0,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path1)
        vx1=   coplot3d(2,nytot/2,1,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path1)/denT1
        vy1=   coplot3d(2,nytot/2,2,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path1)/denT1
        vz1=   coplot3d(2,nytot/2,3,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path1)/denT1
        Etot1= coplot3d(2,nytot/2,4,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path1)
        bx1=   coplot3d(2,nytot/2,5,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path1)
        by1=   coplot3d(2,nytot/2,6,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path1)
        bz1=   coplot3d(2,nytot/2,7,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path=path1)
        Pgas1=  (Etot1-0.5*(denT1*(vx1*vx1+vy1*vy1+vz1*vz1)+bx1*bx1+by1*by1+bz1*bz1))/cv

    plt.clf()

    plt.subplot(3, 3, 1)
    plt.plot(x_ax,denT[1,::],'')
    plt.plot(x_ax,denT1[1,::],'')

    plt.ylabel(r'$\rho$')
    plt.xlabel(r'$x$')
    plt.ylim((0.,1.1))

    plt.subplot(3, 3, 2)
    plt.plot(x_ax,Pgas[1,::],'')
    plt.plot(x_ax,Pgas1[1,::],'')
    plt.ylabel(r'$P_{gas}$')
    plt.xlabel(r'$x$')
    plt.ylim((0.,1.1))

    plt.subplot(3, 3, 3)
    plt.plot(x_ax,Pgas[1,::]/denT[1,::],'')
    plt.plot(x_ax,Pgas1[1,::]/denT1[1,::],'')
    plt.ylabel(r'$P_{th}/\rho$')
    plt.xlabel(r'$x$')
    plt.ylim((0.25,1.75))
    
    plt.subplot(3, 3, 4)
    plt.plot(x_ax,vx[1,::],'')
    plt.plot(x_ax,vx1[1,::],'')
    plt.ylabel(r'$v_{x}$')
    plt.xlabel(r'$x$')
    plt.ylim((-0.4,0.8))
    
    plt.subplot(3, 3, 5)
    plt.plot(x_ax,vy[1,::],'')
    plt.plot(x_ax,vy1[1,::],'')
    plt.ylabel(r'$v_{y}$')
    plt.xlabel(r'$x$')
    plt.ylim((-1.8,0.2))
    
    plt.subplot(3, 3, 6)
    plt.plot(x_ax,vz[1,::],'')
    plt.plot(x_ax,vz1[1,::],'')
    plt.ylabel(r'$v_{z}$')
    plt.xlabel(r'$x$')
    plt.ylim((-0.4,0.8))
    
    plt.subplot(3, 3, 7)
    plt.plot(x_ax,bx[1,::],'')
    plt.plot(x_ax,bx1[1,::],'')
    plt.ylabel(r'$B_{x}$')
    plt.xlabel(r'$x$')
    plt.ylim((-1.2,1.2))

    plt.subplot(3, 3, 8)
    plt.plot(x_ax,by[1,::],'')
    plt.plot(x_ax,by1[1,::],'')
    plt.ylabel(r'$B_{y}$')
    plt.xlabel(r'$x$')
    plt.ylim((-1.2,1.2))

    plt.subplot(3, 3, 9)
    plt.plot(x_ax,bz[1,::],'')
    plt.plot(x_ax,bz1[1,::],'')
    plt.ylabel(r'$B_{z}$')
    plt.xlabel(r'$x$')
    plt.ylim((-1.2,1.2))

'''
    plt.subplot(3, 3, 1)
    plt.plot(x_ax,denT1[1,::],'')
    plt.ylabel(r'$\rho$')
    plt.xlabel(r'$x$')
    plt.ylim((0.9,1.7))

    plt.subplot(3, 3, 2)
    plt.plot(x_ax,Pgas1[1,::],'')
    plt.ylabel(r'$P_{gas}$')
    plt.xlabel(r'$x$')
    #plt.ylim((0.,1.1))

    plt.subplot(3, 3, 3)
    plt.plot(x_ax,Pgas1[1,::]/denT1[1,::],'')
    plt.ylabel(r'$P_{th}/\rho$')
    plt.xlabel(r'$x$')
    #plt.ylim((0.5,2.5))
    
    plt.subplot(3, 3, 4)
    plt.plot(x_ax,vx1[1,::],'')
    plt.ylabel(r'$v_{x}$')
    plt.xlabel(r'$x$')
    #plt.ylim((-0.4,0.8))
    
    plt.subplot(3, 3, 5)
    plt.plot(x_ax,vy1[1,::],'')
    plt.ylabel(r'$v_{y}$')
    plt.xlabel(r'$x$')
    #plt.ylim((-1.8,0.2))
    
    plt.subplot(3, 3, 6)
    plt.plot(x_ax,vz1[1,::],'')
    plt.ylabel(r'$v_{z}$')
    plt.xlabel(r'$x$')
    #plt.ylim((-0.4,0.8))
    
    plt.subplot(3, 3, 7)
    plt.plot(x_ax,bx1[1,::],'')
    plt.ylabel(r'$B_{x}$')
    plt.xlabel(r'$x$')
    #plt.ylim((-1.2,1.2))

    plt.subplot(3, 3, 8)
    plt.plot(x_ax,by1[1,::],'')
    plt.ylabel(r'$B_{y}$')
    plt.xlabel(r'$x$')
    #plt.ylim((-1.2,1.2))

    plt.subplot(3, 3, 9)
    plt.plot(x_ax,bz1[1,::],'')
    plt.ylabel(r'$B_{z}$')
    plt.xlabel(r'$x$')
    #plt.ylim((-1.2,1.2))

'''


    

#plt.ioff()

