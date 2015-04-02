from guachoUtils import *

def add_inner_title(ax, title, loc, size=None, **kwargs):
    from matplotlib.offsetbox import AnchoredText
    from matplotlib.patheffects import withStroke
    if size is None:
        size = dict(size=plt.rcParams['ytick.labelsize'])
    at = AnchoredText(title, loc=loc, prop=size,
                      pad=0., borderpad=0.5,
                      frameon=False, **kwargs)
    ax.add_artist(at)
    #at.txt._text.set_path_effects([withStroke(foreground="w", linewidth=4)])
    return at

nxtot=600
nytot=150
nztot=600
mpiX=6
mpiY=1
mpiZ=6
neqs=7


for nout in range(0,1) :

    #map=coplot3d(2,nytot/2,0,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path='../M1/',endian='<')
    #map3d=readbin3d(0,neqs,nxtot,nytot,nztot,mpiX,mpiY,mpiZ,nout,path='../M1/',endian='<')

    images=[map,map3d[nxtot/2,::,::],map3d[::,nytot/2,::],map3d[::,::,nztot/2]]
    maxV=1.E10#map3d.max()
    minV=100. #map3d.min()


    F = plt.figure(1, (6, 6))
    F.clf()
    grid = ImageGrid(F, 111,
                 nrows_ncols = (2, 2),
                 direction="row",
                 axes_pad = 0.0,
                 add_all=True,
                 label_mode = "1",
                 share_all = True,
                 cbar_location="top",
                 cbar_mode="single",
                 cbar_size="1.8%",
                 cbar_pad=".5%" )

    extent=(-25,25,-25,25)
        
    for i in range(0,4):

        
        im=grid[i].imshow(images[i],cmap='jet',origin='lower',vmin=minV,vmax=maxV
                      ,extent=extent , norm=LogNorm() )

    lvls = np.logspace(np.log10(minV),np.log10(maxV),5)
    grid[2].cax.colorbar(im, extend='none',norm=LogNorm() , ticks=lvls )
    grid[2].cax.toggle_label(True)

    grid[0].cax.text(550.,9.5,r'Density [cm$^{-3}$]')
    grid[0].cax.text(75000.,9.5,r't='+str(nout*5000.).format('f')+' yr')

    for ax, im_title in zip(grid, [r"(a) XZ", r"(b) YZ",r"(c) XZ",r"(d) XY"]):
        t = add_inner_title(ax, im_title,loc=2 )
        t.patch.set_ec("none")
        t.patch.set_alpha(1.0)

#fig = plt.gcf()
#fig.set_size_inches(8.,8.)

#    plt.savefig('./fig-M1-'+str(nout).zfill(3)+'.pdf',dpi=300,transparent=True,bbox_inches='tight')
#plt.savefig('./fig.eps',dpi=300,transparent=True,bbox_inches='tight',facecolor='gray')

plt.draw()
plt.show()
