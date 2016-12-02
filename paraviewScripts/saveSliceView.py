#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
#paraview.simple._DisableFirstRenderCameraReset()
import os
import numpy as np
import getpass
user=getpass.getuser()

# get active view
view = GetActiveViewOrCreate('RenderView')
# Set parallel projection
view.CameraParallelProjection = 1
# current camera placement for view
view.CameraPosition = [-2.15430058882296, 0.0528290681540966, 4.72033740400665e-16]
view.CameraFocalPoint = [-0.00748521089553833, 0.0528290681540966, 0.0]
view.CameraParallelScale = 0.16974135720845
# uncomment following to set a specific view size
# view.ViewSize = [1461, 836]

# get slice object
slc=FindSource('Slice4')

# Initial slice postion
x0=-0.48
# End slice postion
x1=-0.37
# Number of slices
nx=12
# set slice positions
x=np.linspace(x0,x1,nx)

print '%d slices from %5.2f to %5.2f' % (nx,x0,x1)

for i in range(nx):
    xp=x[i]
    slc.SliceType.Origin=[xp,0,0]
    RenderAllViews()
    name='x%-2d' % (xp*100)
    print 'saving slice: '+name
    # export view
    path='/home/'+user+'/Desktop/jfmpaper/figures/fig14/'+name+'.pdf'
    ExportView(path, view=view, Plottitle='ParaView GL2PS Export',
        Compressoutputfile=0,
        Drawbackground=0,
        Cullhiddenprimitives=1,
        Linewidthscalingfactor=1.0,
        Pointsizescalingfactor=1.0,
        GL2PSdepthsortmethod='Simple sorting (fast, good)',
        Rasterize3Dgeometry=0,
        Dontrasterizecubeaxes=1,
        Rendertextaspaths=1)

