import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import matplotlib as mpl
import os

## Custom colormaps ##

kbw = LinearSegmentedColormap.from_list('kbw', [(0,0,0), (0,0,0.501961), (0,0.501961,1), (1,1,1)])
kbw_r = LinearSegmentedColormap.from_list('kbw_r', [(1,1,1), (0,0.501961,1), (0,0,0.501961), (0,0,0)])

## Funtions ##

def getFig(title=None,layout=111):
    fig=plt.figure()
    if (title!=None):
        fig.canvas.set_window_title(title)
        fig.set_tight_layout('tight')
    ax=fig.add_subplot(layout)
    return fig,ax
    
def getFig3d(title=None,layout=111):
    fig=plt.figure()
    if (title!=None):
        fig.canvas.set_window_title(title)
        fig.set_tight_layout('tight')
    ax=fig.add_subplot(layout,projection='3d')
    return fig,ax

def getLabels(ax=None,ncol=1,fontsize=16,loc='best',sep=0,hspace=0.2):
    """Get labels of plot and set legend"""
    if ax==None:
        ax=plt.gca()

    handle,labels=ax.get_legend_handles_labels()
    legend=ax.legend(handle,labels,ncol=ncol,fontsize=fontsize,
                    loc=loc,borderaxespad=sep,columnspacing=hspace)
    return handle,labels,legend

def addAx(fig=None,layout=111):
    if(fig==None):
        fig=plt.gcf()
    ax=fig.add_subplot(layout)
    return ax
    
def cll(ax=None):
    if (ax==None):
        ax=plt.gca()
    n=len(ax.lines)    
    for i in range(n):
        ax.lines.pop(0)

def clt(ax=None):
    if (ax==None):
        ax=plt.gca()
    n=len(ax.texts)    
    for i in range(n):
        ax.texts.pop(0)
        
def axShow (axes=None):
    if (axes==None):
        axes=plt.gca()
    axes.figure.canvas.draw()
    axes.figure.show()
    
def fit (axes=None,spg=(0.0,0.2)):
    if (axes==None):
        axes=plt.gca()
    i=0
    xlim=np.zeros(2);ylim=xlim.copy()
    xXtra=np.zeros(2);yXtra=xXtra.copy()
    for l in axes.lines:
        x=l.get_xdata();y=l.get_ydata()
        if (i==0):
            xlim[0]=min(x)
            ylim[0]=min(y)
            xlim[1]=max(x)
            ylim[1]=max(y) 
        i+=1
        xlim[0]=min(xlim[0],min(x))
        ylim[0]=min(ylim[0],min(y))
        xlim[1]=max(xlim[1],max(x))
        ylim[1]=max(ylim[1],max(y))
    xXtra[0]=-spg[0]*(xlim[1]-xlim[0])
    xXtra[1]=+spg[0]*(xlim[1]-xlim[0])
    yXtra[0]=-spg[1]*(ylim[1]-ylim[0])
    yXtra[1]=+spg[1]*(ylim[1]-ylim[0])
    print(xlim,ylim)
    axes.set_xlim(xlim+xXtra)
    axes.set_ylim(ylim+yXtra)
    axes.figure.canvas.draw()
    axes.figure.canvas.toolbar.update()
    axes.figure.canvas.toolbar.push_current()

def savePlotFile(path=None,ax=None,varx=None,vary=None,name=None,nan=True,sameX=False):
    if (ax==None):
        ax=plt.gca()
    #n=len(ax.lines[0].get_xdata())
    if  name==None:
        name=ax.figure.canvas.get_window_title()+'.dat'
    if path==None:
        path='pgfPlots/'
    if not os.path.exists(path):
        os.makedirs(path)
    if (varx==None):
        m=len(ax.lines);
        varx=['x'+str(j) for j in range(m)]
    if (vary==None):
        handle,labels,legend=getLabels(ax)
        if len(labels)==0:
            m=len(ax.lines);
            vary=['y'+str(j) for j in range(m)]
        else:
            vary=labels
            m=len(vary)
    else:
        m=len(vary)
        
    daty=[col.get_ydata() for col in ax.lines]
    if sameX:
        datx=ax.lines[0].get_xdata()
        s=varx[0]+'\t'
        nn=[]
        for j in range(m):
            s+=vary[j]+'\t'
            nn.append(len(ax.lines[j].get_xdata()))
        s+='\n'
        n=max(nn)
        for i in range(n):
            s+='%e \t' % (datx[i])
            for j in range(m):
                if (i<len(ax.lines[j].get_xdata())):
                    s+='%e \t' % (daty[j][i])
                else:
                    if nan:
                        space='nan'+' '*9
                    else:
                        space=' '*12
                    s+=(space+' \t')
            s+=' \n'

    else:
        datx=[col.get_xdata() for col in ax.lines]
        s='' 
        nn=[]
        for j in range(m):
            s+=varx[j]+'\t'+vary[j]+'\t'
            nn.append(len(ax.lines[j].get_xdata()))
        s+='\n'
        n=max(nn)
        for i in range(n):
            for j in range(m):
                if (i<len(ax.lines[j].get_xdata())):
                    s+='%e \t %e \t' % (datx[j][i], daty[j][i])
                else:
                    if nan:
                        space='nan'+' '*9
                    else:
                        space=' '*12
                    s+=(space+' \t'+space+' \t')
            s+=' \n'
    
    path=path+name
    fh=open(path,'w')
    print('Writing to {}'.format(path))
    fh.write(s)
    fh.close()
    
def saveFigOnly(path=None,fig=None,ax=None,name=None,dpi=600,ext='.png'):
    if (fig==None):
        fig=plt.gcf()
    if (ax==None):
        ax=plt.gca()
    if  name==None:
        name=ax.figure.canvas.get_window_title()+ext
    else:
        name+=ext
    if path==None:
        path='pgfPlots/'
        
    ax.axes.get_yaxis().set_visible(False)
    ax.axes.get_xaxis().set_visible(False)
    extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(path+name,bbox_inches='tight', pad_inches=0,dpi=dpi)

#def saveColorBarOnly(path=None,im=None,fig=None,ax=None,name=None,hv=0,dpi=600,ext='.png'):
#    """docstring for saveColorBarOnly"""
#    if (im==None):
#        im=plt.gci()
#    if (fig==None):
#        fig=plt.gcf()
#    if (ax==None):
#        ax=plt.gca()
#    if  name==None:
#        name=ax.figure.canvas.get_window_title()+ext
#    else:
#        name+=ext
#    if path==None:
#        path='pgfPlots/'
#
#    if hv==0:
#        orientation='horizontal'
#    elif hv==1:
#        orientation='vertical'
#        
#    ax.set_visible(False)
#    cbar=fig.colorbar(im,orientation=orientation)
#    cbar.ax.yaxis.set_visible(False)
#    fig.savefig(path+name,bbox_inches='tight', pad_inches=0,dpi=dpi)
#    pass
    

class Streamlines(object):
    def __init__(self, X, Y, U, V, res=0.125, spacing=4, maxLen=50, detectLoops=True):
        """
        Compute a set of streamlines covering the given velocity field.

        X and Y - 1D or 2D (e.g. generated by np.meshgrid) arrays of the
                  grid points. The mesh spacing is assumed to be uniform
                  in each dimension.
        U and V - 2D arrays of the velocity field.
        res - Sets the distance between successive points in each
              streamline (same units as X and Y)
        spacing - Sets the minimum density of streamlines, in grid points.
        maxLen - The maximum length of an individual streamline segment.
        detectLoops - Determines whether an attempt is made to stop extending
                      a given streamline before reaching maxLen points if
                      it forms a closed loop or reaches a velocity node.

        Plots are generated with the 'plot' or 'plotArrows' methods.
        """

        self.spacing = spacing
        self.detectLoops = detectLoops
        self.maxLen = maxLen
        self.res = res

        xa = np.asanyarray(X)
        ya = np.asanyarray(Y)
        self.x = xa if xa.ndim == 1 else xa[0]
        self.y = ya if ya.ndim == 1 else ya[:,0]
        self.u = U
        self.v = V
        self.dx = (self.x[-1]-self.x[0])/(self.x.size-1) # assume a regular grid
        self.dy = (self.y[-1]-self.y[0])/(self.y.size-1) # assume a regular grid
        self.dr = self.res * np.sqrt(self.dx * self.dy)

        # marker for which regions have contours
        self.used = np.zeros(self.u.shape, dtype=bool)
        self.used[0] = True
        self.used[-1] = True
        self.used[:,0] = True
        self.used[:,-1] = True

        # Don't try to compute streamlines in regions where there is no velocity data
        for i in range(self.x.size):
            for j in range(self.y.size):
                if self.u[j,i] == 0.0 and self.v[j,i] == 0.0:
                    self.used[j,i] = True

        # Make the streamlines
        self.streamlines = []
        while not self.used.all():
            nz = np.transpose(np.logical_not(self.used).nonzero())
            # Make a streamline starting at the first unrepresented grid point
            self.streamlines.append(self._makeStreamline(self.x[nz[0][1]],
                                                         self.y[nz[0][0]]))

    def plot(self, lw=1, ax=None):
        """
        Draw the computed streamline segments.

        Optional keyword arguments:
            lw - line width
            ax - axes to use for plotting

        """
        if ax is None:
            ax = plt.gca()

        for streamline in self.streamlines:
            ax.plot(streamline[0], streamline[1], 'k', lw=lw)

        ax.axis('tight')

        return ax

    def plotArrows(self, lw=1, ax=None, size=16):
        """
        Draw the computed streamline segments with arrows indicating flow direction.

        Optional keyword arguments:
            lw - line width
            size - size of the arrow head
            ax - axes to use for plotting

        """
        if ax is None:
            ax = plt.gca()

        for streamline in self.streamlines:
            path = mpl.path.Path(np.asarray((streamline[0], streamline[1])).T)
            patch = mpl.patches.FancyArrowPatch(path=path, arrowstyle='->',
                                                mutation_scale=size, lw=lw)
            ax.add_patch(patch)

        ax.axis('tight')

        return ax

    def _interp(self, x, y):
        """ Compute the velocity at point (x,y) """
        i = (x-self.x[0])/self.dx
        ai = i % 1

        j = (y-self.y[0])/self.dy
        aj = j % 1

        # Bilinear interpolation
        u = (self.u[j,i]*(1-ai)*(1-aj) +
             self.u[j,i+1]*ai*(1-aj) +
             self.u[j+1,i]*(1-ai)*aj +
             self.u[j+1,i+1]*ai*aj)

        v = (self.v[j,i]*(1-ai)*(1-aj) +
             self.v[j,i+1]*ai*(1-aj) +
             self.v[j+1,i]*(1-ai)*aj +
             self.v[j+1,i+1]*ai*aj)

        self.used[j:j+self.spacing,i:i+self.spacing] = True

        return u,v

    def _makeStreamline(self, x0, y0):
        """
        Compute a streamline extending in both directions from the given point.
        """

        sx, sy = self._makeHalfStreamline(x0, y0, 1) # forwards
        rx, ry = self._makeHalfStreamline(x0, y0, -1) # backwards

        rx.reverse()
        ry.reverse()

        return rx+[x0]+sx, ry+[y0]+sy

    def _makeHalfStreamline(self, x0, y0, sign):
        """
        Compute a streamline extending in one direction from the given point.
        """

        xmin = self.x[0]
        xmax = self.x[-1]
        ymin = self.y[0]
        ymax = self.y[-1]

        sx = []
        sy = []

        x = x0
        y = y0
        i = 0
        while xmin < x < xmax and ymin < y < ymax:
            u, v = self._interp(x, y)
            theta = np.arctan2(v,u)

            x += sign * self.dr * np.cos(theta)
            y += sign * self.dr * np.sin(theta)
            sx.append(x)
            sy.append(y)

            i += 1

            if self.detectLoops and i % 10 == 0 and self._detectLoop(sx, sy):
                break

            if i > self.maxLen / 2:
                break

        return sx, sy

    def _detectLoop(self, xVals, yVals):
        """ Detect closed loops and nodes in a streamline. """
        x = xVals[-1]
        y = yVals[-1]
        D = np.array([np.hypot(x-xj, y-yj)
                      for xj,yj in zip(xVals[:-1],yVals[:-1])])
        return (D < 0.9 * self.dr).any()

