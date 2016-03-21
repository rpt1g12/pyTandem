import matplotlib.pyplot as plt
import numpy as np

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
        

def fit (axes=None):
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
    xXtra[0]=-0.1*(xlim[1]-xlim[0])
    xXtra[1]=+0.1*(xlim[1]-xlim[0])
    yXtra[0]=-0.1*(ylim[1]-ylim[0])
    yXtra[1]=+0.1*(ylim[1]-ylim[0])
    print(xlim,ylim)
    axes.set_xlim(xlim+xXtra)
    axes.set_ylim(ylim+yXtra)
    axes.figure.canvas.draw()
    axes.figure.canvas.toolbar.update()
    axes.figure.canvas.toolbar.push_current()

def savePlotFile(path='plt.dat',ax=None,varx=None,vary=None):
    if (ax==None):
        ax=plt.gca()
    n=len(ax.lines[0].get_xdata())
    if (varx==None):
        m=len(ax.lines);
        varx=['x'+str(j) for j in range(m)]
        
    if (vary==None):
        m=len(ax.lines);
        vary=['y'+str(j) for j in range(m)]
    else:
        m=len(vary)
        
    daty=[col.get_ydata() for col in ax.lines]
    datx=[col.get_xdata() for col in ax.lines]
    s=''    
    for j in range(m):
        s+=varx[j]+'\t'+vary[j]+'\t' 
    s+='\n'
    n=len(ax.lines[j].get_xdata())
    for i in range(n):
        for j in range(m):
            s+='%.5f \t %.5f \t' % (datx[j][i], daty[j][i])
        s+=' \n'
    
    fh=open(path,'w')
    fh.write(s)
    fh.close()
    
