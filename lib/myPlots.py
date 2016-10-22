import matplotlib.pyplot as plt
import numpy as np

def getFig(title=None,layout=111):
    fig=plt.figure()
    if (title!=None):
        fig.canvas.set_window_title(title)
        fig.set_tight_layout('tight')
    ax=fig.add_subplot(layout)
    return fig,ax

def getLabels(ax=None,ncol=1,pos=(1,1)):
    """Get labels of plot and set legend"""
    if ax==None:
        ax=plt.gca()

    handle,labels=ax.get_legend_handles_labels()
    legend=ax.legend(handle,labels,bbox_to_anchor=pos,ncol=ncol)
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

def savePlotFile(path=None,ax=None,varx=None,vary=None,name=None):
    if (ax==None):
        ax=plt.gca()
    #n=len(ax.lines[0].get_xdata())
    if  name==None:
        name=ax.figure.canvas.get_window_title()
    if path==None:
        path='pgfPlots/'+name+'.dat'
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
                space='nan'+' '*9
                s+=(space+' \t'+space+' \t')
        s+=' \n'
    
    fh=open(path,'w')
    print('Writing to {}'.format(path))
    fh.write(s)
    fh.close()
    
