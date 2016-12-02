import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
pi=np.pi
sin=np.sin
cos=np.cos
plt.close('all')
#%%
def newboundaryF(n,gamma,b0,b1,b2):
    w=np.linspace(0,pi,n)
    A=np.zeros((3,n))
    B=A.copy()
    C=A.copy()
    D=A.copy()
    
    for t in range(n):
        A[0,t]=1+gamma[0,1]*cos(w[t])+gamma[0,2]*cos(2*w[t])
        B[0,t]=gamma[0,1]*sin(w[t])+gamma[0,2]*sin(2*w[t])
        for m in range(7):
            if (m!=0):
                C[0,t]+=b0[m]*(cos((m)*w[t])-1)
                D[0,t]+=b0[m]*sin((m)*w[t])
        A[1,t]=1+(gamma[1,0]+gamma[1,2])*cos(w[t])+gamma[1,3]*cos(2*w[t])
        B[1,t]=(gamma[1,2]-gamma[1,0])*sin(w[t])+gamma[1,3]*sin(2*w[t])
        for m in range(7):
            if (m!=1):
                C[1,t]+=b1[m]*(cos((m-1)*w[t])-1)
                D[1,t]+=b1[m]*sin((m-1)*w[t])
        A[2,t]=1+(gamma[2,1]+gamma[2,3])*cos(w[t])+(gamma[2,0]+gamma[2,4])*cos(2*w[t])
        B[2,t]=(gamma[2,3]-gamma[2,1])*sin(w[t])+(gamma[2,4]-gamma[2,0])*sin(2*w[t])
        for m in range(7):
            if (m!=2):
                C[2,t]+=b2[m]*(cos((m-2)*w[t])-1)
                D[2,t]+=b2[m]*sin((m-2)*w[t])
    return A,B,C,D
    
def nbpwave(n,i,A,B,C,D):
    rwp=[];iwp=[]
    for t in range(n):
        rwp.append((A[i,t]*D[i,t]-B[i,t]*C[i,t])/(A[i,t]*A[i,t]+B[i,t]*B[i,t])/pi)
        iwp.append(-(A[i,t]*C[i,t]+B[i,t]*D[i,t])/(A[i,t]*A[i,t]+B[i,t]*B[i,t])/pi)
    return rwp,iwp
#%%
n=1000;invpi=1/pi
w=np.linspace(0,pi,n)
nw=invpi*w
iw=np.zeros(n)

gamma00=0; gamma01=5.912678614078549;  gamma02=3.775623951744012; gamma03=0; gamma04=0;
gamma10=0.08360703307833438; gamma11=0;  gamma12=2.058102869495757; gamma13=0.9704052014790193; gamma14=0;
gamma20=0.03250008295108466; gamma21=0.3998040493524358;  gamma22=0; gamma23=0.7719261277615860; gamma24=0.1626635931256900;

gamma=np.array([
    [gamma00, gamma01, gamma02, gamma03, gamma04],
    [gamma10, gamma11, gamma12, gamma13, gamma14],
    [gamma20, gamma21, gamma22, gamma23, gamma24]])
    
b00=0; b01=-3.456878182643609; b02=5.839043358834730; b03=1.015886726041007;
b04=-0.226526470654333; b05=0.08564940889936562; b06=-0.01836710059356763;

b10=-0.3177447290722621; b11=0; b12=-0.02807631929593225; b13=1.593461635747659;
b14=0.2533027046976367; b15=-0.03619652460174756; b16=0.004080281417108407;

b20=-0.1219006056449124; b21=-0.6301651351188667; b22=0; b23=0.6521195063966084;
b24=0.3938843551210350; b25=0.01904944407973912; b26=-0.001057260523947668;

b0 = [b00, b01, b02, b03, b04, b05, b06];
b1 = [b10, b11, b12, b13, b14, b15, b16];
b2 = [b20, b21, b22, b23, b24, b25, b26];
#%%
A,B,C,D=newboundaryF(n,gamma,b0,b1,b2)
fi,ar=getFig('boundary_Scheme_R')
fr,ai=getFig('boundary_Scheme_I')
for i in range(3):
    rwp,iwp=nbpwave(n,i,A,B,C,D)
    ar.plot(nw,rwp,label='i={:d}'.format(i))
    ai.plot(nw,iwp,label='i={:d}'.format(i))
#%%
save=1
ax=ai
if save:
    handle,labels,legend=getLabels(ax=ax)
    savePlotFile(ax=ax)