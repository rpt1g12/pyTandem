# -*- coding: utf-8 -*-
import os
import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
from lib.myPlot3dOperator import *
import lib.Plot3DClasses as p3d
from scipy import stats
from scipy.interpolate import griddata
from scipy.signal import argrelextrema as locmaxmin
pi=np.pi
import getpass
user=getpass.getuser()
#from lib.matplotlib2tikz import save as tikz_save
plt.close('all')
figs=[];axs=[];hdls=[];lbls=[];lgds=[];nfig=-1;fs=18
import importlib
#%%
importlib.reload(p3d)
#%% Options
read=False
save=True;reconstruct=False
sPattern='solT*.q'
vnames=['r','u','v','w','p']; # Variable names
plane=2
autoMinMax=False
bar=True
option=1 #1 Scale with POD projection, 2 use leastDamped modes, 3 use DMD amplitude
plotOptions=2 #1 Only POD, 2 Only DMD, 3 ALL
reducedSpectrum=False
if option==1:
    usePOD=True
    useLD=False;useB=False
elif option==2:
    useLD=True
    usePOD=False;useB=False
elif option==3:
    useB=True
    useLD=False;usePOD=False
if usePOD:
    scale='POD'
elif useLD:
    scale='LD'
elif useB:
    scale='B0'
#%% Simulation Options
vname='p'

A=15 #WLE Amplitude, if SLE A=0
AoA=10 #Angle of Attack
block=0 #Block to look at
kk=0 #Spanwise slice to look at
nwave=8 #Number of LE wavelengths
xiRange=np.asarray(range(0,150));xrange=(-0.7,-0.05)
etRange=np.asarray(range(0,100));yrange=(0,0.35)
#xiRange=np.asarray(range(0,120));xrange=(-0.55,-0.25)
#etRange=np.asarray(range(0,60));yrange=(0,0.15)
zeRange=np.asarray(range(kk,kk+1))
#%% Paths set-up
if A>0:
    wavy=True
    sfolder='{}WLEk{:d}'.format(nwave,kk)
    subpath='heaving/ss005/'
else:
    wavy=False
    sfolder='{}SLEk{:d}'.format(nwave,kk)
    subpath='heaving/ss003/'

simfolder='{:1d}A{:02d}W11AoA{:02d}'.format(nwave,A,AoA)
path="/media/{}/sonyHDD/post/{}/{}".format(user,simfolder,subpath)
spath='/home/rperezt/Documents/thesis/figures/nearStall/troughROM{}2/'.format(sfolder)

dpath='/home/rperezt/Documents/thesis/data/nearStall/troughROM{}2/'.format(sfolder)
if not os.path.exists(spath) and save:
    os.makedirs(spath)
if not os.path.exists(dpath) and save:
    os.makedirs(dpath)
print('Reading data from:\n {}'.format(path))
#%%
gfile='grid.xyz' # Grid file name
files=p3d.getFileNames(path=path,pattern=sPattern)
sfile=files[0]#'solTavgCf+tw+Cp.q'# Solution file name

fl=p3d.flow(path,gfile,sfile) # Check-out flow object
fl.rdHdr() # Read header from grid file and set-up flow object properties
fl.rdGrid() # Read Grid
fl.rdSol(vnames=vnames) # Read solution file
flInfo=fl.rdFlowInfo() # Read solution file info, i.e. Mach, AoA, Re and time
mach=flInfo[0]
nt=len(files)
nxi,net,nze=fl.blk[block].size
#%%
flROM=fl.getSubsets(xRanges=[xiRange],yRanges=[etRange],zRanges=[zeRange],ssName='ROM',sspath=path+'ROM/')
temp=flROM.blk[block].var['x'].getValues()

xiD=len(xiRange)
etD=len(etRange)
xiD0=150
etD0=100
D0=np.zeros((xiD0*etD0,nt))
D=np.zeros((xiD*etD,nt))
t=np.zeros(nt)
#%%
if read:
    for n,file in enumerate(files):
        print('{:3.2f}%'.format(100*n/nt))
        fl.rdSol(vnames=vnames,sfile=file)
        t[n]=fl.rdFlowInfo(file)[3]
        for i,ii in enumerate(xiRange):
            for j,jj in enumerate(etRange):
                D0[i+j*xiD,n]=fl.blk[block].var[vname].getValues()[ii,jj,kk]
    
    D0.tofile(dpath+'Dkk{:d}.dat'.format(kk))
    t.tofile(dpath+'tkk{:d}.dat'.format(kk))
    for n,file in enumerate(files):
        temp=np.reshape(D0[:,n],(etD0,xiD0))
        for i,ii in enumerate(xiRange):
                for j,jj in enumerate(etRange):
                    D[i+j*xiD,n]=temp[jj,ii]
else:
    print('Reading previously stored data...')
    t=np.fromfile(dpath+'tkk{:d}.dat'.format(kk))
    nt=len(t)
    D0=np.reshape(np.fromfile(dpath+'Dkk{:d}.dat'.format(kk)),(xiD0*etD0,nt))
    D=D0
#    for n,file in enumerate(files):
#        print('filling {:d}th column..'.format(n))
#        temp=np.reshape(D0[:,n],(etD0,xiD0))
#        for i,ii in enumerate(xiRange):
#                for j,jj in enumerate(etRange):
#                    D[i+j*xiD,n]=temp[jj,ii]
#    
del D0
#%% Remove average
Dmean=np.zeros(xiD*etD)
Dvar=Dmean.copy()
Dnorm=D.copy()
for i in range(xiD*etD):
    Dmean[i]=np.mean(D[i,:])
    Dvar[i]=(np.var(D[i,:]))**0.5
    D[i,:]-=Dmean[i]
    Dnorm[i,:]=D[i,:]/Dvar[i]
del Dvar,Dmean
#%%
nmodes=22;
#rPro=[1,2,4,5,6,7,8,9,11,12] ;t0,t1=0,600
#rPro=[2,3,4,5,6,7,8,9,10,11,12] ;t0,t1=2800,3800;maxE=0.9 #this is good
rPro=range(2,8) ;t0,t1=2800,3800;maxE=0.9 #this is best for second ROM
#rPro=range(2,9); t0,t1=0,1728;maxE=0.90#;D=Dnorm.copy() #used in the report
   
plt.close('all')
tstep=1
print('using from t0={:d} to t1={:d}'.format(t0,t1))
tRange=range(t0,t1,tstep)
U,s,V=POD(D[:,tRange[:-1]],r=-1)
# Normalise energy
sn=s/sum(s)
i=1;tot=0
while tot<maxE:
    tot=sum(sn[:i])
    i+=1
rank=i
print('Rank is {:d}'.format(rank))
eV,eVec,Phi,Psi,PhiPOD=DMD(D[:,tRange],U,s,V,mode=1,r=rank,rPro=rPro)
#%Identify the unstable modes   
ntmode=Phi.shape[1]
unst_mode_index = []
eVNorm=np.zeros(ntmode)
for i in range(ntmode):
    temp = np.linalg.norm(eV[i])
    eVNorm[i]=temp

#% Calculate Frequency and Growth rate
tn=t[tRange][:-1]-230
dt=tstep/64.0
freqs = np.zeros((ntmode))
grate = np.zeros((ntmode))
for i in range(ntmode):
    temp = np.log(eV[i])/(dt*mach)
    freqs[i] = temp.imag/(2*pi)
    grate[i] = temp.real
    if grate[i]>=0 and freqs[i] >= 0:
        unst_mode_index.append(i)

inst_max=grate.argmax()

print ("\nUnstable modes indeces\n", unst_mode_index)
totUnst=len(unst_mode_index)
totPosFreq=len(np.asarray(np.where(freqs>=0))[0])
print ("\nNumber of unstable modes: {:d} out of {:d} total".format(totUnst,totPosFreq))  
print("{:2.2f}%".format(100*totUnst/totPosFreq))

print('Most unstable mode is:'+str(grate.argmax()))
print('Growth rate='+str(grate[inst_max]))
print('Freq='+str(freqs[inst_max]))

#% Get least damped modes
leastDamp=[]
temp=list(grate)
minG=grate.min()
while len(leastDamp)<nmodes:
    ii=np.argmax(temp)
    if freqs[ii]>=2:
        leastDamp.append(ii)
    temp[ii]=minG
print('\nLeast damped modes are:\nn\tmode\tG\tf')
for ii,i in enumerate(leastDamp):
    print('{:03d}\t{:03d}\t{:+02.4f}\t{:02.4f}'.format(ii,i,grate[i],freqs[i]))
    
#Get modes with largest POD projection
mostPOD=[]
temp=list(PhiPOD)
minPOD=min(PhiPOD)
while len(mostPOD)<nmodes:
    ii=np.argmax(temp)
    if freqs[ii]>=0:
        mostPOD.append(ii)
    temp[ii]=minPOD
print('\nModes with biggest POD projection are:\nn\tmode\tPOD\tf')
for ii,i in enumerate(mostPOD):
    print('{:03d}\t{:03d}\t{:02.4f}\t{:02.4f}'.format(ii,i,PhiPOD[i],freqs[i]))
    
#Get modes with largest DMD amplitude
b=np.abs(Psi[:,0])
bigB=[]
temp=list(b)
minPOD=min(b)
while len(bigB)<nmodes:
    ii=np.argmax(temp)
    if freqs[ii]>=0:
        bigB.append(ii)
    temp[ii]=minPOD
print('\nModes with biggest DMD amplitude are:\nn\tmode\tB\tf')
for ii,i in enumerate(bigB):
    print('{:03d}\t{:03d}\t{:02.4f}\t{:02.4f}'.format(ii,i,b[i],freqs[i]))

# Plot DMD coefficients for least-damped modes
if useLD:   
    plotRange=leastDamp
elif usePOD:
    plotRange=mostPOD
elif useB:
    plotRange=bigB
else:
    plotRange=unst_mode_index

# Plot POD energy
f,a=getFig('Sigma')
figs.append(f);axs.append(a);nfig+=1 # Append them to the figures and axes arrays
nper=np.linspace(0,1,len(sn))
axs[nfig].loglog(range(1,len(sn)+1),sn,lw=2) 
if save:
    savePlotFile(path=dpath,ax=axs[nfig])

sSum=s.copy() 
for i in range(len(s)):
    sSum[i]=sum(sn[:i+1])
f,a=getFig('SigmaSum')
figs.append(f);axs.append(a);nfig+=1 # Append them to the figures and axes arrays
axs[nfig].plot(nper,sSum,lw=2)
if not save:
    axs[nfig].axhline(y=maxE,lw=2,color='k')
if save:
    savePlotFile(path=dpath,ax=axs[nfig])

if plotOptions==1 or plotOptions==3:
    # Plot POD time coefficients
    #plt.close('all')
    nrange=len(plotRange)
    f,a=getFig('PODtime');figs.append(f),axs.append(a);nfig+=1
    for i in range(rPro[-1]+1):#rPro:
        axs[nfig].plot(tn,s[i]*V[:,i],lw=2,label='m{:02d}'.format(i))
    hdl,lbl,lgd=getLabels(ax=a,ncol=3,fontsize=15)
    hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
    fit(a)
    if save:
        savePlotFile(path=dpath,ax=axs[nfig],sameX=True)
    
    # Plot POD
    for nn in range(rPro[-1]+1):#rPro:
        if save:
            bar=False
        temp=np.reshape(U[:,nn],(xiD,etD,1),order='F')
        flROM.blk[block].setData(vname='POD{:d}'.format(nn),val=temp/np.max(temp))
        vmax=1
        vmin=-vmax
        f,a,im=flROM.blk[block].contourf(varname='POD{:d}'.format(nn),vmin=vmin,vmax=vmax,nlvl=21,cmap=plt.cm.jet,bar=bar);nfig+=1
        figs.append(f);axs.append(a)
        axs[nfig].set_xlim(xrange)
        axs[nfig].set_ylim(yrange)
        axs[nfig].set_aspect('equal')
        if save:
            saveFigOnly(path=spath,fig=figs[nfig],ax=axs[nfig],name='POD{:04d}'.format(nn),ext='.pdf') 

if plotOptions==2 or plotOptions==3:    
    # Plot DMD 
    f,a=getFig('DMDtime');figs.append(f),axs.append(a);nfig+=1
    for i in plotRange:
        axs[nfig].plot(tn,Psi[i,:].real,lw=2,label='m{:02d}R'.format(i))
    hdl,lbl,lgd=getLabels(ax=a,ncol=3,fontsize=15)
    hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
    fit(a)    
    

    for ii,nmode in enumerate(plotRange):
        if save:
            bar=False
        temp=np.reshape(Phi[:,nmode],(xiD,etD,1),order='F')
        flROM.blk[block].setData(vname='DMD{:d}'.format(nmode),val=temp/np.max(temp))
        f,a,im=flROM.blk[block].contourf(varname='DMD{:d}'.format(nmode),vmin=-1,vmax=1,nlvl=21,cmap=plt.cm.jet,bar=bar);nfig+=1
        figs.append(f);axs.append(a)
        #flROM.blk[block].contour(varname='DMD{:d}'.format(nmode),vmin=-1,vmax=1,nlvl=11,ax=a)
        if not save:    
            ftag='f={:3.2f}\ng={:3.2f}'.format(freqs[nmode],grate[nmode])  
            a.text(0.1,0.75,ftag,ha='center',va='bottom',transform=a.transAxes,fontsize=28)    
        axs[nfig].set_xlim(xrange)
        axs[nfig].set_ylim(yrange)
        axs[nfig].set_aspect('equal')
        if save:
            saveFigOnly(path=spath,fig=figs[nfig],ax=axs[nfig],name='DMD{:04d}_s{}'.format(ii,scale),ext='.pdf') 
    

    # DMD spectrum
    if reducedSpectrum:
        spectrumRange=plotRange
        extra=scale
    else:
        spectrumRange=range(rank)
        extra=''
    f,a=getFig('DMDspectrum'+extra);figs.append(f),axs.append(a);nfig+=1
    a.scatter(freqs[spectrumRange],grate[spectrumRange],s=(PhiPOD[spectrumRange]*10)**2.5,c=eVNorm[spectrumRange],cmap=plt.cm.jet)
    if reducedSpectrum:    
        a.scatter(-freqs[spectrumRange],grate[spectrumRange],s=(PhiPOD[spectrumRange]*10)**2.5,c=eVNorm[spectrumRange],cmap=plt.cm.jet)
    if not save:
        a.axhline(y=0,lw=2,linestyle='--')
    fabs=np.max(np.abs(freqs))
    #a.set_xlim(-100,100)  
    #a.set_ylim(-100,10)
    a.set_xlim(-15,15)  
    a.set_ylim(-20,5)
    if save:
        saveFigOnly(path=spath,fig=figs[nfig],ax=axs[nfig],name='DMDspectrum',ext='.pdf') 

#%% POD reconstruction
if reconstruct:
    print('\nRemove previous reconstruction...\n')
    linuxCommand='rm {}*'.format(flROM.path)
    os.system(linuxCommand)
    flROM.wrGrid()
    valROM=np.zeros((xiD,etD,1,rank))
    temp=np.zeros((xiD,etD,1))
    for nmode in range(rank):
        valROM[:,:,:,nmode]=np.reshape(U[:,nmode],(xiD,etD,1),order='F')
    for i,ifile in enumerate(tRange[:-1]):
        file=files[ifile]
        temp=0    
        for nmode in range(1,rank):
            temp+=valROM[:,:,:,nmode]*V[i,nmode]*s[nmode]
        flROM.blk[block].setData(vname='r',val=temp)
        flROM.wrSol(vnames=vnames,sfile=file)
 
#%%
print ("\nUnstable modes indeces\n", unst_mode_index)
print ("\nNumber of unstable modes: {:d} out of {:d} total".format(len(unst_mode_index),ntmode))  

print('Most unstable mode is:'+str(grate.argmax()))
print('Growth rate='+str(grate[inst_max]))
print('Freq='+str(freqs[inst_max]))

print('\nLeast damped modes are:\nn\tmode\tG\tf')
for ii,i in enumerate(leastDamp):
    print('{:03d}\t{:03d}\t{:+02.4f}\t{:02.4f}'.format(ii,i,grate[i],freqs[i]))
print('\nModes with biggest POD projection are:\nn\tmode\tPOD\tG\tf')
for ii,i in enumerate(mostPOD):
    print('{:03d}\t{:03d}\t{:02.4f}\t{:02.4f}\t{:02.4f}'.format(ii,i,PhiPOD[i],grate[i],freqs[i]))
print('\nModes with biggest DMD amplitude are:\nn\tmode\tB\tf')
for ii,i in enumerate(bigB):
    print('{:03d}\t{:03d}\t{:02.4f}\t{:02.4f}'.format(ii,i,b[i],freqs[i]))
    
print('\n{:2.1f}% of the energy is contained in ther fist {:d} modes'.format(maxE*100,rank))

#%% Resample
fsam=64;nsam=len(tn)
deltaT=dt*nsam
tO=np.floor(tn[0])

tnew=np.linspace(tO,tO+deltaT-dt,nsam)

#plt.close('all')
nw=62;ovlp=0.7;sclg='density'
nw=2;ovlp=0.9;sclg='density'
sgn,tn0,nsam,fsam=rsample(s[0]*V[:,0],tnew)
nseg,novlp,ntt,fmax,fmin=defWin(tn0,sgn,nw,ovlp,verbose=False)
rfmode=range(rPro[-1]+1)
fdata=np.zeros((len(rfmode),int(nseg/2+1)))
for nn in rfmode:#rPro:
    sgn,tn0,nsam,fsam=rsample(s[nn]*V[:,nn],tnew,verbose=True,rmAvg=True)
    ff,fdata[nn,:]=psdw(sgn,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)
    
#%Plot and save frequencies histories
st=ff/mach


f,a=getFig('freqPOD');nfig+=1
figs.append(f);axs.append(a)
for nn in rfmode:
    spsgn=fdata[nn,:]
    axs[nfig].loglog(st,spsgn,lw=2,label='m{:03d}'.format(nn))
    axs[nfig].set_xlim(0.5,100)
    
hdl,lbl,lgd=getLabels(ax=axs[nfig])
hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
if save:
    savePlotFile(path=dpath,ax=axs[nfig],sameX=True,vary=lbl)
    
#%%
save=True
plt.close('all')
tempx=flROM.blk[block].var['x'].getValues()
xN=tempx[:,-1,0]
xS=tempx[:,0,0]
xW=tempx[0,:,0]
xE=tempx[-1,:,0]
tempy=flROM.blk[block].var['y'].getValues()
yE=tempy[-1,:,0]
yW=tempy[0,:,0]
yS=tempy[:,0,0]
yN=tempy[:,-1,0]

fl.rdSol(vnames=vnames,sfile=files[t1])
f,a,im=fl.blk[block].contourf(varname='p',vmin=0.6,vmax=0.7,nlvl=21,cmap=plt.cm.hot,bar=False);nfig+=1
figs.append(f);axs.append(a)    
axs[nfig].set_xlim(-0.75,0)
axs[nfig].set_ylim(0,0.5)
axs[nfig].set_aspect('equal')
if save:
    saveFigOnly(path=spath,fig=figs[nfig],ax=axs[nfig],name='boundaries'.format(ii,scale),ext='.pdf')

axs[nfig].plot(xN,yN,label='n',lw=2)
axs[nfig].plot(xS,yS,label='s',lw=2)
axs[nfig].plot(xE,yE,label='e',lw=2)
axs[nfig].plot(xW,yW,label='w',lw=2)

hdl,lbl,lgd=getLabels(ax=axs[nfig])
hdls.append(hdl);lbls.append(lbl);lgds.append(lgd)
xlbl=['x'+c for c in lbl]
if save:
    savePlotFile(path=dpath,ax=axs[nfig],sameX=False,varx=xlbl,vary=lbl)

