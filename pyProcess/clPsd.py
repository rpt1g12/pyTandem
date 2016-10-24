import numpy as np
from scipy.signal import welch as psdw
import matplotlib.pyplot as plt
from lib.stats import *
from lib.myPlots import *
from scipy import stats
pi=np.pi
#from lib.matplotlib2tikz import save as tikz_save
#plt.close('all')
#%%
folder='6blocks/';sclcd='SCl'
osim=8;aoa=20;iaoah=20;freq=0.19
if(osim==0):
    sim='A00'
    odata=0
    prdl=1
elif(osim==80):
    sim='8A00'
    odata=1
    prdl=1
elif(osim==40):
    sim='4A00'
    odata=1
    prdl=1
elif(osim==2):
    sim='A15'
    odata=0
    prdl=0
elif(osim==3):
    sim='3A15'
    odata=0
    prdl=0
elif(osim==4):
    sim='A4A15'
    odata=0
    prdl=0
elif(osim==41):
    sim='4A15'
    odata=0
    prdl=0
elif(osim==8):
    sim='8A15'
    odata=1
    prdl=0
    
sim+='W11AoA'+str(aoa)

dataset='/home/rperezt/anaconda3/pyTandem/clData/'+folder+sim+'.dat';
if(odata==1):    
    n,tin,clin,cdin,taoa,tmach=np.loadtxt(dataset,skiprows=1,unpack=True)
elif(odata==0):
    n,tin,clin,cdin=np.loadtxt(dataset,skiprows=1,unpack=True);taoa=np.array([20.0 for i in range(len(n))]);tmach=np.array([0.3 for i in range(len(n))]);
tdum,clha,cdha=np.loadtxt('/home/rperezt/anaconda3/pyTandem/clData/HansenClCd.dat',skiprows=1,unpack=True)
idum=np.where(tdum==iaoah)[0][0]
clh=clha[idum];cdh=cdha[idum]
#clh=0.54;cdh=0.31


#%%
#fOriginal=plt.figure()
#fOriginal.canvas.set_window_title('Original')
#axOriginal=fOriginal.add_subplot(111)


fTime,axTime=getFig('Time '+sim)
name=sim.split('W')[0]


show=False;showAoA=False;getAvg=False
plt.sca(axTime)
axTime.lines.clear()
cll();clt();ns=1024


#tmin=95.0/0.3
tmin=tin[-1]-(30/0.3)       #in jfm paper its -(25/0.3)
n0=np.where(tin>tmin)[0][0]
t0=tin[n0:]*0.3#-tin[n0];
if prdl==1:
    fctr=np.sqrt(1-0.3**2)
else:
    fctr=1
cl0=clin[n0:]*fctr;
cd0=cdin[n0:]*fctr
taoa0=taoa[n0:]
cln,tn,nsam,fsam=rsample(cl0,t0,verbose=True,nsample=ns)
cdn,tn,nsam,fsam=rsample(cd0,t0,nsample=ns)
dtaoa,tn,nsam,fsam=rsample(taoa0,t0,nsample=ns)
s=np.sqrt(np.var(cln))
clM=cln.mean()
cds=np.sqrt(np.var(cdn))
cdM=cdn.mean()


if (show):
    axTime.plot(tin[:n0],clin[:n0],'b');lClin=plt.gca().lines[-1]
    axTime.plot(tin[:n0],cdin[:n0],'r');lCdin=plt.gca().lines[-1]
    axTime.axvline(x=tmin,color='green',linewidth=2,linestyle='--')
    
if (showAoA):
    ax2=axTime.twinx()
    ax2.plot(tn,dtaoa,color='lime',linewidth=2,label=r'$\alpha$')
    ax2.set_ylabel(r'$\alpha(t^*)$',fontsize=20)
    
axTime.plot(tn,cln,'b',linewidth=2,label='Cl'+name)
axTime.plot(tn,cdn,'r',linewidth=2,label='Cd'+name)
axTime.plot([tn[0],tn[-1]],[clh,clh],'k',linewidth=2,label=r'exp $\alpha=$'+str(iaoah))
axTime.plot([tn[0],tn[-1]],[cdh,cdh],'k',linewidth=2)

if (getAvg):
    axTime.plot([tn[0],tn[-1]],[clM,clM],'g-.',linewidth=2,label='avg');lClnm=plt.gca().lines[-1]
    axTime.plot([tn[0],tn[-1]],[clM-s,cln.mean()-s],'r-.',linewidth=2,label='s0');lCls0=plt.gca().lines[-1]
    axTime.plot([tn[0],tn[-1]],[clM+s,clM+s],'r-.',linewidth=2,label='s1');lCls1=plt.gca().lines[-1]
    axTime.plot([tn[0],tn[-1]],[cdM,cdM],'g-.',linewidth=2,label='avg');lCdnm=plt.gca().lines[-1]
    axTime.plot([tn[0],tn[-1]],[cdM-cds,cdn.mean()-cds],'r-.',linewidth=2,label='s0');lCds0=plt.gca().lines[-1]
    axTime.plot([tn[0],tn[-1]],[cdM+cds,cdM+cds],'r-.',linewidth=2,label='s1');lCds1=plt.gca().lines[-1]


axTime.set_xlabel(r'$t^*$',fontsize=20)
axTime.set_ylabel(r'$C_l$',fontsize=20)
handle,labels=axTime.get_legend_handles_labels()
legend=axTime.legend(handle,labels,bbox_to_anchor=(1,1),ncol=3)
scl=r'$\sigma_L=$'+str(np.around(s,3))+'\t'+r'$\bar{C_L}=$'+str(np.around(clM,3))
scd=r'$\sigma_D=$'+str(np.around(cds,3))+'\t'+r'$\bar{C_D}=$'+str(np.around(cdM,3))
axTime.text(0.25,-0.15,scl,ha='center',va='bottom',transform=axTime.transAxes)
axTime.text(0.75,-0.15,scd,ha='center',va='bottom',transform=axTime.transAxes)
#xstart=np.where(tn>45)[0][0];xend=np.where(tn>165)[0][0]
#ystart=cdn.min();yend=cln.max()
#axTime.fill_between(tn[xstart:xend], ystart, yend, facecolor='gray', alpha=0.5)
#xstart=np.where(tn>195)[0][0];xend=np.where(tn>314.8)[0][0]
#ystart=cdn.min();yend=cln.max()
#axTime.fill_between(tn[xstart:xend], ystart, yend, facecolor='gray', alpha=0.5)
fit(axTime)

if (showAoA):
    fit(ax2,spg=(0,0))


#%%

scale=False;sclg='density'
nw=4;ovlp=0.5;

if sclcd=='SCl':
    sgnl=cln
else:
    sgnl=cdn

fFreq,axFreq=getFig('Frequency '+sim)

nseg,novlp,ntt,fmax,fmin=defWin(tn,sgnl,nw,ovlp,verbose=False)
#sgnl=myFilter(sgnl,0.25/(fmax))
ff,pcl=psdw(sgnl,fs=fsam,nperseg=nseg,noverlap=novlp,scaling=sclg)
if (scale):
    pcl/=np.var(cln)

plt.sca(axFreq)
axFreq.lines.clear()
#cll();clt()
sin20=np.sin(np.deg2rad(20))
#st=(sin20/0.3)*ff
st=(sin20)*ff
stdif=abs(st-freq)
sts=np.argmin(stdif)
axFreq.loglog(st,pcl,label=sclcd)

axFreq.grid(b=True, which='major', color='k', linestyle='--')
axFreq.grid(b=True, which='minor', color='k', linestyle=':')
axFreq.set_xlabel(r'$St$',fontsize=20)
axFreq.set_ylabel(r'$'+sclcd+'$',fontsize=20)
axFreq.axvline(x=st[sts],color='black',linewidth=2,linestyle='--')
sfreq=r'$f^*=$'+'{:4.2f}'.format(st[sts])+'->{:8.3e}'.format(pcl[sts])
axFreq.text(0.25,-0.15,sfreq,ha='center',va='bottom',transform=axTime.transAxes,fontsize=16)
fit(axFreq)
axFreq.set_xlim(0,1)
#%%
save=False;
name=sim.split('W')[0]+sclcd
if (save==True):
    plt.sca(axFreq)
    path='pgfPlots/'
    savePlotFile(path=path+name+'.dat',vary=[sclcd],varx=['f'])

#%%
#fig,ax=getFig('vsAoA')
#tmin=0.0#tin[-1]-(0.0*0.3)
#n0=np.where(tin>tmin)[0][0]
#t0=tin[n0:]-tin[n0];taoa0=taoa[n0:]
#taoa2,tn,nsam,fsam=rsample(taoa0,t0,verbose=True,nsample=ns)
#idum=np.where(tdum==10)[0][0]
#idum2=np.where(tdum==20)[0][0]+1
#clhn=clha[idum:idum2];aoah=tdum[idum:idum2]
#
#n1=np.where(taoa2>13)[0][0]
#ms,iss,rs,_,_=stats.linregress(taoa2[:n1],cln[:n1])
#n2=np.where(aoah>12)[0][0]
#mh,ih,rh,_,_=stats.linregress(aoah[:n2],clhn[:n2])
#ys=cln.copy()
#yh=clhn.copy()
#for i in range(len(ys)):
#    ys[i]=ms*taoa2[i]+iss
#for i in range(len(yh)):
#    yh[i]=mh*aoah[i]+ih
#ax.plot(taoa2,cln,'b',linewidth=2,label='Sim')
#ax.plot(aoah,clhn,'r--o',linewidth=2,label='Exp')
#ax.plot(taoa2,ys,'k-.',linewidth=2,label='Sim slope='+'{:4.3f}'.format(ms))
#ax.plot(aoah,yh,'k-.',linewidth=2,label='Exp slope='+'{:4.3f}'.format(mh))
#handle,labels=ax.get_legend_handles_labels()
#legend=ax.legend(handle,labels,bbox_to_anchor=(1,1),ncol=2)
#ax.set_xlim(10,20)
#ax.set_ylabel(r'$C_L$',fontsize=20)
#ax.set_xlabel(r'$\alpha$',fontsize=20)
#print(ms,mh)
#line=ax.axvline(x=taoa2[0],color='green',linewidth=2,linestyle='--')
#
#axShow(ax)
#%%
#lineT=axTime.axvline(x=85,color='green',linewidth=2,linestyle='--')
#for i in range(641):
#for i in [286]:
#    ts=(85+0.125*i-45)/120
#    tt=85+0.125*i
#    rdum=10+5*(1-np.cos(pi*ts))
#    adum=[rdum,rdum]
#    line.set_xdata(adum)
#    adum=[tt,tt]
#    lineT.set_xdata(adum)
#    axShow(ax)
#    axShow(axTime)
#    #jobplt.pause(0.01)
    #fig.savefig('/home/rpt1g12/Desktop/post/4A15W11AoA10/hplane2/pics/img.'+'{:04d}'.format(i)+'.png')