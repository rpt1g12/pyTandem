import numpy as np
from scipy import interpolate
from scipy.signal import butter, filtfilt
import matplotlib.pyplot as plt
def acorr(x):
    n=len(x)
    r=np.array([0.0 for i in range(0,n)])
    for h in range(0,n):
        for i in range(0,n-h):
            r[h]+=(x[i+h]*x[i])
        r[h]/=((n)+1)
    r[:]/=r[0]
    return r


def acorrp(x):
	n=len(x)
	r=np.array([0.0 for i in range(0,n)])
	for h in range(0,n):
		for i in range(0,n):
			if (i+h)>n-1:
				s=i+h-n
			else:
				s=i+h
			r[h]+=(x[s]*x[i])/(float(n)+1.0)
	r[:]/=r[0]
	return r

def avg(x,t):
	n=len(x)-1
	delt=np.array([0.0 for i in range(0,n+1)])
	fctr=0.5/(t[n]-t[0])
	delt[0]=fctr*(t[1]-t[0])
	delt[n]=fctr*(t[n]-t[n-1])
	for i in range(1,n):
		delt[i]=fctr*(t[i+1]-t[i-1])
	r=0.0
	for i in range(0,n+1):
		r+=delt[i]*x[i]
	return r

def find(lst, tar):
	array=abs(lst-tar)
	n=len(lst)
	for i in range(0,n):
		x=array[i]
		if x==min(array):
			result=i
	return result

def myDFT(x,time,fmax=1.0,lmf=50,mav=1,mfilt=1):
    lmt=len(x)-1
    slc=lmt//mav
    varf=np.array([0.0 for i in range(0,lmf+1)])
    freq=np.array([0.0 for i in range(0,lmf+1)])
    filt=np.array([0.0 for i in range(0,lmf+7)])
    dt=np.array([0.0 for i in range(0,slc+1)])
    vdt=np.array([0.0 for i in range(0,slc+1)])
    t=np.array([0.0 for i in range(0,slc+1)])
    nslc=2*mav-1
    mwin=min(abs(mav-1),1)
    avgp=0.0
    pi=np.pi
    print('Windowing...')
    for m in range(1,nslc+1):
        print('Averaging: %i of %i'%(m,nslc))
        ls=(m-1)*slc/2
        le=ls+slc
        period=time[le]-time[ls]
        tpop=2*pi/period
        fctr=0.5*(time[ls]+time[le])
        avgp=avgp+period
        for l in range(0,slc+1):
            t[l]=time[l+ls]-time[ls]-fctr
        dt[0]=0.5*(t[1]-t[0])
        dt[slc]=0.5*(t[slc]-t[slc-1])
        for l in range(1,slc):
            dt[l]=0.5*(t[l+1]-t[l-1])
        adt=sum(dt)/len(dt)
        print('Average dt = %f'%(adt))
        vdt=dt*x[ls:le+1]*(mwin*((np.sin(tpop*t))**2-1)+1)
        vmean=sum(vdt)/period
        vdt=vdt-vmean*t

        ra1=fmax*2*pi-tpop
        fctr=(ra1-tpop)/lmf
        ra0=tpop
        for j in range(0,lmf+1):
            freq[j]=j*fctr+ra0
            varr=0.0
            vari=0.0
            for i in range(0,slc+1):
                fnt=freq[j]*time[i]
                cosf=np.cos(fnt)
                sinf=np.sin(fnt)
                varr=varr+cosf*vdt[i]
                vari=vari-sinf*vdt[i]
            varf[j]+=varr**2+vari**2
    avgp=avgp/nslc
    tpop=2*pi/avgp
    fctr=2/(nslc*avgp)
    varf*=fctr

    ls=0
    le=lmf
    lsf=3
    lef=lmf+3
    ra0=0.5
    ra1=9.0/32.0
    ra2=0
    ra3=-1.0/32.0
    print('Filtering...')
    for m in range(1,mfilt+1):
        filt[lsf:lef+1]=varf
        filt[0:4]=varf[0:4]
        filt[lef-3:lef+1]=varf[le-3:le+1]
        for l in range(ls,le+1):
            lf=l+3
            varf[l]=( ra0*filt[lf]
            +ra1*(filt[lf-1]+filt[lf+1])
            +ra2*(filt[lf-2]+filt[lf+2])
            +ra3*(filt[lf-3]+filt[lf+3]) )

    ra0=0.5/pi
    freq*=ra0
    print('Done!')

    return freq,varf

def nextpow2(i):
    """returns next power of 2"""
    n=1
    while n<i: n*=2
    return n

def rsample(x,t,nsample=0,verbose=False,rmAvg=False):
    """docstring for rsample"""
    xspln=interpolate.splrep(t,x,s=0)
    if (nsample==0):
        nsample=len(t)
    else:
        nsample=nextpow2(nsample)
        if (nsample>len(t)):
            nsample/=2
            if (verbose): print('# of samples modified to: ',nsample)
    tnew=np.linspace(t[0],t[-1],nsample)
    fsam=1/(tnew[1]-tnew[0])
    fmax=fsam/2; fmin=fsam/nsample
    if (verbose): print('fmax:',fmax,' fmin:',fmin)
    xnew=interpolate.splev(tnew,xspln,der=0)
    if (rmAvg):
        xnew=xnew-xnew.mean()
    return xnew,tnew,nsample,fsam

def defWindows(zin,nwin=2,ovlp=0.0,plot=0,zname='signal',verbose=True):
    """docstring for defWindows"""

    if (len(zin.shape)<2):
        z=np.zeros((1,len(zin)))
        z[0,:]=zin[:]
    else:
        z=zin
        
    ntt=0;ntotal=z.shape[1]
    start=[];end=[]
    ofst=1.0-ovlp
    nt=ntotal-1
    nseg=round(float(nt)/((nwin-1)*ofst+1.0));print('nseg0=',nseg)
    while(ntt!=ntotal):
        if(nwin==1):
            ofst=1.0
        else:
            ofst=(float(nt)/nseg-1)/(nwin-1)

        nseg=round(float(nt)/((nwin-1)*ofst+1.0));

        if(nwin==1):
            ofst=1.0
        else:
            ofst=(float(nt)/nseg-1)/(nwin-1);

        ovlp=1.0-ofst
        iovlp=int(np.ceil(ovlp*nseg));iofst=nseg-iovlp

        flag=0
        while (ofst>1.0):
            nseg+=1
            ofst=(float(nt)/nseg)/(nwin-1)
            flag=1

        if(flag==1):
            if(verbose): print('Offset was more than segment length!!')

        ovlp=1.0-ofst 

        if(verbose): print('Overlap modified to ',round(ovlp,3),'to fit data')
        iovlp=int(np.ceil(ovlp*nseg));iofst=nseg-iovlp

        for i in range(nwin):
            start.append((i)*iofst)
            end.append((i)*iofst+nseg)
        ntt=1+end[-1]-start[0]
        if(ntt!=ntotal):
            nseg+=1


    if(plot>0):  
       ctitle=('ntotal='+str(ntt)+' nseg='+str(nseg)+' iovlp='+str(iovlp))
       ctitle2=('\n Final overlap of segments is '+str(round(100.0*ovlp,3))+'%')
       plt.title(ctitle+ctitle2)
       plt.ylabel(zname)
       plt.xlabel('t')
       for i in range(nwin):
           plt.plot(z[0,start[i]:end[i]+1],z[plot,start[i]:end[i]+1],label='s'+str(i+1))

       #ax1.legend(bbox_to_anchor=(0., -0.115, 1., -.15), loc=3,
       #           ncol=nwin, mode="expand", borderaxespad=0.)

       plt.show()
    fmax=1/((z[0,1]-z[0,0])*2); fmin=1/(z[0,end[0]]-z[0,start[0]])
    if(verbose):
        for i in range(nwin):
            print('Range',i+1,start[i],end[i])
    if(verbose): print('fmax:',fmax,' fmin:',fmin)
    return nseg,iovlp,ntt,fmax,fmin
    
    
def defWin(t,f,nwin=2,ovlp=0.0,zname='signal',verbose=True):
    """docstring for defWindows"""
      
    if (nwin!=1):    
          
        if(verbose):    
            fWin=plt.figure()
            fWin.canvas.set_window_title('Windowing')
            axWin=fWin.add_subplot(111)
            
        ntt=0;ntotal=len(t)
    
        ofst=1.0-ovlp
        nt=ntotal;flag=0
        nseg=round(float(nt)/((nwin-1)*ofst+1.0));print('nseg0=',nseg)
        while(ntt!=ntotal):
            start=[];end=[]
            if(flag==1):
                nseg+=1
            if(nwin==1):
                ofst=1.0
            else:
                ofst=(float(nt)/nseg-1)/(nwin-1)
    
            nseg=round(float(nt)/((nwin-1)*ofst+1.0));
    
            if(nwin==1):
                ofst=1.0
            else:
                ofst=(float(nt)/nseg-1)/(nwin-1);
    
            ovlp=1.0-ofst
            iovlp=int(np.ceil(ovlp*nseg));iofst=nseg-iovlp
    
            flag=0
            while (ofst>1.0):
                nseg+=1
                ofst=(float(nt)/nseg)/(nwin-1)
                flag=1
    
            if(flag==1):
                if(verbose): print('Offset was more than segment length!!')
    
            ovlp=1.0-ofst 
    
            if(verbose): print('Overlap modified to ',round(ovlp,3),'to fit data')
            iovlp=int(np.ceil(ovlp*nseg));iofst=nseg-iovlp
    
            for i in range(nwin):
                start.append((i)*iofst)
                end.append((i)*iofst+nseg)
            ntt=1+end[-1]-start[0]
            if(ntt!=ntotal):
                flag=1
                print(ntt,ntotal)
    
        nseg=end[0]-start[0]+1
        iovlp=iovlp+1
        if(verbose):  
           ctitle=('ntotal='+str(ntt)+' nseg='+str(nseg)+' iovlp='+str(iovlp))
           ctitle2=('\n Final overlap of segments is '+str(round(100.0*ovlp,3))+'%')
           axWin.set_title(ctitle+ctitle2)
           axWin.set_ylabel(zname)
           axWin.set_xlabel('t')
           for i in range(nwin):
               axWin.plot(t[start[i]:end[i]+1],f[start[i]:end[i]+1],label='w'+str(i+1))
    
           #ax1.legend(bbox_to_anchor=(0., -0.115, 1., -.15), loc=3,
           #           ncol=2, mode="expand", borderaxespad=0.)
    
           axWin.figure.canvas.draw()
        fmax=1/((t[1]-t[0])*2); fmin=1/(t[end[0]]-t[start[0]])
        if(verbose):
            for i in range(nwin):
                print('Range',i,start[i],end[i],end[i]-start[i]+1)
                
            print('fmax:',fmax,' fmin:',fmin)
            print('nseg='+str(nseg)+'; iovlp='+str(iovlp)+'; ntot='+str(ntt))
        
    else:
        fmax=1/((t[1]-t[0])*2); fmin=1/(t[-1]-t[0])
        ntt=len(t);iovlp=0;nseg=ntt
            
    return nseg,iovlp,ntt,fmax,fmin

def myFilter(x,cut):    
    b,a=butter(4,cut,analog=False)
    f=filtfilt(b,a,x)
    return f