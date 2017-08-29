import numpy as np
from lib.myPlots import *
from scipy import interpolate
from scipy import fftpack
from scipy.signal import butter, filtfilt
import matplotlib.pyplot as plt

def normalise(u,v):
    m=np.sqrt(u**2+v**2)
    u=u/m
    v=v/m
    return u,v

def acorr(x):
    n=len(x)
    r=np.array([0.0 for i in range(0,n)])
    for h in range(0,n):
        for i in range(0,n-h):
            r[h]+=(x[i+h]*x[i])
        r[h]/=((n)+1)
    r[:]/=r[0]
    return r

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

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

def nextpow2(i):
    """returns next power of 2"""
    n=1
    while n<i: n*=2
    return n

def rsample(x,t,nsample=0,verbose=False,rmAvg=False,force=False,tnew=None):
    """docstring for rsample"""
    xspln=interpolate.splrep(t,x,s=0)
    if tnew==None:
        if (nsample==0):
            nsample=len(t)
            if(force):
                nsample=nextpow2(nsample)
        else:
            nsample=nextpow2(nsample-1)+1
            if (nsample>len(t) and force==False):
                nsample=(nsample-1)/2+1
                if (verbose): print('# of samples modified to: ',nsample)
        tnew=np.linspace(t[0],t[-1],nsample)
        fsam=1/(tnew[1]-tnew[0])
        fmax=fsam/2; fmin=fsam/nsample
        if (verbose): print('fmax:',fmax,' fmin:',fmin)
        xnew=interpolate.splev(tnew,xspln,der=0)
        if (rmAvg):
            xnew=xnew-xnew.mean()
        return xnew,tnew,nsample,fsam
    else:
        xnew=interpolate.splev(tnew,xspln,der=0)
        nsample=len(tnew)
        fsam=1/(tnew[1]-tnew[0])
        fmax=fsam/2; fmin=fsam/nsample
        if (verbose): print('fmax:',fmax,' fmin:',fmin)
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

def fcbFD(x,fctr=1):
    """Derivative scheme 2nd Order
    forward in left boundary
    central in middle 
    backward in right boundary"""
    n=max(x.shape); dx=np.zeros(n)
    a=np.zeros(n-2);b=np.zeros(n-1);c=np.zeros(n)
    d=np.zeros(n-1);e=np.zeros(n-2);
    a[-1]=1;e[0]=-1
    b[0:-1]=-1;b[-1]=-4;d[1:]=1;d[0]=4;
    c[0]=-3;c[-1]=3
    m=(1/(2*fctr))*np.mat(np.diag(a,-2)+np.diag(b,-1)+np.diag(c,0)+np.diag(d,1)+np.diag(e,2))
    dx=m*np.mat(x).T
    dx=np.asarray(dx.T)[0]
    return dx

def fcbFD2(x,fctr=1):
    """2st Order 2nd Derivative central finite difference"""
    n=len(x);dx=np.zeros(n)
    dx[0]=x[2]-2*x[1]+x[0]
    for i in range(1,n-1):
        dx[i]=x[i+1]-2*x[i]+x[i-1]
    dx[-1]=x[-1]-2*x[-2]+x[-3]
    if fctr != 1:
        dx/=fctr**2
    return dx

def myIntegral (u,l):
    n=len(u)
    U=np.zeros_like(u)
    for i in range(1,n):
        U[i]=0.5*(l[i]-l[i-1])*(u[i]+u[i-1])+U[i-1]
    return U

def rmvLS(x,y,o):
    """Removes the least-squares polynomial approximation of y(x) of order 'o' from y(x)"""
    p=np.polyfit(x,y,o)
    f=np.polyval(p,x)
    yn=y-f
    return yn

def fourierFilter(y,mode,fsam,cutOff,width,verbose=False):
    """Cut-off fourier filter"""
    n=len(y)
    Y=fftpack.fft(y)
    ff=np.linspace(-fsam/2,fsam/2,n)
    Y_shift=fftpack.fftshift(Y)

    flag=True
    changed=False
    while(flag):
        hWidth=int(np.ceil(width/2))
        fRlim=hWidth

        if width%2==0:
            fLlim=hWidth
        else:
            fLlim=hWidth-1

        window=np.hanning(width)
        windowR = np.hanning(width)[0:fRlim]
        windowL = np.hanning(width)[fLlim:]
        fWindow=np.ones(n)
        if mode=='high':
            fWindow[np.abs(ff)<cutOff]=0
            nL=np.where(fWindow==0)[0][0]    
            nR=np.where(fWindow==0)[0][-1]+1
            fWindow[nL:nL+hWidth]=windowL
            fWindow[nR-hWidth:nR]=windowR
            if (nL+hWidth-1>nR-hWidth):
                flag=True
                width-=1
                changed=True
            else:
                flag=False
        elif mode=='low':
            fWindow[np.abs(ff)>cutOff]=0
            nL=np.where(fWindow==1)[0][0]+1    
            nR=np.where(fWindow==1)[0][-1]
            fWindow[nR:nR+hWidth]=windowL
            fWindow[nL-hWidth:nL]=windowR
            if (nR+hWidth-1>nL-hWidth):
                flag=True
                width-=1
                changed=True
            else:
                flag=False

    Y_filt=Y_shift*fWindow
    Y_filt_shift=fftpack.ifftshift(Y_filt)
    y_filt=fftpack.ifft(Y_filt_shift)

    if verbose:
        if changed:
            print('\n###### width changed to {} ######\n'.format(width))
        rFctr=np.real(Y_shift).max()
        iFctr=np.imag(Y_shift).max()
        f,a=getFig('Frequency Filter');
        a.plot(ff,np.real(Y_shift)/rFctr,color='blue',lw=2,label='Re{y}');
        a.plot(ff,np.real(Y_filt)/rFctr,color='cyan',lw=2,linestyle='--',label='Re{y_filt}');
        a.plot(ff,np.imag(Y_shift)/iFctr,color='red',lw=2,label='Im{y}');
        a.plot(ff,np.imag(Y_filt)/iFctr,color='orange',lw=2,linestyle='--',label='Im{y_filt}');
        a.plot(ff,fWindow,color='green',lw=2,linestyle='-.',label='window');
        fit(a)
        hdl,lbl,lgd=getLabels(ax=a,ncol=3,fontsize=15)


    return np.real(y_filt)

def ROM(D,r=-1,mode=1,kfit=0):
    """Performs POD and DMD decomposition of rank r"""
    nt=len(D[0,:])-1
    A1=D[:,:-1]
    A2=D[:,1:]
    # POD decomposition
    print('Performing SVD decomposition...')
    U,s,Vt=np.linalg.svd(A1,False)
    sn=s/(s[0])
    if mode==1:
        if r==-1:
            r=len(s)
        #Build Atilde matrix (A2=A*A1)
        # Atilde=Ut*A2*V*S^(-1)
        # Atilde=Ut*Mwork
        Ur=U[:,:r]
        S=np.diag(s)[:r,:r]
        V=Vt.conj().T[:,:r]
        print('Matrix multiplication (M=A*V*S**-1)...')
        Mwork=np.dot(np.dot(A2,V),np.linalg.inv(S))
        print('Build Atilde (Atilde=Ut*M)...')
        Atilde = np.dot(Ur.conj().T, Mwork)

        print('Performing Eigen decomposition...')
        # Eigen decomposition of "Atilde"
        eV,eVec=np.linalg.eig(Atilde)

        print('Build DMD modes...')
        # Build DMD mode (space)
        # Phi=A2*V*S^(-1)*eVec
        Phi= np.dot(Mwork, eVec)
        print('Compute DMD weights...')
        # Compute mode amplitude "b"
        # A1[:,n]=Phi*b*eV^(n)
        # b=Phi^(-1)*A1[:,n]*eV^(-n)
        if kfit==-1:
            k=nt-1
        else:
            k=kfit
        b=np.dot(np.dot(np.linalg.pinv(Phi),A1[:,k]),np.diag(eV**(-k)))
        # Compute "Psi", the time variation of Phi
        Psi=np.zeros((r,nt),dtype='complex')
        for i in range(nt):
            Psi[:,i]=np.multiply(np.power(eV,i),b)

        print('Compute DMD projection on POD...')
        # DMD projection onto POD
        PhiPOD=np.zeros(r)
        for i in range(r):
            temp=[]
            for j in range(r):
                temp.append(np.dot(Phi[:,i],Ur[:,j]*sn[j]))
            PhiPOD[i]=np.linalg.norm(temp)

        return Ur,s,V,eV,eVec,Phi,Psi,PhiPOD
    elif mode==0:
        return Ur,s,V







##%% Produce log bar levels
#expM=-1;expm=-5
#levs=[10**(expM)]
#for n in range(-expM,-expm+1):
#    fctr=10**(-n)
#    for i in range(2,10):
#        levs.append(i*fctr)
#    fctr=10**(-n-1)
#    levs.append(fctr)
