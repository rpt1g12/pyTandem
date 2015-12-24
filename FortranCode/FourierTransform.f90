program test
   implicit none
   
 integer,parameter :: nr=kind(0.0d0),lmf=800,mav=8,mfilt=1,nvar=4,nvs=0
 character(18),parameter :: csgnl='inflowsignal.dat'
 character(19),parameter :: cspct='inflowspectra.dat'

 real(nr),parameter :: pi=3.141592653589793_nr,half=0.5_nr,amach=0.4_nr
 real(nr),parameter :: turbi=1.0_nr*0.025_nr,turbl=0.04_nr,turblo=1.33898*turbl

 integer :: i,is,ie,j,k,lmt,m,mwin,stat,nslc,slc

 real(nr),dimension(:,:),allocatable :: vart,vdt
 real(nr),dimension(:),allocatable :: time,dt,tke,t

 real(nr),dimension(0:lmf,nvar) :: varf
 real(nr),dimension(0:lmf) :: freq,ek11,ek22
 real(nr),dimension(-3:lmf+3) :: filt
 real(nr),dimension(nvar) :: vmean,varr,vari
 real(nr) :: period,tpop,fctr,fnt,cosf,sinf,ra0,ra1,ra2,ra3,res,fmax,avgp,avgdt

 fmax=100_nr
 open(0,file='signal.dat')
 lmt=-2
 do while (stat.ge.0)
    read(0,*,IOSTAT=stat) 
    lmt=lmt+1
 end do
 rewind(0)
    slc=lmt/mav
    nslc=2*mav-1

 allocate(vart(0:lmt,nvar),time(0:lmt))
 allocate(vdt(0:slc,nvar),dt(0:slc),t(0:slc))

 do i = 0, lmt
    read(0,*) vart(i,1:nvar-1),time(i)
 end do
 close(0)

    vart(:,1:nvar-1)=vart(:,1:nvar-1)/amach
    time(:)=time(:)*amach

 do i = 1, nvar-1
   vmean(i)=avg(vart(:,i),time(:))
 end do

 vart(:,nvar)=((vart(:,1)-vmean(1))**2+&
        (vart(:,2)-vmean(2))**2+&
        (vart(:,3)-vmean(3))**2)*half
 
 
 mwin=min(abs(mav-1),1); varf(:,:)=0
 avgp=0.0_nr
 do m=1,nslc
       write(*,"('Averaging: ',i2,' of ',i2)") m,nslc
       is=(m-1)*slc/2
       ie=is+slc
       period=time(ie)-time(is); tpop=2*pi/period; fctr=half*(time(is)+time(ie))
       avgp=avgp+period
    do i=0,slc
       t(i)=time(i+is)-time(is)-fctr
    end do
       dt(0)=half*(t(1)-t(0)); dt(slc)=half*(t(slc)-t(slc-1))
    do i=1,slc-1
       dt(i)=half*(t(i+1)-t(i-1))
    end do
    do k=1,nvar
       vdt(:,k)=dt(:)*vart(is:ie,k)*(mwin*(sin(tpop*t(:))**2-1)+1)
       vmean(k)=sum(vdt(:,k))/period
       vdt(:,k)=vdt(:,k)-vmean(k)*dt(:)
    end do
 
       !ra0=tpop; ra1=ra0+lmf*tpop; fctr=(ra1-ra0)/lmf
       ra1=fmax*2*pi-tpop; fctr=(ra1-tpop)/(lmf); ra0=tpop
    do j=0,lmf
          freq(j)=j*fctr+ra0; varr(:)=0; vari(:)=0
       do i=0,slc
          fnt=freq(j)*time(i); cosf=cos(fnt); sinf=sin(fnt)
          varr(:)=varr(:)+cosf*vdt(i,:)
          vari(:)=vari(:)-sinf*vdt(i,:)
       end do
          varf(j,:)=varf(j,:)+varr(:)**2+vari(:)**2
    end do
 end do
 
 avgp=avgp/nslc; tpop=2*pi/avgp;

 fctr=2/((nslc)*avgp); varf(:,:)=fctr*varf(:,:)
 
 is=0; ie=lmf
 ra0=half; ra1=9.0_nr/32; ra2=0; ra3=-1.0_nr/32
 do k=1,nvar
    if(k==nvs) then
       varf(:,k)=log10(varf(:,k));
       filt(is-(/1,2,3/))=varf(is,k);
       filt(ie+(/1,2,3/))=varf(ie,k)
    end if
    do m=1,mfilt
       filt(is:ie)=varf(is:ie,k)
       if(k/=nvs) then
          filt(is-(/1,2,3/))=varf(is,k); filt(ie+(/1,2,3/))=varf(ie,k)
       end if
       do i=is,ie
          varf(i,k)=ra0*filt(i)&
          +ra1*(filt(i-1)+filt(i+1))&
          +ra2*(filt(i-2)+filt(i+2))&
          +ra3*(filt(i-3)+filt(i+3))
       end do
    end do
    if(k==nvs) then
       varf(:,k)=10**varf(:,k)
    end if
 end do


    ra0=half/(pi); ra1=(1.01325e5/(2e-5))**2
    open(0,file='spectra.dat')
 do j=0,lmf
    write(0,'(es15.7)',advance='no') ra0*freq(j)
    write(0,"(es15.7)",advance='no') varf(j,1)
    write(0,"(es15.7)",advance='no') varf(j,2)
    write(0,"(es15.7)",advance='no') varf(j,3)
    write(0,"(es15.7)") varf(j,4)
 end do
    close(0)

  contains

  function avg(x,t) result(r)
    implicit none 
    real(nr) :: r
    real(nr) :: fctr
    real(nr), dimension(0:lmt), intent(in) :: x,t
    real(nr), dimension(0:lmt) :: delt
    integer :: n,i

    n=lmt
    fctr=half/(t(n)-t(0))
    delt(0)=fctr*(t(1)-t(0))
    delt(n)=fctr*(t(n)-t(n-1))
    do i = 1, n-1
       delt(i)=fctr*(t(i+1)-t(i-1))
    end do
    r=0_nr
    do i = 0, n
       r=r+delt(i)*x(i)
    end do
  end function avg

end program test
