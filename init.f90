module initialize
use global
use diff
USE EMPIRICAL_COEFFICIENTS
contains

!*****************************************
!To read  the dados.in file
!*****************************************
subroutine init_read()
implicit none
integer::i,j
real(kind=ip):: alpham,c1,c2,c3 
!real(kind=ip):: alpha


!................................ input data ..........................
namelist /grid/ imax, jmax, maxit,dt,sigmaf,cgd,cfd,ctd
namelist /freq/ r0_acous,amplitude,m_pt_acous,n_pt_acous,wfi,wff,wf,iterw,pxfft1,pxfft2,pxfft3,pxfft4,pxfft5
namelist /freeadm/  mmin, mmax, nmin, nmax,pxr1,pyr1,pxr2,pyr2,pxr3,pyr3
namelist /flow/ U1,U2,Re,Pr,T1,T2,Tref
namelist /grid/ D,c1b,c2b,c3b
namelist /mpi/num_procs_x,num_procs_y 
!................................ read input data ......................
open(unit=1,file='dados.in',status='old')

!................................ grid .................................
read (1,nml=grid)
write(2,nml=grid) 
!................................ grid .................................
read (1,nml=freq)
write(2,nml=freq) 

!................................ free stream conditions ...............
read (1,nml=freeadm)
write(2,nml=freeadm)

!................................ flow condicions ...............
read (1,nml=flow)
write(2,nml=flow)

!.................................PML GRIDe..............................
read (1,nml=grid)
write(2,nml=grid)
!.................................PML GRIDe..............................
read (1,nml=mpi)
write(2,nml=mpi)

close(unit=1)


end subroutine init_read

!******************************************************************

subroutine init_dados()

real(kind=ip):: Dx,Dy,Aa
real(kind=ip):: x1,x2,y1,y2
real(kind=ip):: xs,xe,ys,ye
real(kind=ip),allocatable,dimension(:):: x
real(kind=ip),allocatable,dimension(:):: y

!MESH CREATION
!Vertical  and Horizontal step, dm and dn 
!.......................................................................
!horizontal step 
 dm       =    (mmax-mmin)/real(imax-1)
!vertical step 
 dn       =    (nmax-nmin)/real(jmax-1)
  
! if ( procs_id == master ) then

    write(*,*)'dx=',dm
    write(*,*)'dy=',dn

! end if

!Point on the mesh to write de values of the variables 
!.......................................................................
 px1=int(-mmin/dm+pxr1/dm+1)          
 py1=int(-nmin/dn+pyr1/dn+1)
 px2=int(-mmin/dm+pxr2/dm+1)          
 py2=int(-nmin/dn+pyr2/dn+1)
 px3=int(-mmin/dm+pxr3/dm+1)          
 py3=int(-nmin/dn+pyr3/dn+1)


! Defini��o do sigma da BUFFER ZONE(NRBC) 
!.......................................................................
 Dx       =    dm*D
 Dy       =    dn*D
!.......................................................................
 xli= mmin+Dx
 xld= mmax-Dx
 yli= nmin+Dy
 yls= nmax-Dy

! Horizontal and vertical grid points
!.......................................................................
allocate (m(1:imax))

!movement in i, is a movement in x direction 
do i =1,imax
!do i = 1,imax 
   m(i)  = mmin + real(i-1)*dm
end do

allocate (n(1:jmax))
!allocate (n(1:jmax))

!movement in j, is a movement in y direction 
do j = 1,jmax 
   n (j) = nmin + real(j-1)*dn
end do
 
!Grid compress, only to the NRBC(non reflecting boundary condition)
!.......................................................................
Aa      = 2.d0
alpham  = 2.d0

!*********x
allocate (meshx(1:imax))
!allocate (meshx(1:imax))

meshx   = 1.d0

do i=1,imax
!do i=1,imax
   if ((i>=1).and.(i<=D+1))then
     meshx(i)   =  1.d0/(1.d0+Aa*(dabs((m(i)-(xli))/(Dx)))**alpham)
   end if

   if((i>=imax-D).and.(i<=imax))then
    meshx(i)   =  1.d0/(1.d0+Aa*(dabs((m(i)-xld)/(Dx)))**alpham)
   end if
end do

allocate (meshy(1:jmax))
!allocate (meshy(1:jmax))

meshy   = 1.d0

do j=1,jmax
!do j=1,jmax
  if ((j>=1).and.(j<=D+1))then
      meshy(j)   =  1.d0/(1.d0+Aa*(dabs((n(j)-(yli))/(Dy)))**alpham)
  end if

  if ((j>=jmax-D+1).and.(j<=jmax))then
      meshy(j)   =  1.d0/(1.d0+Aa*(dabs((n(j)-(yls))/(Dy)))**alpham)
  end if
end do

!*********************************************
!sigma buffer only to NRBC

allocate (sigmax (1:imax))
allocate (sigmaxc(1:imax))
allocate (x(1:imax))

!Initilized
sigmax  =0.d0
sigmaxc =0.d0

do i =1,imax
!do i =1,imax
   if ((i>=1).and.(i<=imaxpml))then
      xs= xli
      xe= mmin
      x(i) =(m(i)-xs)/(xe-xs)
      sigmax(i)=(1.d0-c1b*x(i)**2)*(1.d0-(1.d0-dexp(c2b*x(i)**2))/(1.d0-dexp(c2b)))
   end if 
   
   if ((i>=imax-D+1).and.(i<=imax))then
      xs= xld
      xe= mmax
      x(i) =(m(i)-xs)/(xe-xs)
      sigmax(i)=(1.d0-c1b*x(i)**2)*(1.d0-(1.d0-dexp(c2b*x(i)**2))/(1.d0-dexp(c2b)))
   end if 
end do 

allocate (sigmay (1:jmax))
allocate (sigmayc(1:jmax))
allocate (y(1:jmax))
!
sigmay  =0.d0
sigmayc =0.d0

do j =1,jmax
!do j =1,jmax
   if ((j>=1).and.(j<=jmaxpml))then
    ys= yli
    ye= nmin
    y(j) =(n(j)-ys)/(ye-ys)
    sigmay(j)=(1.d0-c1b*y(j)**2)*(1.d0-(1.d0-dexp(c2b*y(j)**2))/(1.d0-dexp(c2b)))
    !write(*,*)y(j),sigmay(j)
   end if 

   if ((j>=jmax-D+1).and.(j<=jmax))then
    ys= yls
    ye= nmax
    y(j) =(n(j)-ys)/(ye-ys)
    sigmay(j)=(1.d0-c1b*y(j)**2)*(1.d0-(1.d0-dexp(c2b*y(j)**2))/(1.d0-dexp(c2b)))
   end if 
end do 

!*********************************************

end subroutine init_dados
!****************************************************************

subroutine init_variables(U,nk)
implicit none

integer::i,j,k,nk
real(kind=ip),dimension(1:nk,1:imax,1:jmax)::U
!real(kind=ip),dimension(1:nk,1:imax,1:jmax)::U
real(kind=ip):: u1b,u2b,u3b,deltamx,delta1,delta2,psl,alpha 
real(kind=ip):: A_mi,B_mi,C_mi,D_mi
!...........................................................................................
!Allocate vicous  variables 
!ALLOCATE MEMORY TO THE VISCOSITY COEFFICIENT
ALLOCATE(MI (1:imax,1:jmax))
ALLOCATE(MIB(1:imax,1:jmax))
!ALLOCATE(MiBase(1:jmax))
MI = 0.0D0
MIB= 0.0D0
!MiBase(:) = 0.0D0
!ALLOCATE MEMORY TO THE THERMAL CONDUCTION COEFFICIENT
ALLOCATE(K1(1:imax,1:jmax))
K1 = 0.0D0
!ALLOCATE(K1B(1:imax,1:jmax))
!K1B(:,:) = 0.0D0

!ALLOCATE MEMORY TO THE TEMPERATURE AND IT'S GRADIENT
ALLOCATE(TEMP    (1,1:imax,1:jmax))
ALLOCATE(DTEMPDM (1,1:imax,1:jmax))
ALLOCATE(DTEMPDN (1,1:imax,1:jmax))
ALLOCATE(D2TEMPDM(1,1:imax,1:jmax))
ALLOCATE(D2TEMPDN(1,1:imax,1:jmax))
TEMP     = 0.0D0
DTEMPDM  = 0.0D0
DTEMPDN  = 0.0D0
D2TEMPDM = 0.0D0
D2TEMPDN = 0.0D0

!ALLOCATE MEMORY TO THE SPECIES DERIVATIVES
ALLOCATE(psiDiff    (1,1:imax,1:jmax))
ALLOCATE(DpsiFDM (1,1:imax,1:jmax))
ALLOCATE(DpsiFDN (1,1:imax,1:jmax))
ALLOCATE(D2psiFDM(1,1:imax,1:jmax))
ALLOCATE(D2psiFDN(1,1:imax,1:jmax))
PsiDiff     = 0.0D0
DpsiFDM  = 0.0D0
DpsiFDN  = 0.0D0
D2psiFDM = 0.0D0
D2psiFDN = 0.0D0

ALLOCATE(DpsiODM (1,1:imax,1:jmax))
ALLOCATE(DpsiODN (1,1:imax,1:jmax))
ALLOCATE(D2psiODM(1,1:imax,1:jmax))
ALLOCATE(D2psiODN(1,1:imax,1:jmax))
DpsiODM  = 0.0D0
DpsiODN  = 0.0D0
D2psiODM = 0.0D0
D2psiODN = 0.0D0

!ALLOCATE MEMORY TO THE TENSOR
!TAU(1,...) = TAU_XX, // TAU(2,...) = TAU_XY // TAU(3,...) = TAU_YY, 
ALLOCATE(   TAU(1:3,1:imax,1:jmax))
ALLOCATE(DTAUDM(1:3,1:imax,1:jmax))
ALLOCATE(DTAUDN(1:3,1:imax,1:jmax))

TAU   (:,1:imax,1:jmax) = 0.0D0
DTAUDM(:,1:imax,1:jmax) = 0.0D0
DTAUDN(:,1:imax,1:jmax) = 0.0D0
!...........................................................................................
!Allocate base variables 

allocate ( uba (1:jmax))
allocate ( Tb  (1:jmax))
allocate ( rhob(1:jmax))
allocate ( psib(1:jmax))

ALLOCATE(a_dudx(3)), a_dudx4(7), a_dudx5(7), a_dudx6(7), a_dudx4b(7), a_dudx5b(7), a_dudx6b(7))

a_dudx(1)= 0.770882380518d0   
a_dudx(2)= -0.166705904415d0   
a_dudx(3)= 0.0208431427703d0 

a_dudx4(1)= 0.049041958000d0   
a_dudx4(2)=-0.46884035700d0  
a_dudx4(3)=-0.47476091400d0   
a_dudx4(4)= 1.27327473700d0  
a_dudx4(5)=-0.51848452600d0
a_dudx4(6)= 0.16613853300d0
a_dudx4(7)=-0.026369431000d0 

a_dudx4b(7)=-0.049041958000d0 
a_dudx4b(6)= 0.46884035700d0  
a_dudx4b(5)= 0.47476091400d0  
a_dudx4b(4)=-1.27327473700d0  
a_dudx4b(3)= 0.51848452600d0
a_dudx4b(2)=-0.16613853300d0
a_dudx4b(1)= 0.026369431000d0

a_dudx5(1)=-0.20933762200d0    
a_dudx5(2)=-1.08487567600d0   
a_dudx5(3)= 2.14777605000d0   
a_dudx5(4)=-1.38892832200d0   
a_dudx5(5)= 0.76894976600d0 
a_dudx5(6)=-0.28181465000d0 
a_dudx5(7)= 0.048230454000d0 

a_dudx5b(7)= 0.20933762200d0    
a_dudx5b(6)= 1.08487567600d0    
a_dudx5b(5)=-2.14777605000d0    
a_dudx5b(4)= 1.38892832200d0    
a_dudx5b(3)=-0.76894976600d0  
a_dudx5b(2)= 0.28181465000d0  
a_dudx5b(1)=-0.048230454000d0

a_dudx6(1)=-2.19228033900d0    
a_dudx6(2)= 4.74861140100d0   
a_dudx6(3)=-5.10885191500d0   
a_dudx6(4)= 4.46156710400d0   
a_dudx6(5)=-2.83349874100d0 
a_dudx6(6)= 1.12832886100d0 
a_dudx6(7)=-0.20387637100d0 

a_dudx6b(7)=  2.19228033900d0    
a_dudx6b(6)= -4.74861140100d0   
a_dudx6b(5)=  5.10885191500d0   
a_dudx6b(4)= -4.46156710400d0   
a_dudx6b(3)=  2.83349874100d0 
a_dudx6b(2)= -1.12832886100d0 
a_dudx6b(1)=  0.20387637100d0 

!inicialized the variables 
 uba = 0.d0
  Tb = 0.d0
rhob = 0.d0

!parametro pml
!beta=-1.d0/c0
!...........................................................................................
!parameters to defined the  mixing layer 

!up velocity
u1b     = U1
!down velocity
u2b     = U2
!up temperature
!At: dados.in  T1      = 1.0d0
!down temperature
!At dados.in   T2      = 0.1d0
!mixing layer thickness 
deltamx = 0.4d0
!gamma=cp/vc
gamma   = 1.4d0
psi1	= 1.0d0
psi2	= 0.0d0

!CARACTERISC TEMPERATURE
!At dados.in Tref = 273.15d0 

do j=1,jmax
   uba(j)   = 0.5d0*(u1b+u2b+(u1b-u2b)*tanh(2.0d0*n(j)/deltamx))
   Tb(j)    = T1*(uba(j)-u2b)/(u1b-u2b)+T2*(u1b-uba(j))/(u1b-u2b)+(gamma-1.d0)/2.d0*(u1b-Uba(j))*(uba(j)-u2b)
   rhob(j)  = 1.d0/(Tb(j))
   psib(j)   = 0.5d0*(psi1+psi2+(psi1-psi2)*tanh(2.0d0*n(j)/deltamx))
end do

!...........................................................................................

!!!............................zerando as variaveis......................
!U = 0.d0

!!........................... inicialização ............................

do j=1,jmax
  do i=1,imax
       U(1,i,j)  = rhob(j) 
       U(2,i,j)  = uba(j)
       U(3,i,j)  = 0.d0
       U(4,i,j)  = 1.d0/gamma 
       U(5,i,j)  = psib(j) 
  end do 
end do

allocate (Ub(1:nk,1:imax,1:jmax))
ub = 0.0d0

do j=1,jmax
  do i=1,imax
    Ub(1,i,j) = rhob(j) 
    Ub(2,i,j) = uba(j)
    Ub(3,i,j) = 0.0d0
    Ub(4,i,j) = 1.0d0/gamma 
    Ub(5,i,j)  = psib(j) 

  end do 
end do

!! Allocation derivative variables 

allocate (dUdm  (nk,1:imax,1:jmax))
allocate (dUdn  (nk,1:imax,1:jmax))
allocate (dUdmb (nk,1:imax,1:jmax))
allocate (dUdnb (nk,1:imax,1:jmax))

dUdm  = 0.0d0
dUdn  = 0.0d0
dUdmb = 0.0d0
dUdnb = 0.0d0

CALL GET_COEFFICIENTS()

end subroutine init_variables


end module
