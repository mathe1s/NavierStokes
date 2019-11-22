module global

      implicit none
      integer            :: im,jm
      parameter(im=1000,jm=1200)      
      integer, parameter :: ip = selected_real_kind(15, 307) 
      integer            ::  maxit
      real(kind=ip)      :: mmin,mmax,nmin,nmax
      integer            :: imax,jmax
      real(kind=ip)      :: dm, dn, dt,pi
      parameter(pi=acos(-1.d0))

      !Mesh streching
      real(kind=ip),allocatable,dimension(:) :: meshx,d2meshx
      real(kind=ip),allocatable,dimension(:) :: meshy,d2meshy


      !Number of equations used for Runge-Kutta
      integer, parameter :: nVar = 5


!FLow parameterssss
!-------------------------------------------------------------------
      !VISCOSITY COEFFICIENT MI(I,J)
      !REAL(KIND=IP), ALLOCATABLE, DIMENSION(:,:,:) :: Ub

      !Base constants 
      real(kind=ip) :: gamma, gammaB 

      !Characteristics Scales
      REAL(KIND=IP) :: l_0, a_0, rho_0, mi_0, TEMP_0

      !Flow Namelist
      REAL(KIND=IP) :: U1,U2,T1,T2,Tref,Amplitude

      !Non dimensional numbers
      REAL(KIND=IP) :: RE,PR,ST,MA
      REAL(KIND=IP),parameter :: Pe = 1000.0d0
      REAL(KIND=IP),parameter :: alpha_cond = 0.5d0

!-------------------------------------------------------------------
      !BASE FLOW(I,J)
      REAL(KIND=IP), ALLOCATABLE, DIMENSION(:,:,:) :: Ub, UBase

      !-------------------------------------------------------------------
      !VISCOSITY COEFFICIENT MI(I,J)
      REAL(KIND=IP), ALLOCATABLE, DIMENSION(:,:) :: MI,MI2, MIB
      REAL(KIND=IP), ALLOCATABLE, DIMENSION(:) :: MiBase
      
      !Array for use in finite-differences
      REAL(KIND=IP), ALLOCATABLE, DIMENSION(:) :: a_dudx, a_dudx4, a_dudx5, a_dudx6, a_dudx4b, a_dudx5b, a_dudx6b
      
      !COEFFICIENTS TO CALCULE MI(I,J)
      REAL(KIND=IP) :: A_MI, B_MI, C_MI, D_MI, MIB0, MI0

      !VISCOUS TERMS
      !CHANGED THE THREE VARIABLES FOR ONE MULTIDIMENTIONAL MATRIX TAU, DTAUDM, DTAUDN 
      !TAU(1,...) = TAU_XX, // TAU(2,...) = TAU_XY // TAU(3,...) = TAU_YY,
      REAL(KIND=IP),ALLOCATABLE,DIMENSION(:,:,:) :: TAU, TAUB              
      REAL(KIND=IP),ALLOCATABLE,DIMENSION(:,:,:) :: DTAUDM,DTAUDN,DTAUDMB,DTAUDNB

      !THERMAL CONDUCTION COEFFICIENT K(I,J)
      REAL(KIND=IP), ALLOCATABLE, DIMENSION(:,:) :: K1,K1B

      !TO OBTAIN THE THERMAL CONDUCTION COEFFICIENT K AND REFERENCE VALUE FOR K
      REAL(KIND=IP) :: A_K, B_K, C_K, D_K, K0, KB0

      !TEMPERATURE GRADIENT IN I AND J DIRECTION
      REAL(KIND=IP),ALLOCATABLE, DIMENSION(:,:,:) :: DTEMPDM, DTEMPDN,DTEMPDMB, DTEMPDNB, D2TEMPDM, D2TEMPDN

      !SPECIES
      REAL(KIND=IP), ALLOCATABLE, DIMENSION(:,:,:) :: dPsiFdm, dPsiOdm, dPsiFdn, dPsiOdn, psiDiff
      REAL(KIND=IP), ALLOCATABLE, DIMENSION(:,:,:) :: d2PsiFdm, d2PsiOdm, d2PsiFdn, d2PsiOdn
      real(kind=ip),allocatable,dimension(:) :: psib
      REAL(KIND=IP) :: psi1, psi2

      !TEMPERATURE
      REAL(KIND=IP), ALLOCATABLE, DIMENSION(:,:,:) :: TEMP, TEMPB

      !Base varibles, only depend of u(y)
      real(kind=ip),allocatable,dimension(:) :: drhobdn,dubadn
      real(kind=ip),allocatable,dimension(:) :: uba,rhob,Tb,mibj,mib2
    
      integer :: impml,jmpml
      parameter(impml=60,jmpml=60)
      integer :: imaxpml,jmaxpml,D

      real(kind=ip):: xli,xld,yli,yls,sigmamx,delta,sigmamy
      real(kind=ip):: alpha,beta,sigmamxc,sigmamyc,c0,c1b,c2b,c3b

      !controle do filtro 
      real(kind=ip):: sigmaf

      !controle da pml
      real(kind=ip),allocatable,dimension(:) :: sigmax
      real(kind=ip),allocatable,dimension(:) :: sigmaxc
      real(kind=ip),allocatable,dimension(:) :: sigmay
      real(kind=ip),allocatable,dimension(:) :: sigmayc

      !acoutics font
      real(kind=ip),allocatable,dimension(:,:) :: s_acou

      !Frequency variables
      real(kind=ip)::r0_acous
      !Ponto de inflexao 
      real(kind=ip)::m_pt_acous,n_pt_acous
      integer::iterw
      real(kind=ip)::wf,wfi,wff
      real(kind=ip)::pxfft1,pxfft2,pxfft3,pxfft4,pxfft5 
      integer:: ix1,ix2,ix3,ix4,ix5
   

      !Size domain variables
      real(kind=ip),allocatable,dimension(:) :: m
      real(kind=ip),allocatable,dimension(:) :: n

      real(kind=ip):: alphafx,alphafy

      real(kind=ip)::pxr1,pxr2,pxr3,pyr1,pyr2,pyr3 
      integer :: px1,px2,px3,py1,py2,py3

      integer :: cg,cf,ct 
      integer :: cgd,cfd,ctd 

      !Euler equation varibles
      real(kind=ip),allocatable,dimension(:,:,:,:) :: A,B,C,Aa,Bb 

      !MPI Variables
      !real(kind=ip),dimension(4,4,im,jm)::C
      integer ndim
      integer dim_procs
      parameter(ndim =2)
      integer num_procs
      integer num_procs_x
      integer num_procs_y
      integer :: imax_procs
      integer :: jmax_procs
      integer :: procs_id
      integer :: error
      integer :: master
      integer,dimension(ndim)::dims,coords
      logical,dimension(ndim)::isperiodic
      logical:: reorder
      integer :: cartesian_comm
      integer :: cartesian_id
      integer :: cartesian_master
      parameter(cartesian_master = 0)

      integer :: dim
      
      !size of the subdomains
      !left to right
      integer ::  i_s,d_s
      !bottom to top
      integer ::  s_s,e_s
      ! ghost points to exchange bethew the subdomains
      integer ::  g_p
      parameter(g_p=6)
      !vector to save the size of the other process
      !in the master process,only the master allocate it.
      !In module mpiown.90
      integer,allocatable,dimension(:)::com_vector 

      !NEIBORHOT OF THE DIFFERENT PROCESS DEFINED BY THE SUBROUTINE 
      !MPI_CART_SHIFT (CARTESIAN_COMM, 0, SHIFT, LEFT, RIGHT, ERROR) &
      !MPI_CART_SHIFT (CARTESIAN_COMM, 1, SHIFT, BOT, TOP, ERROR)
      !LOCATION: EULER.F90
      integer ::  shift
      integer ::  bot
      integer ::  top
      integer ::  left 
      integer ::  right

      !PML q variables 
      real(kind=ip),allocatable,dimension(:,:,:) :: q1,q2
      !real(kind=ip),dimension(4,im,jm):: q1dtp0,q1dtp1,q1dtp2,q1dtp3
      !real(kind=ip),dimension(4,im,jm):: q2dtp0,q2dtp1,q2dtp2,q2dtp3

      !Derived of the variables 
      real(kind=ip),allocatable,dimension(:,:,:) :: dUdm,dUdn, dUdmb,dUdnb


      real(kind=ip),allocatable,dimension(:,:,:) :: dq1dm,dq1dn,dq2dm,dq2dn

      !Datatype
      integer :: band3type

      character(40) :: efile1,efile2,efile3

end module



