!***************************************************************************** 
!*                                                                           *
!* NONLINEAR EULER EQUATION-AEROACUSTIC- 2D PARALLED DOMAIN DECOMPOSITION    *           
!*                                                                           *
!*   Program to resolv the non linear euler equation using a two dimensional *
!*   domain decomposition and resolv using MPI.                              *
!*                                                                           *
!*    Autor: Jhonatan Andres Aguirre Manco                                   * 
!*    Date : 10-04/2017                                                      *
!*                                                                           *
!*                                                                           *
!*                                                                           *
!* UPDATE            : 10-04/2018                                            *
!* MODIFICATION      : SPECIES EQUATION FOR FUTURE COMBUSTION                * 
!*                     MATHEUS SILVA                                         *
!*                                                                           *
!*                                                                           *
!* UPDATE            : 20-08/2013                                            *
!* MODIFICATION      : IMPLEMENTAÇÃO DA PML NAS  EQUAÇÕES NAO LINEARES       *
!*                     DE EULER EM VARIAVEIS PRIMITIVAS.                     *
!*                                                                           *
!* DATE              : 16-08/2013                                            *
!* MODIFY  BY        : JHONATAN                                              *
!* BASED ON          : PML ABC FOR NONLINEAR EULER EQUATIONS IN              *
!*                     PRIMITIVE VARIABLES,AIAA 2009-6 - HU, LIN, LI         *
!* GRID              : GENERAL GRID(M,N)                                     *
!* ALGORITHMS        : EXPLICIT METHOD- COMPLETE EULER EQUATION              *
!*                     2D-STRECHING NA PML EN X E Y, FILTRO 10 ORDEM         *
!* RESOLT METHOD     : RUNGE-KUTA 4 ORDER                                    *
!* BOUNDARY CONDITION:                                                       * 
!* FORTRAN 90        : MODULE , THOMA'S ALGORITHM                            *
!*                                                                           *
!*****************************************************************************

program EULER

!module with the global variables 
use global
!module to inicialized all the variables, this read 'dados.in' 
use initialize
!module to do a Runge-Kutta 
use rk
!module to arranged the results 
use outt
!module to aplied a numerical filter 
use filtering


!DEFINITIONS
implicit none
integer:: iter,icount,i,j,f,tama
real(kind=ip)::Tw,tec,tecif
real(kind=ip)::timei,timef
real(kind=ip)::t,dw
!Kinetic energy
real(kind=ip),dimension(10)::ecx
real(kind=ip),dimension(im):: ec, decdm, ekk, ekka, ekkb
!files
character(40) :: temporalfile1,temporalfile2,temporalfile3
!ALLOCATABLE
!........................................................................
!VARIABLES, U(1,:,:)=rho, U(2,:,:)=u, U(3,:,:)=v, U(4,:,:)=p
!Completed domain 
real(kind=ip),allocatable,dimension(:,:,:):: U
!........................................................................
!character(40) :: filename1,filename2,filename3,filename4,filename5,filename6,filename7
write (temporalfile1,'("temp1_D60.dat")') 
open(unit=1000,file=temporalfile1,status='unknown')

write (temporalfile2,'("temp2_D60.dat")') 
open(unit=1001,file=temporalfile2,status='unknown')
    
write (temporalfile3,'("temp3_D60.dat")') 
open(unit=1002,file=temporalfile3,status='unknown')

!INITILIZE DADES FOR EACH PROCESS
!Location: initialize  module, init.f90
  call init_read()
!........................................................................
!Print a message using the master process
  print*,  "Tref", Tref 
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'NAVIER STOKES MPI - Master process:'
  write ( *, '(a)' ) '  FORTRAN90 serial version'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  A program to solve the Non Linear Navier Stokes equations.'
  write ( *, '(a)' ) ' '

 !Subroutine to count the time of execution of the program
       call CPU_TIME(timei)
 !........................................................................
 !Contadores
 !........................................................................
 !Contador temporal
        icount = 0
 !Contador Graficas
        cg=0
 !Contador filtro 
        cf=0
!...........................................................................
!  Using to inizilized some variables, and defining the size of the mesh 
!  Location: initialize module init.f90
   call init_dados()
!...........................................................................
!exp: a(1:5)= 1,2,3,4,5---a(5:10)=5,6,7,8,9,10
   !for all domain  
   allocate (U(1:nVar,1:imax,1:jmax))
   !Base flow
   !allocate (UBase(1:4,1:imax,1:jmax))
!...........................................................................
! Fonte acustica para gerar o pulso de pressão que instabiliza o escoamento  
   allocate (s_acou(imax,jmax))
!   allocate (s_acou(i_s:d_s,s_s:e_s))
!...........................................................................
! Initial condition for the solution. It is in inicial module
! Location: initialize module init.f90

   call init_variables(U,size(U,1))
!   call init_variables(U_procs,size(U_procs,1))
 
!Loop temporal begining, using a runge kutta solution

    do iter = 0,maxit
!.........................................................................
!Real time
    t=(iter+1)*dt
!.........................................................................
!Subroutine to make the temporal iteration,using runge kuta, it is rk module 
!Location: rk module  rk.f90
    call rk_euler(U,t,iter,size(U,1))
!    call rk_euler(U_procs,t,iter,size(U_procs,1))
!.........................................................................
      if(cg==cgd)then
!.........................................................................
!Subroutine to created the output files
!Location out.f90 module
          call output(U,icount,f,size(U,1))
!          call output(U_procs,U,icount,f,size(U_procs,1))
!Write in the screen to control the iteration process, only the master process
          if(procs_id==master)then
              write(*,*)t,U(2,px1,py1),U(4,px1,py1)
!.........................................................................
!to calculated the time to make cgd iterations
              call CPU_TIME(timef) 
              write(*,*)'Tempo :',timef-timei,'segundos'
!.........................................................................
          endif
          icount = icount+1
          cg = 0
!.........................................................................
!to calculated the time to make cgd iterations
      end if
!.........................................................................
!If to write a temporal file in 3 differents  points a each ctd 
      
if(ct==ctd)then
    write(1000,*)t,U(1,px1,py1),U(2,px1,py1),U(3,px1,py1),U(4,px1,py1)!,U(5,px1,py1),U(6,px1,py1)
    write(1001,*)t,U(1,px2,py2),U(2,px2,py2),U(3,px2,py2),U(4,px2,py2)!,U(5,px2,py2),U(6,px2,py2)
    write(1002,*)t,U(1,px3,py3),U(2,px3,py3),U(3,px3,py3),U(4,px3,py3)!,U(5,px3,py3),U(6,px3,py3)
  ct = 0
 endif
        
  !!.........................................................................
  cf=cf+1
  cg=cg+1
  ct=ct+1
  !!.........................................................................
!End of the temporal loop
end do 

  close(unit=1000)
  close(unit=1001)
  close(unit=1002)

!.........................................................................
! Free memory.

deallocate ( dUdm      )
deallocate ( dUdn      )
deallocate ( U         )
deallocate ( uba       )
deallocate ( Tb        )
deallocate ( rhob      )
deallocate ( drhobdn   )
deallocate ( Ub        )
deallocate ( s_acou    )
deallocate (   A       )
deallocate (   B       )
deallocate (   C       )
deallocate (  MI       )
deallocate (  K1       )
deallocate (  DTEMPDM  )
deallocate (  DTEMPDN  )
deallocate (  D2TEMPDM )
deallocate (  D2TEMPDN )
deallocate (  DTAUDM   )
deallocate (  DTAUDN   )
deallocate ( K1        )

    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Euler_MPI:'
    write ( *, '(a)' ) '  Normal end of execution.'
    write ( *, '(a)' ) ' '

end program EULER

