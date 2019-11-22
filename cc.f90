module cc
use global
contains
!************Boundary conditions***********

subroutine ccpml(U,nk)

implicit none
integer :: nk
real(kind=ip),dimension(1:nk,1:imax,1:jmax)::U 

!  rho=G(1,:,:)
!  u=G(2,:,:)
!  v=G(3,:,:)
!  P=G(4,:,:)

!Inflow Boundary!!!!!!!
!On perfectly matched layer as an absorbing boundary condition
!Fang Q. Hu
!AIAA Paper 96-1664
!"odd" extrapolation of function is carried out for simplicity

  U(1,1,1:jmax) = rhob(1:jmax)!2.d0*U(1,2,1:jmax)-U(1,3,1:jmax)
  U(2,1,1:jmax) = 2.d0*U(2,2,1:jmax)-U(2,3,1:jmax)
  U(3,1,1:jmax) = 0.d0!2.d0*U(3,2,1:jmax)-U(3,3,1:jmax)
  U(4,1,1:jmax) = 1.d0/gamma!2.d0*U(4,2,1:jmax)-U(4,3,1:jmax) 
  U(5,1,1:jmax) = psib(1:jmax)!2.d0*U(5,2,1:jmax)-U(5,3,1:jmax)

!outflow
! outflow boundary

  U(1,imax,1:jmax)   = 2.d0*U(1,imax-1,1:jmax)-U(1,imax-2,1:jmax)
  U(2,imax,1:jmax)   = 2.d0*U(2,imax-1,1:jmax)-U(2,imax-2,1:jmax)
  U(3,imax,1:jmax)   = 2.d0*U(3,imax-1,1:jmax)-U(3,imax-2,1:jmax)
  U(4,imax,1:jmax)   = 1.d0/gamma!2.d0*U(4,imax-1,1:jmax)-U(4,imax-2,1:jmax) 
  U(5,imax,1:jmax)   = 2.d0*U(5,imax-1,1:jmax)-U(5,imax-2,1:jmax)

! upper boundary 
     U(1,1:imax,jmax)   = 2.d0*U(1,1:imax,jmax-1)-U(1,1:imax,jmax-2)
     U(2,1:imax,jmax)   = 2.d0*U(2,1:imax,jmax-1)-U(2,1:imax,jmax-2) 
     U(3,1:imax,jmax)   = 2.d0*U(3,1:imax,jmax-1)-U(3,1:imax,jmax-2)
     U(4,1:imax,jmax)   = 1.d0/gamma!2.d0*U(4,1:imax,jmax-1)-U(3,:,jmax-2) 
     U(5,1:imax,jmax)   = 1.0d0

! All perturbations all zero in -+inf
!   U(1,:,jmax)   = 0.d0
!   U(2,:,jmax)   = 0.d0 
!   U(3,:,jmax)   = 0.d0
!   U(4,:,jmax)   = 0.d0 
 

! lower boundary

! All perturbations all zero in -+inf
    U(1,1:imax,1)   = 2.d0*U(1,1:imax,1+1)-U(1,1:imax,1+2)
    U(2,1:imax,1)   = 2.d0*U(2,1:imax,1+1)-U(2,1:imax,1+2)
    U(3,1:imax,1)   = 2.d0*U(3,1:imax,1+1)-U(3,1:imax,1+2)
    U(4,1:imax,1)   = 1.d0/gamma!2.d0*U(4,1:imax,1+1)-U(4,:,1+2) 
    U(5,1:imax,1)   = 0.0d0

!    U(1,:,1)     = 0.d0 
!    U(2,:,1)     = 0.d0 
!    U(3,:,1)     = 0.d0 
!    U(4,:,1)     = 0.d0 

    
   
end subroutine     

subroutine ccq(q1,nk)

implicit none
integer :: nk
real(kind=ip),dimension(nk,im,jm)::q1

!Inflow Boundary!!!!!!!

  q1(1,1,:)    = 2.d0*q1(1,2,:)-q1(1,3,:)
  q1(2,1,:)    = 2.d0*q1(2,2,:)-q1(1,3,:)
  q1(3,1,:)    = 2.d0*q1(3,2,:)-q1(1,3,:)
  q1(4,1,:)    = 2.d0*q1(4,2,:)-q1(1,3,:)


! outflow boq2ndary

   q1(1,imax,:)   = 2.d0*q1(1,imax-1,:)-q1(1,imax-2,:)
   q1(2,imax,:)   = 2.d0*q1(2,imax-1,:)-q1(2,imax-2,:)
   q1(3,imax,:)   = 2.d0*q1(3,imax-1,:)-q1(3,imax-2,:)
   q1(4,imax,:)   = 2.d0*q1(4,imax-1,:)-q1(4,imax-2,:) 


   
! upper boq2ndary 
   q1(1,:,jmax)   = 2.d0*q1(1,:,jmax-1)-q1(1,:,jmax-2)
   q1(2,:,jmax)   = 2.d0*q1(2,:,jmax-1)-q1(2,:,jmax-2) 
   q1(3,:,jmax)   = 2.d0*q1(3,:,jmax-1)-q1(3,:,jmax-2)
   q1(4,:,jmax)   = 2.d0*q1(4,:,jmax-1)-q1(3,:,jmax-2) 


! lower boq2ndary

    q1(1,:,1)     = 2.d0*q1(1,:,2)-q1(1,:,3)
    q1(2,:,1)     = 2.d0*q1(2,:,2)-q1(2,:,3)
    q1(3,:,1)     = 2.d0*q1(3,:,2)-q1(3,:,3)
    q1(4,:,1)     = 2.d0*q1(4,:,2)-q1(4,:,3) 
    
   
end subroutine     
end module
