module outt

!***************************************************
!This module is use to create the dates to make 
!the differents plots using GNUPLOT.
!These subrotines are usning to create a matrix 
!that join the different process in a completed 
!matrix in the master process.   

use global
use diffall

!***************************************************

contains
 
subroutine output(U,icount,f,nk)
!**********************************************************************************
!Principal subroutine in this module, is used to call 
!other subroutines. 

! Input 
!       matrix of each process: U_procs,
!       timer: icount 

! Other information is passed by the global.

! Output 
!       completed matrix: U

! Created: Jhonatan Aguirre 
! Date   : 08-04-2017 
! Update : 08-04-2017 
!**********************************************************************************

implicit none
integer::i,j,icount,f,nk
real(kind=ip),dimension(1:nk,1:imax,1:jmax):: U
real(kind=ip),dimension(1:nk,1:imax,1:jmax):: dUcdm,dUcdn
real(kind=ip),dimension(1:imax):: m_all
real(kind=ip),dimension(1:jmax):: n_all
character(40) :: contornofile

integer::ierr,tama
integer::gsize

!subroutine to form a complete matrix, using the 
!submatrixs of the differente process 
     call deropm(dUcdm,U,1,imax,1,jmax,size(U,1)) 
     call deropn(dUcdn,U,1,imax,1,jmax,size(U,1)) 
     write (contornofile,'("contorno_"I0".dat" )' )icount 
     open(unit=icount+201,file=contornofile,status='unknown')
        do i=1,imax
          do j=1,jmax
             write(201+icount,*)m(i), n(j), U(1,i,j), U(2,i,j), U(3,i,j), U(4,i,j), dUcdm(3,i,j)-dUcdn(2,i,j), U(5,i,j)
          end do
          write(201+icount,*)
       end do
      close(unit=icount+201)

end subroutine


subroutine output_energy(U,icount,time,f,tec,nk)
integer::i,j,f,icount,nk
real(kind=ip),dimension(1:nk,1:imax,1:jmax):: U
real(kind=ip),dimension(1:imax)       :: ecx
real(kind=ip)                         :: ekka,ekkb,ekk
real(kind=ip)                         :: time,tec 

                  do i=1,10
                         call kineticenergy_x(U(2,ix1+i-1,:),U(3,ix1+i-1,:),ecx(i))
                  end do
                 
                  write(f+100000,*)time,(ecx(i),i=1,10)
                 
                  do i=1,10
                         call kineticenergy_x(U(2,ix2+i-1,:),U(3,ix2+i-1,:),ecx(i))
                  end do
                 
                  write(f+200000,*)time,(ecx(i),i=1,10)
                 
                  do i=1,10
                         call kineticenergy_x(U(2,ix3+i-1,:),U(3,ix3+i-1,:),ecx(i))
                  end do
                 
                  write(f+300000,*)time,(ecx(i),i=1,10)
end subroutine

subroutine kineticenergy_x(u,v,ec)

integer::i,j
real(kind=ip),dimension(1:jmax):: u,v
real(kind=ip)                         :: ec,decdm
real(kind=ip)                         :: ekka,ekkb,ekk
character(40)                         :: filename

   ec  = 0.d0

   do j=D+1,jmax-D-1
      ekk = 0.5d0*(u(j)**2.d0   + v(j)**2.d0)
      ec  = ec  + ekk
   end do

   ekka=0.5d0*(u(D)**2.d0+v(D)**2.d0) 
   ekkb=0.5d0*(u(jmax-D)**2.d0+v(jmax-D)**2.d0)
   ec =(2.d0*ec +(ekka+ekkb))*dn/2.d0

end subroutine

end module
