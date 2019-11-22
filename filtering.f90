module filtering


use global
use mpi

contains

subroutine lddfilterx(U,nk)
implicit none
!Paper: A famlily of low dispersive and low dissapative 
!explicit schemes  for computing the aerodynamics noise 
!Autor: Bogey ann Bailly

integer::i,j,k,s,nk
integer::ii,fi,ij,fj,i1,f1,i2,f2
real(kind=ip),dimension(1:nk,1:imax,1:jmax)::U
real(kind=ip),dimension(-6:6)::df
real(kind=ip),dimension(-1:5)::dfb
real(kind=ip)::Du

!right  and left  boundaries 
!coeficients 

dfb(-1) =-0.057717512738d0
dfb(0)  = 0.199278374994d0
dfb(1)  =-0.292668277650d0
dfb(2)  = 0.244537361546d0
dfb(3)  =-0.134605018019d0
dfb(4)  = 0.056184263460d0
dfb(5)  =-0.015009191593d0

ij=1
fj=jmax

  i1=2
  i2=7

  do k=1,nk
    do j=ij,fj
       do i=i1,i2
          Du=0.d0
          do s=-1,5
             Du=Du+dfb(s)*u(k,i+s,j)
          end do
          u(k,i,j) = u(k,i,j)-sigmaf*Du
      end do 
    end do
  end do  

  ii=7
  fi=imax-6

  df(-6)= 0.001254597714d0
  df(-5)=-0.008520738659d0
  df(-4)= 0.029662754736d0
  df(-3)=-0.069975429105d0
  df(-2)= 0.123632891797d0
  df(-1)=-0.171503832236d0
  df(0) = 0.190899511506d0
  df(1) =-0.171503832236d0
  df(2) = 0.123632891797d0
  df(3) =-0.069975429105d0
  df(4) = 0.029662754736d0
  df(5) =-0.008520738659d0
  df(6) = 0.001254597714d0
  Du=0.d0
  
  do k=1,nk
      do i=ii,fi
         do j=ij,fj
          Du=0.d0
            do s=-6,6
               Du=Du+df(s)*u(k,i+s,j)
            end do 
            u(k,i,j) = u(k,i,j)-sigmaf*Du
         end do 
      end do
  end do  

  i1=imax-1
  i2=imax-6

  do k=1,nk
    do j=ij,fj
         do i=i1,i2,-1
          Du=0.d0
          do s=-1,5
             Du=Du+dfb(s)*u(k,i-s,j)
          end do
          u(k,i,j) = u(k,i,j)-sigmaf*Du
      end do 
    end do
  end do  
 
end subroutine

subroutine lddfiltery(U,nk)
implicit none
!Paper: A famlily of low dispersive and low dissapative 
!explicit schemes  for computing the aerodynamics noise 
!Autor: Bogey ann Bailly

integer::i,j,k,s,nk
integer::ii,fi,ij,fj,j1,j2
real(kind=ip),dimension(1:nk,1:imax,1:jmax)::U
real(kind=ip),dimension(-6:6)::df
real(kind=ip),dimension(-1:5)::dfb
real(kind=ip)::Du

!coeficients to bottom and top boundaries

!we need a point to the one side and 6 ponit to another side
dfb(-1)  =-0.057717512738
dfb( 0)  = 0.199278374994
dfb( 1)  =-0.292668277650
dfb( 2)  = 0.244537361546
dfb( 3)  =-0.134605018019
dfb( 4)  = 0.056184263460
dfb( 5)  =-0.015009191593

ii=1
fi=imax

  j1=2
  j2=7

  do k=1,nk
    do i=ii,fi
      do j=j1,j2

          Du=0.d0
          do s=-1,5
             Du=Du+dfb(s)*u(k,i,j+s)
          end do
          u(k,i,j) = u(k,i,j)-sigmaf*Du

      end do 
    end do
  end do  

  ij=7
  fj=jmax-6

  df(-6) = 0.001254597714d0
  df(-5) =-0.008520738659d0
  df(-4) = 0.029662754736d0
  df(-3) =-0.069975429105d0
  df(-2) = 0.123632891797d0
  df(-1) =-0.171503832236d0
  df( 0) = 0.190899511506d0
  df( 1) =-0.171503832236d0
  df( 2) = 0.123632891797d0
  df( 3) =-0.069975429105d0
  df( 4) = 0.029662754736d0
  df( 5) =-0.008520738659d0
  df( 6) = 0.001254597714d0

  Du=0.d0
  
  do k=1,nk
     do i=ii,fi
        do j=ij,fj
           Du=0.d0
           do s=-6,6
              Du=Du+df(s)*u(k,i,j+s)
           end do 
           u(k,i,j) = u(k,i,j)-sigmaf*Du
        end do 
    end do
  end do  

  j1=jmax-1
  j2=jmax-6

  do k=1,nk
    do i=ii,fi
      do j=j1,j2,-1

          Du=0.d0
          do s=-1,5
             Du=Du+dfb(s)*u(k,i,j-s)
          end do
          u(k,i,j) = u(k,i,j)-sigmaf*Du

      end do 
    end do
  end do  

end subroutine

end module
