module derivs

use global
use diff
USE MPIOWN

contains

  subroutine deriv(U,nk)

    implicit none
    integer :: i,j,nk
    real(kind=ip),dimension(1:nk,1:imax,1:jmax):: U 

  !!!! Derivada otimizada em x, o dudm � passado pelo modulo global  
  !!!! Location: diff.f90
    call deropm(dUdm,U,1,imax,1,jmax,size(u,1)) 
  !!!! Derivada otimizada em y, o dudn � passado pelo modulo global
  !!!! Location: diff.f90
    call deropn(dUdn,U,1,imax,1,jmax,size(u,1)) 
  !CONSTRUCTING \NABLA(K \NABLA (T))
    !To make the comunication of the gosth point between the process 
    !location:mpiown.f90
    !The temperature is calculated in rhs.f90
    call deropm(DTEMPDM,TEMP,1,imax,1,jmax,size(temp,1))
    call deropn(DTEMPDN,TEMP,1,imax,1,jmax,size(temp,1))
    !!Applying metrics to make k.dTdm and k.dTdn
    call kdT()
  !!!!Communication to get the ghost points
  !!!!Location:mpiown.f90
    !!Derivates dtempdm and dtempdn again to make: nabla(k nabla T) 
    call deropm(D2TEMPDM,DTEMPDM,1,imax,1,jmax,size(temp,1))
    call deropn(D2TEMPDN,DTEMPDN,1,imax,1,jmax,size(temp,1))

    do j = 1,jmax
      do i = 1,imax
        psiDiff(1,i,j) = U(1,i,j)*alpha_cond*(dudm(5,i,j)+dudn(5,i,j))
      end do
    end do

    !!Derivates dtempdm and dtempdn again to make: nabla(k nabla T) 
    call deropm(D2psiFDM,psiDiff,1,imax,1,jmax,size(dpsiFdm,1))
    call deropn(D2psiFDN,PsiDiff,1,imax,1,jmax,size(dpsiFdn,1))
    
  end subroutine

  !Inserts transformations metrics to temperature derivatives 
  subroutine kdT() 
    integer :: i , j
    do j = 1,jmax
      do i = 1,imax
        DTEMPDN(1,i,j) = DTEMPDN(1,i,j)*k1(i,j)
        DTEMPDM(1,i,j) = DTEMPDM(1,i,j)*k1(i,j)
      end do
    end do
  end subroutine

end module
