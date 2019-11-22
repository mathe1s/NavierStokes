module TransformationMetrics
use global
use mpi
!******************************************************************
!Book: Computational Fluid Mechanics And Heat Transfer
!Author: John Tannehill, Dale Anderson, Richard Pletcher
!Pg: 352 (pdf), 337(book)
subroutine Transformation_3()
implicit none
!.......................................................................

real(kind=ip) :: CoefA, CoefB, CoefC
integer       :: YDomain

YDomain = nmax - nmin

!StrechingParameter: Initialized at dados.in
!YcParameter: Initialized at dados.in

CoefA = dexp(StrechingParameter) - 1.0d0
CoefB = dexp(-StrechingParameter) - 1.0d0
CoefC = YcParameter/YDomain

BParameter = 1.0d0/2.0d0/StrechingParameter &
            * log((1.0d0+CoefA*CoefC)/(1.0d0+CoefB*CoefC))

allocate (Ybar(s_s:e_s))
do j=s_s,e_s
    Ybar(j) = BParameter + 1.0d0/StrechingParameter &
              * 1.0d0/( sinh((n(j)/YcParameter -1.0d0)) &
              * sinh(StrechingParameter*BParameter) )
end do

allocate (Yref(s_s:e_s))
do j=s_s,e_s
    Yref(j) = YcParameter*(1.0d0 + &
               sinh(StrechingParameter*(Ybar(j)-BParameter)) / &
               sinh(StrechingParameter*BParameter) )
end do


meshy   = 1.d0

do j=s_s,e_s
    meshy(j) = sinh(StrechingParameter*BParameter) / &
               StrechingParameter*YcParameter / &
               dsqrt( &
               1.0d0 + ((Yref(j)/YcParameter) - 1.0d0)**2.0d0 / &
               sinh(StrechingParameter*BParameter)/sinh(StrechingParameter*BParameter) &
               )
end do
end subroutine Transformation_3

end module TransformationMetrics
