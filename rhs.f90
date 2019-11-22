module rhs
use global
use derivs
use eqeuler
use mpi 
USE MPIOWN
USE EMPIRICAL_COEFFICIENTS
contains   

subroutine rhs_euler(dUdt,U,nk)

   implicit none
     
   integer::k,i,j,nk
   real(kind=ip),dimension(1:nk,1:imax, 1:jmax)::dUdt,U
   real(kind=ip),dimension(1,1:imax, 1:jmax)::T_inst

!...Escomento instantaneo,varia com o tempo 

    TEMP(1,1:imax, 1:jmax) = (gamma)*(U(4,1:imax, 1:jmax))/(U(1,1:imax, 1:jmax)) 

    T_inst(:,:,:) = Temp(:,:,:)

    CALL MOLECULAR_MI(T_inst)
    CALL K_COEFFICIENT(T_inst)


!.........Derivadas das variaveis usadas no rhs da equa��o Euler,
!.........O Vetor de derivadas foi passado usando o modulo global
!.........location: derivs.f90 .......

    CALL deriv(U,nk) 
!...Create the tau varibles and make the exchanged with other process
!...Location:eqeuler.f90
!...Use the variable tau that it send by the global.

    CALL TAU_IJ()

!... RHS(right hand side) Equações de Euler with and without PML.
!... Using the variables U, necessary to calculate the rhs of euler equations
!... dUdt have the minus because this is the rhs of the euler equation dUdt=-deuler.
!... location: eqeuler.f90 .......

    call eulernonlinear(dUdt,U,size(U,1))

end subroutine



end module

