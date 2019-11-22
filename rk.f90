module rk

use global
use filtering
use rhs
use diff
use cc

contains


subroutine rk_euler(U,t,iter,nk)
      implicit none         
      integer:: k, i, j,irk,iter, nk
      real(kind=ip)::t,dt2, dt6
      real(kind=ip),dimension(6)::alphat,betat,ct
      real(kind=ip),dimension(1:nk,1:imax,1:jmax):: U,dudt,w,wq1,wq2

!......................................................................
     !inicializacao das variaveis do avanco temporal segundo Berland
            w   = 0.d0
            wq1 = 0.d0
            wq2 = 0.d0
!......................................................................
     !coefiecientes do metodo lddrk do Berland de 6 passo e 4 ordem

      alphat(1) =  0.d0
      alphat(2) = -0.737101392796d0
      alphat(3) = -1.634740794341d0
      alphat(4) = -0.744739003780d0
      alphat(5) = -1.469897351522d0
      alphat(6) = -2.813971388035d0
      
      betat(1)  =  0.032918605146d0
      betat(2)  =  0.823256998200d0
      betat(3)  =  0.381530948900d0
      betat(4)  =  0.200092213184d0
      betat(5)  =  1.718581042715d0
      betat(6)  =  0.27d0
     
          ct(1) =  0.0d0
          ct(2) =  0.032918605146d0
          ct(3) =  0.249351723343d0
          ct(4) =  0.466911705055d0
          ct(5) =  0.582030414044d0
          ct(6) =  0.847252983783d0

!...........................................................................
!Runge Kutta loop, six steps for non linear operators and dispersion relation
!preserved 
!paper: 
do irk=1,6
!...........................................................................
!Acoustic font, location: here
    call  acousticfont(ct(irk),iter)
!............................................................................
!Temporal derivs, dudt=right hand side, location: rhs.f90
    call rhs_euler(dUdt ,U,nk)
!............................................................................
!   Steps-Runge-Kutta, location: here
    call  rk_steep(alphat(irk),betat(irk),ct(irk),U,dUdt,w,nk)
!............................................................................
!Buffer Zone non reflecting Boundary conditons
    call bufferxy(U,nk)
!............................................................................
!Filter is applied at each estage of the rk to ETA boundary condition 
    call lddfilterx(U,size(U,1))
    call lddfiltery(U,size(U,1))
!.............................................................................
!Boundary conditons
!Boundary conditons
    call ccpml(U,size(U,1))
end do     
!............................................................................
end subroutine

subroutine bufferxy(U,nk)
integer::k,l,i,j,ii,fi,ij,fj,nk
real(kind=ip),dimension(1:nk,1:imax,1:jmax):: U

   do j=1,jmax
       do i=1,imax
    
!***********************ylayer_left***********************************
    
        if((i>=1).and.(i<=imax))then
           if((j>=1).and.(j<=D))then 
               call buffery(U,i,j,nk)
           end if 
           if((j>=jmax-D+1).and.(j<=jmax))then
               call buffery(U,i,j,nk)
          end if 
        end if 
    
    
!***********************xlayer_bottom********************************* 
    
          if ((j>=1).and.(j<=jmax)) then
            if ((i>=1).and.(i<=D))then
                call bufferx(U,i,j,nk)
           end if 
           if((i>=imax-D+1).and.(i<=imax))then
                call bufferx(U,i,j,nk)
           end if 
         end if 

     end do 
   end do
    

end subroutine

subroutine bufferx(U,i,j,nk)
integer::k,l,i,j,ii,fi,ij,fj,nk
real(kind=ip),dimension(1:nk,1:imax,1:jmax):: U

do k=1,nk
!Paper:Simulation techniques for spacialy evolving instabilities 
!      in compressible flows over a flat plate
!Autor:B. Wasistho
!Computers & Fluids Vol. 26. No. 7. pp. 713-739. 1997 
!*********************************************************************
         U(k,i,j)=Ub(k,i,j)+sigmax(i)*(U(k,i,j)-Ub(k,i,j))
         !U(k,i,j)=sigmax(i)*U(k,i,j)
!*********************************************************************
!Paper:
!The evaluation of non_reflecting bopundary conditions for duct acoustic 
!computation.
!Autor: S. K. Ricahards..
!Jounal of Sound and VIbration 270 (2004) 539-557
!      in compressible flows over a flat plate
         !U(k,i,j)=U(k,i,j)-sigmax(i)*(U(k,i,j)-Ub(k,i,j))
         !U(k,i,j)=U(k,i,j)-sigmax(i)*(U(k,i,j))
!*********************************************************************

end do 

end subroutine

subroutine buffery(U,i,j,nk)
integer::k,l,i,j,ii,fi,ij,fj,nk
real(kind=ip),dimension(1:nk,1:imax,1:jmax):: U

do k=1,nk
         U(k,i,j)=Ub(k,i,j)+sigmay(j)*(U(k,i,j)-Ub(k,i,j))
         !U(k,i,j)=sigmay(j)*U(k,i,j)
         !U(k,i,j)=U(k,i,j)-sigmay(j)*(U(k,i,j)-Ub(k,i,j))
         !U(k,i,j)=U(k,i,j)-sigmay(j)*(U(k,i,j))
end do 

end subroutine

subroutine rk_steep(alphat,betat,ct,U,dudt,w,nk)

     implicit none
     integer::k,i,j,nk
     real(kind=ip)::alphat,betat,ct 
     real(kind=ip),dimension(nk,1:imax,1:jmax):: U,dudt,w

      do k=1,nk 
            do j =2,jmax-1
              do i =2,imax-1 
                  w(k,i,j) = alphat*w(k,i,j)+dt*dudt(k,i,j)
                  U(k,i,j) = U(k,i,j) + betat*w(k,i,j)
              end do 
            end do  
      end do 

end subroutine

subroutine acousticfont(ct,iter)
!...........................................................................
!...........................................................................
!       Fonte acustica para gerar o pulso de pressão que instabiliza o escoamento  
!...........................................................................
!     allocate (s_acou(1:imax,s_s:e_s))

     implicit none
     integer::i,j,iter
     real(kind=ip)::ct 

!                        Parametros Fonte
!********Raio da fonte
     !At dados.in r0  = 0.03d0  
!*******A FONTE ESTA EN X,Y(0.0)
     do  j=1,jmax
        do i =1,imax 
           s_acou(i,j)= amplitude*sin(wf*((iter+ct)*dt))*  &
            & exp(-log(2.d0)*((m(i)-m_pt_acous)**2+(n(j)-n_pt_acous)**2)/(r0_acous*r0_acous))
        end do
     end do

end subroutine

subroutine energiacinetica(u,v,ec,decdm)

integer::i,j
real(kind=ip),dimension(im,jm)::u,v
real(kind=ip),dimension(im)::ec,decdm
real(kind=ip)::ekka,ekkb,ekk

do i=D+1,imax-D

   ec(i) = 0.d0
   do j=D+1,jmax-D
      ekk=.5d0*(u(i,j)**2.d0   + v(i,j)**2.d0)
      ec(i) = ec(i) + ekk
   end do

   ekka=0.5d0*(u(i,D+1)**2.d0+v(i,D+1)**2.d0) 
   ekkb=0.5d0*(u(i,jmax-D)**2.d0+v(i,jmax-D)**2.d0)
   ec(i)=(2.d0*ec(i)+(ekka+ekkb))*dn/2.d0

end do

call dermb(decdm,ec)

end subroutine

subroutine energiacineticax(i,u,v,ec)

integer::i,j
real(kind=ip),dimension(im,jm)::u,v
real(kind=ip)::ekka,ekkb,ekk,ec


   ec  = 0.d0
   do j=D+1,jmax-D
      ekk=.5d0*(u(i,j)**2.d0   + v(i,j)**2.d0)
      ec  = ec  + ekk
   end do

   ekka=0.5d0*(u(i,D+1)**2.d0+v(i,D+1)**2.d0) 
   ekkb=0.5d0*(u(i,jmax-D)**2.d0+v(i,jmax-D)**2.d0)
   ec =(2.d0*ec +(ekka+ekkb))*dn/2.d0


end subroutine

end module
