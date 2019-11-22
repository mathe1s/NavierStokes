module diffall
use global
contains
!.................................................................
subroutine deropm(dU,U,ii,fi,ij,fj,nk)

!.......................................................................
!     subrotina para calculo da derivada primeira de uma funcao u
!     utilizando diferencas centradas de quarta ordem nos pontos
!     internos do dominio, diferencas centradas de segunda ordem
!     nos pontos vizinhos 'a fronteira e diferencas unilateral de
!     segunda ordem nos pontos de fronteira.
!.......................................................................
      implicit none
      
      integer:: k, i, j, ii, fi, ij, fj,nk
      
      real(kind=ip),dimension(nk,1:imax,1:jmax):: U,dU

      do k=1,nk
!....   ...................................................................
         do j=ij,fj
!....   ...................................................................
           do i=ii+3,fi-3      
               du(k,i,j)  = dudx(u,k,i,j,nk)
           end do
           i=ii+2
           du(k,i,j) = dx4(u,k,i,j,1,7,1,nk)
           i=fi-2
           du(k,i,j) = dx4(u,k,i,j,7,1,-1,nk)
!....     ...................................................................
           i=ii+1
           du(k,i,j) = dx5(u,k,i,j,1,7,1,nk)
           i=fi-1
           du(k,i,j) = dx5(u,k,i,j,7,1,-1,nk)
!....     ...................................................................
           i=ii
           du(k,i,j) = dx6(u,k,i,j,1,7,1,nk)
           i=fi
           du(k,i,j) = dx6(u,k,i,j,7,1,-1,nk)
!....   ...................................................................
         end do 
!....   ...................................................................
     end do 
end subroutine

real(kind=ip)function dudx(u,k,i,j,nk)

    integer::k,i,j,z,nk
    real(kind=ip),dimension(nk,1:imax,1:jmax):: U
    real(kind=ip),dimension(3)::a

    !a(1)=!0.766855408972008323 ! 0.770882380518d0    
    !a(2)=!-0.163484327177606692!-0.166705904415d0   
    !a(3)=!0.02003774846106832  ! 0.0208431427703d0   

    a(1)= 0.770882380518d0   
    a(2)=-0.166705904415d0   
    a(3)= 0.0208431427703d0 

    dudx=0.d0

    do z=1,3
         dudx= dudx + (a(z)*(u(k,i+z,j)-u(k,i-z,j)))/dm
    end do

end function dudx

real(kind=ip)function dx4(u,k,i,j,iz,fz,p,nk)

    integer::k,i,j,z,iz,fz,p,nk
    real(kind=ip),dimension(nk,1:imax,1:jmax):: U
    real(kind=ip),dimension(7)::a

    if(p==1)then 

       a(1)= 0.049041958000d0   
       a(2)=-0.46884035700d0  
       a(3)=-0.47476091400d0   
       a(4)= 1.27327473700d0  
       a(5)=-0.51848452600d0
       a(6)= 0.16613853300d0
       a(7)=-0.026369431000d0  

       dx4=0.d0
        
       do z=iz,fz,p
          dx4= dx4+a(z)*u(k,i+z-3,j)/dm
       end do 

    else

       a(7)=-0.049041958000d0 
       a(6)= 0.46884035700d0  
       a(5)= 0.47476091400d0  
       a(4)=-1.27327473700d0  
       a(3)= 0.51848452600d0
       a(2)=-0.16613853300d0
       a(1)= 0.026369431000d0
        
       dx4=0.d0
   
       do z=iz,fz,p
          dx4= dx4+a(z)*u(k,i+z-5,j)/dm
       end do 

    end if

end function dx4

real(kind=ip)function dx5(u,k,i,j,iz,fz,p,nk)

    integer::k,i,j,z,iz,fz,p,nk
    real(kind=ip),dimension(nk,1:imax,1:jmax):: U
    real(kind=ip),dimension(7)::a

    if(p==1)then 

       a(1)=-0.20933762200d0    
       a(2)=-1.08487567600d0   
       a(3)= 2.14777605000d0   
       a(4)=-1.38892832200d0   
       a(5)= 0.76894976600d0 
       a(6)=-0.28181465000d0 
       a(7)= 0.048230454000d0 

       dx5=0.d0
        
       do z=iz,fz,p
          dx5= dx5+a(z)*u(k,i+z-2,j)/dm
       end do 

    else

       a(7)= 0.20933762200d0    
       a(6)= 1.08487567600d0    
       a(5)=-2.14777605000d0    
       a(4)= 1.38892832200d0    
       a(3)=-0.76894976600d0  
       a(2)= 0.28181465000d0  
       a(1)=-0.048230454000d0 
        
       dx5=0.d0
   
       do z=iz,fz,p
          dx5= dx5+a(z)*u(k,i+z-6,j)/dm
       end do 

    end if

end function dx5

real(kind=ip)function dx6(u,k,i,j,iz,fz,p,nk)

    integer::k,i,j,z,iz,fz,p
    real(kind=ip),dimension(nk,1:imax,1:jmax):: U
    real(kind=ip),dimension(7)::a

    if(p==1)then 

       a(1)=-2.19228033900d0    
       a(2)= 4.74861140100d0   
       a(3)=-5.10885191500d0   
       a(4)= 4.46156710400d0   
       a(5)=-2.83349874100d0 
       a(6)= 1.12832886100d0 
       a(7)=-0.20387637100d0 

       dx6=0.d0
        
       do z=iz,fz,p
          dx6= dx6+a(z)*u(k,i+z-1,j)/dm
       end do 

    else

       a(7)=  2.19228033900d0    
       a(6)= -4.74861140100d0   
       a(5)=  5.10885191500d0   
       a(4)= -4.46156710400d0   
       a(3)=  2.83349874100d0 
       a(2)= -1.12832886100d0 
       a(1)=  0.20387637100d0 
        
       dx6=0.d0
   
       do z=iz,fz,p
          dx6= dx6+a(z)*u(k,i+z-7,j)/dm
       end do 

    end if

end function dx6

subroutine deropn(dU,U,ii,fi,ij,fj,nk)

!.......................................................................
!
!.......................................................................
      implicit none
      
      integer:: k, i, j, ii, fi, ij, fj,nk
      
      real(kind=ip),dimension(nk,1:imax,1:jmax):: U,dU

      do k=1,nk
!....   ...................................................................
         do i=ii,fi
!....   ...................................................................
           do j=ij+3,fj-3      
               du(k,i,j)  = dudy(u,k,i,j,nk)
           end do
           j=ij+2
           du(k,i,j) = dy4(u,k,i,j,1,7,1,nk)
           j=fj-2
           du(k,i,j) = dy4(u,k,i,j,7,1,-1,nk)
!....     ...................................................................
           j=ij+1
           du(k,i,j) = dy5(u,k,i,j,1,7,1,nk)
           j=fj-1
           du(k,i,j) = dy5(u,k,i,j,7,1,-1,nk)
!....     ...................................................................
           j=ij
           du(k,i,j) = dy6(u,k,i,j,1,7,1,nk)
           j=fj
           du(k,i,j) = dy6(u,k,i,j,7,1,-1,nk)
!....   ...................................................................
         end do 
!....   ...................................................................
     end do 

end subroutine

real(kind=ip)function dudy(u,k,i,j,nk)

    integer::k,i,j,z,nk
    real(kind=ip),dimension(nk,1:imax,1:jmax):: U
    real(kind=ip),dimension(3)::a

    !a(1)=!0.766855408972008323 ! 0.770882380518d0    
    !a(2)=!-0.163484327177606692!-0.166705904415d0   
    !a(3)=!0.02003774846106832  ! 0.0208431427703d0   

    a(1)= 0.770882380518d0   
    a(2)=-0.166705904415d0   
    a(3)= 0.0208431427703d0 

    dudy=0.d0

    do z=1,3
         dudy= dudy + (a(z)*(u(k,i,j+z)-u(k,i,j-z)))/dn
    end do

end function dudy

real(kind=ip)function dy4(u,k,i,j,iz,fz,p,nk)

    integer::k,i,j,z,iz,fz,p
    real(kind=ip),dimension(nk,1:imax,1:jmax):: U
    real(kind=ip),dimension(7)::a

    if(p==1)then 

       a(1)= 0.049041958000d0   
       a(2)=-0.46884035700d0  
       a(3)=-0.47476091400d0   
       a(4)= 1.27327473700d0  
       a(5)=-0.51848452600d0
       a(6)= 0.16613853300d0
       a(7)=-0.026369431000d0  

       dy4=0.d0
        
       do z=iz,fz,p
          dy4= dy4+a(z)*u(k,i,j+z-3)/dn
       end do 

    else

       a(7)=-0.049041958000d0 
       a(6)= 0.46884035700d0  
       a(5)= 0.47476091400d0  
       a(4)=-1.27327473700d0  
       a(3)= 0.51848452600d0
       a(2)=-0.16613853300d0
       a(1)= 0.026369431000d0
        
       dy4=0.d0
   
       do z=iz,fz,p
          dy4= dy4+a(z)*u(k,i,j+z-5)/dn
       end do 

    end if

end function dy4

real(kind=ip)function dy5(u,k,i,j,iz,fz,p,nk)

    integer::k,i,j,z,iz,fz,p,nk
    real(kind=ip),dimension(nk,1:imax,1:jmax):: U
    real(kind=ip),dimension(7)::a

    if(p==1)then 

       a(1)=-0.20933762200d0    
       a(2)=-1.08487567600d0   
       a(3)= 2.14777605000d0   
       a(4)=-1.38892832200d0   
       a(5)= 0.76894976600d0 
       a(6)=-0.28181465000d0 
       a(7)= 0.048230454000d0 

       dy5=0.d0
        
       do z=iz,fz,p
          dy5= dy5+a(z)*u(k,i,j+z-2)/dn
       end do 

    else

       a(7)= 0.20933762200d0    
       a(6)= 1.08487567600d0    
       a(5)=-2.14777605000d0    
       a(4)= 1.38892832200d0    
       a(3)=-0.76894976600d0  
       a(2)= 0.28181465000d0  
       a(1)=-0.048230454000d0 
        
       dy5=0.d0
   
       do z=iz,fz,p
          dy5= dy5+a(z)*u(k,i,j+z-6)/dn
       end do 

    end if

end function dy5

real(kind=ip) function dy6(u,k,i,j,iz,fz,p,nk)

    integer::k,i,j,z,iz,fz,p,nk
    real(kind=ip),dimension(nk,1:imax,1:jmax):: U
    real(kind=ip),dimension(7)::a

    if(p==1)then 

       a(1)=-2.19228033900d0    
       a(2)= 4.74861140100d0   
       a(3)=-5.10885191500d0   
       a(4)= 4.46156710400d0   
       a(5)=-2.83349874100d0 
       a(6)= 1.12832886100d0 
       a(7)=-0.20387637100d0 

       dy6=0.d0
        
       do z=iz,fz,p
          dy6= dy6+a(z)*u(k,i,j+z-1)/dn
       end do 

    else

       a(7)=  2.19228033900d0    
       a(6)= -4.74861140100d0   
       a(5)=  5.10885191500d0   
       a(4)= -4.46156710400d0   
       a(3)=  2.83349874100d0 
       a(2)= -1.12832886100d0 
       a(1)=  0.20387637100d0 
        
       dy6=0.d0
   
       do z=iz,fz,p
          dy6= dy6+a(z)*u(k,i,j+z-7)/dn
       end do 

    end if

end function dy6

end module
