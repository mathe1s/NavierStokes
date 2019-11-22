module diff
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
      integer  ::   k, i, j, ii, fi, ij, fj, nk
      real(kind=ip),dimension(nk,1:imax,1:jmax)   ::   U,dU

      do k=1,nk
         do j=1,jmax
            i=1
            du(k,i,j) = dx6(u,k,i,j,nk)
            i=2
            du(k,i,j) = dx5(u,k,i,j,nk)
            i=3
            du(k,i,j) = dx4(u,k,i,j,nk)
            do i=4,imax-3
               du(k,i,j) = dudx(u,k,i,j,nk)
            end do
            i=imax-2
            du(k,i,j) = dx4I(u,k,i,j,nk)
            i=imax-1
            du(k,i,j) = dx5I(u,k,i,j,nk)
            i=imax
            du(k,i,j) = dx6I(u,k,i,j,nk)
         end do 
      end do  
      call metricsX(du,size(du,1))
   end subroutine

   real(kind=ip) function dudx(u,k,i,j,nk)
      integer :: k,i,j,nk
      real(kind=ip),dimension(nk,1:imax,1:jmax) :: U
      dudx=0.0d0
      dudx = (a_dudx(1)*(u(k,i+1,j)-u(k,i-1,j)))/dm & 
            +(a_dudx(2)*(u(k,i+2,j)-u(k,i-2,j)))/dm &
            +(a_dudx(3)*(u(k,i+3,j)-u(k,i-3,j)))/dm &
   end function dudx

   real(kind=ip) function dx4(u,k,i,j,nk)
      integer  ::  k,i,j,z,nk
      real(kind=ip),dimension(nk,1:imax,1:jmax)  ::   U
      dx4 = 0.d0        
      dx4 = (a_dudx4(1)*u(k,i+1-3,j))/dm &
         + (a_dudx4(2)*u(k,i+2-3,j))/dm &
         + (a_dudx4(3)*u(k,i+3-3,j))/dm &
         + (a_dudx4(4)*u(k,i+4-3,j))/dm &
         + (a_dudx4(5)*u(k,i+5-3,j))/dm &
         + (a_dudx4(6)*u(k,i+6-3,j))/dm &
         + (a_dudx4(7)*u(k,i+7-3,j))/dm &
   end function dx4

   real(kind=ip) function dx4I(u,k,i,j,nk)
      integer  ::  k,i,j,z,nk
      real(kind=ip),dimension(nk,1:imax,1:jmax)  ::   U
      dx4I = 0.d0        
      dx4I = (a_dudx4b(7)*u(k,i+7-5,j))/dm &
         + (a_dudx4b(6)*u(k,i+6-5,j))/dm &
         + (a_dudx4b(5)*u(k,i+5-5,j))/dm &
         + (a_dudx4b(4)*u(k,i+4-5,j))/dm &
         + (a_dudx4b(3)*u(k,i+3-5,j))/dm &
         + (a_dudx4b(2)*u(k,i+2-5,j))/dm &
         + (a_dudx4b(1)*u(k,i+1-5,j))/dm &
   end function dx4I

   real(kind=ip) function dx5(u,k,i,j,nk)
      integer  ::  k,i,j,nk
      real(kind=ip),dimension(nk,1:imax,1:jmax)  ::   U
      dx5 = 0.d0        
      dx5 = (a_dudx5(1)*u(k,i+1-2,j))/dm &
         + (a_dudx5(2)*u(k,i+2-2,j))/dm &
         + (a_dudx5(3)*u(k,i+3-2,j))/dm &
         + (a_dudx5(4)*u(k,i+4-2,j))/dm &
         + (a_dudx5(5)*u(k,i+5-2,j))/dm &
         + (a_dudx5(6)*u(k,i+6-2,j))/dm &
         + (a_dudx5(7)*u(k,i+7-2,j))/dm &
   end function dx5

   real(kind=ip) function dx5I(u,k,i,j,nk)
      integer  ::  k,i,j,nk
      real(kind=ip),dimension(nk,1:imax,1:jmax)  ::   U
      dx5I = 0.d0        
      dx5I = (a_dudx5b(7)*u(k,i+7-6,j))/dm &
         + (a_dudx5b(6)*u(k,i+6-6,j))/dm &
         + (a_dudx5b(5)*u(k,i+5-6,j))/dm &
         + (a_dudx5b(4)*u(k,i+4-6,j))/dm &
         + (a_dudx5b(3)*u(k,i+3-6,j))/dm &
         + (a_dudx5b(2)*u(k,i+2-6,j))/dm &
         + (a_dudx5b(1)*u(k,i+1-6,j))/dm &
   end function dx5I

   real(kind=ip) function dx6(u,k,i,j,nk)
      integer  ::  k,i,j,nk
      real(kind=ip),dimension(nk,1:imax,1:jmax)  ::   U
      dx6 = 0.d0        
      dx6 = (a_dudx6(1)*u(k,i+1-1,j))/dm &
         + (a_dudx6(2)*u(k,i+2-1,j))/dm &
         + (a_dudx6(3)*u(k,i+3-1,j))/dm &
         + (a_dudx6(4)*u(k,i+4-1,j))/dm &
         + (a_dudx6(5)*u(k,i+5-1,j))/dm &
         + (a_dudx6(6)*u(k,i+6-1,j))/dm &
         + (a_dudx6(7)*u(k,i+7-1,j))/dm &
   end function dx6

   real(kind=ip) function dx6I(u,k,i,j,nk)
      integer  ::  k,i,j,nk
      real(kind=ip),dimension(nk,1:imax,1:jmax)  ::   U
      dx6I = 0.d0        
      dx6I = (a_dudx6b(7)*u(k,i+7-7,j))/dm &
         + (a_dudx6b(6)*u(k,i+6-7,j))/dm &
         + (a_dudx6b(5)*u(k,i+5-7,j))/dm &
         + (a_dudx6b(4)*u(k,i+4-7,j))/dm &
         + (a_dudx6b(3)*u(k,i+3-7,j))/dm &
         + (a_dudx6b(2)*u(k,i+2-7,j))/dm &
         + (a_dudx6b(1)*u(k,i+1-7,j))/dm &
   end function dx6I

   subroutine deropn(dU,U,ii,fi,ij,fj, nk)
      implicit none
      integer  ::   k, i, j, ii, fi, ij, fj, nk
      real(kind=ip),dimension(nk,1:imax,1:jmax) :: dU, U
      do k=1,nk
         do i=ii,fi
            j=1
            du(k,i,j) = dy6(u,k,i,j,nk)
            j=2
            du(k,i,j) = dy5(u,k,i,j,nk)
            j=3
            du(k,i,j) = dy4(u,k,i,j,nk)
            do j=4,jmax-3 
               du(k,i,j)  = dudy(u,k,i,j,nk)
            end do
            j=jmax-2
            du(k,i,j) = dy4I(u,k,i,j,nk)
            j=jmax-1
            du(k,i,j) = dy5I(u,k,i,j,nk)
            j=jmax
            du(k,i,j) = dy6I(u,k,i,j,nk)
         end do 
      end do 
      call metricsY(du,size(du,1))
   end subroutine

   real(kind=ip) function dudy(u,k,i,j,nk)
      integer :: k,i,j,z,nk
      real(kind=ip),dimension(nk,1:imax,1:jmax) :: U
      dudy = 0.d0
      dudy =  (a_dudx(1)*(u(k,i,j+1)-u(k,i,j-1)))/dn &
            + (a_dudx(2)*(u(k,i,j+2)-u(k,i,j-2)))/dn &
            + (a_dudx(3)*(u(k,i,j+3)-u(k,i,j-3)))/dn
   end function dudy

   real(kind=ip) function dy4(u,k,i,j,nk)
      integer  ::  k,i,j,z,nk
      real(kind=ip),dimension(nk,1:imax,1:jmax) :: U
      dy4 = 0.d0
      dy4 = a_dudx4(1)*u(k,i,j+1-3)/dn &
          + a_dudx4(2)*u(k,i,j+2-3)/dn &
          + a_dudx4(3)*u(k,i,j+3-3)/dn &
          + a_dudx4(4)*u(k,i,j+4-3)/dn &
          + a_dudx4(5)*u(k,i,j+5-3)/dn &
          + a_dudx4(6)*u(k,i,j+6-3)/dn &
          + a_dudx4(7)*u(k,i,j+7-3)/dn &
   end function dy4

   real(kind=ip) function dy4I(u,k,i,j,nk)
      integer  ::  k,i,j,z,nk
      real(kind=ip),dimension(nk,1:imax,1:jmax) :: U
      dy4I = 0.d0
      dy4I = a_dudx4b(7)*u(k,i,j+7-5)/dn &
           + a_dudx4b(6)*u(k,i,j+6-5)/dn &
           + a_dudx4b(5)*u(k,i,j+5-5)/dn &
           + a_dudx4b(4)*u(k,i,j+4-5)/dn &
           + a_dudx4b(3)*u(k,i,j+3-5)/dn &
           + a_dudx4b(2)*u(k,i,j+2-5)/dn &
           + a_dudx4b(1)*u(k,i,j+1-5)/dn &
   end function dy4I

   real(kind=ip) function dy5(u,k,i,j,nk)
      integer  ::  k,i,j,nk
      real(kind=ip),dimension(nk,1:imax,1:jmax) :: U
      dy5 = 0.d0
      dy5 = + a_dudx5(1)*u(k,i,j+1-2)/dn &
            + a_dudx5(2)*u(k,i,j+2-2)/dn &
            + a_dudx5(3)*u(k,i,j+3-2)/dn &
            + a_dudx5(4)*u(k,i,j+4-2)/dn &
            + a_dudx5(5)*u(k,i,j+5-2)/dn &
            + a_dudx5(6)*u(k,i,j+6-2)/dn &
            + a_dudx5(7)*u(k,i,j+7-2)/dn &
   end function dy5

   real(kind=ip) function dy5I(u,k,i,j,nk)
      integer  ::  k,i,j,z,nk
      real(kind=ip),dimension(nk,1:imax,1:jmax) :: U
      dy5I = 0.d0
      dy5I = a_dudx5b(7)*u(k,i,j+7-6)/dn &
           + a_dudx5b(6)*u(k,i,j+6-6)/dn &
           + a_dudx5b(5)*u(k,i,j+5-6)/dn &
           + a_dudx5b(4)*u(k,i,j+4-6)/dn &
           + a_dudx5b(3)*u(k,i,j+3-6)/dn &
           + a_dudx5b(2)*u(k,i,j+2-6)/dn &
           + a_dudx5b(1)*u(k,i,j+1-6)/dn &
   end function dy5I

   real(kind=ip) function dy6(u,k,i,j,nk)
      integer  ::  k,i,j,nk
      real(kind=ip),dimension(nk,1:imax,1:jmax) :: U
      dy6 = 0.d0
      dy6 = + a_dudx6(1)*u(k,i,j+1-1)/dn &
            + a_dudx6(2)*u(k,i,j+2-1)/dn &
            + a_dudx6(3)*u(k,i,j+3-1)/dn &
            + a_dudx6(4)*u(k,i,j+4-1)/dn &
            + a_dudx6(5)*u(k,i,j+5-1)/dn &
            + a_dudx6(6)*u(k,i,j+6-1)/dn &
            + a_dudx6(7)*u(k,i,j+7-1)/dn &
   end function dy6

   real(kind=ip) function dy6I(u,k,i,j,nk)
      integer  ::  k,i,j,z,nk
      real(kind=ip),dimension(nk,1:imax,1:jmax) :: U
      dy6I = 0.d0
      dy6I = a_dudx6b(7)*u(k,i,j+7-7)/dn &
           + a_dudx6b(6)*u(k,i,j+6-7)/dn &
           + a_dudx6b(5)*u(k,i,j+5-7)/dn &
           + a_dudx6b(4)*u(k,i,j+4-7)/dn &
           + a_dudx6b(3)*u(k,i,j+3-7)/dn &
           + a_dudx6b(2)*u(k,i,j+2-7)/dn &
           + a_dudx6b(1)*u(k,i,j+1-7)/dn &
   end function dy6I

!..........................................................................

   subroutine dermb(dudm,u)
         implicit none
         integer  ::   i, j
         real(kind=ip)  ::   dv4th, dv2nd, onesp, onesm
         real(kind=ip)  ::   ap2, ap1, am1, am2, del, a 
         real(kind=ip),dimension(1:imax)  ::   u, dudm
   !......................................................function statement
         dv4th(ap2,ap1,am1,am2,del) =                                      &
      &          (-ap2 + 8.d0*ap1 - 8.d0*am1 + am2)     / 12.d0 / del
         dv2nd(ap1,am1,del)   = ( ap1 - am1 )             * 0.5d0 / del
         onesp(a,ap1,ap2,del) = (-3.d0*a + 4.d0*ap1 - ap2) * 0.5d0 / del
         onesm(a,am1,am2,del) = ( 3.d0*a - 4.d0*am1 + am2) * 0.5d0 / del
   !.................. Dentro do Dominio ..................................
         do i=3,imax-2
            dudm(i) = dv4th(u(i+2),u(i+1),u(i-1),u(i-2),dm)
         enddo
   !..........................Linha Especiais .............................
         i=2
            dudm(i) = dv2nd(u(i+1),u(i-1),dm)
   !.......................................................................
         i=imax-1
            dudm(i) = dv2nd(u(i+1),u(i-1),dm)

   !.......................................................................
         i=imax
            dudm(i) = onesm(u(i),u(i-1),u(i-2),dm)

   !.......................................................................
         i=1
            dudm(i) = onesp(u(i),u(i+1),u(i+2),dm)
         
   end subroutine

   !Inserts transformations metrics
   subroutine metricsY(dU,nk)
   real(kind=ip),dimension(nk,1:imax,1:jmax),intent(inout) :: dU
   integer :: i , j, k, nk
   do k=1,nk
      do j = 1,jmax
         do i = 1,imax
            !***********************ylayer_left***********************************
            if((i>=D+1).and.(i<=imax-D).and.(j>=1).and.(j<=D))then
               du(k,i,j) = du(k,i,j)*meshy(j)
            end if
            !***********************ylayer_right***********************************
            if((i>=D+1).and.(i<=imax-D).and.(j>=jmax-D+1).and.(j<=jmax))then
               du(k,i,j) = du(k,i,j)*meshy(j)
            end if
            !***********************corners********************************* 
         !***********************1********************************* 
            if ((i>=1).and.(i<=D).and.(j>=jmax-D+1).and.(j<=jmax))then
               du(k,i,j) = du(k,i,j)*meshy(j)
            end if
      !***********************2********************************* 
            if ((i>=1).and.(i<=D).and.(j>=1).and.(j<=D)) then
               du(k,i,j) = du(k,i,j)*meshy(j)
            end if
      !***********************3********************************* 
            if ((i>=imax-D+1).and.(i<=imax).and.(j>=1).and.(j<=D)) then
               du(k,i,j) = du(k,i,j)*meshy(j)
            end if
      !**********************4********************************* 
            if ((i>=imax-D+1).and.(i<=imax).and.(j>=jmax-D+1).and.(j<=jmax)) then
               du(k,i,j) = du(k,i,j)*meshy(j)
            end if
      !******************************************************** 
         end do
      end do
   end do
   end subroutine


   !Inserts transformations metrics to temperature derivatives 
   subroutine metricsX(dU,nk)
   real(kind=ip),dimension(nk,1:imax,1:jmax),intent(inout)   ::   dU
   integer :: i , j, k, nk
   do k=1,nk
      do j = 1,jmax
         do i = 1,imax
      !***********************xlayer_bottom********************************* 
            if ((i>=1).and.(i<=D).and.(j>=D+1).and.(j<=jmax-D))then
               du(k,i,j) = du(k,i,j)*meshx(i)
            end if
      !***********************xlayer_top********************************* 
            if((i>=imax-D+1).and.(i<=imax).and.(j>=D+1).and.(j<=jmax-D))then
               du(k,i,j) = du(k,i,j)*meshx(i)
            end if
      !***********************corners********************************* 
      !***********************1********************************* 
            if ((i>=1).and.(i<=D).and.(j>=jmax-D+1).and.(j<=jmax))then
               du(k,i,j) = du(k,i,j)*meshx(i)
            end if
      !***********************2********************************* 
            if ((i>=1).and.(i<=D).and.(j>=1).and.(j<=D)) then
               du(k,i,j) = du(k,i,j)*meshx(i)
            end if
      !***********************3********************************* 
            if ((i>=imax-D+1).and.(i<=imax).and.(j>=1).and.(j<=D)) then
               du(k,i,j) = du(k,i,j)*meshx(i)
            end if
      !**********************4********************************* 
            if ((i>=imax-D+1).and.(i<=imax).and.(j>=jmax-D+1).and.(j<=jmax)) then
               du(k,i,j) = du(k,i,j)*meshx(i)
            end if
      !******************************************************** 
         end do
      end do
   end do
   end subroutine

end module
