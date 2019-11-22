Module rhsq
use global

contains

subroutine  rhsq1(u,qq,q1dt)

      implicit none
      integer::ii,fi,ij,fj,k,i,j
      real(kind=ip),dimension(4,i_s:d_s,s_s-g_p:e_s+g_p):: u,qq,q1dt

!Loop to calculate the rhsq of q variable for pml to the Euler equations (d/dt=rhsq) with PML 

do k=1,4 

   do i=i_s,d_s 
   
      do j=s_s,e_s
   
        if( (i>=1).and.(i<=imaxpml))then
   
           if((j>=jmaxpml+1).and.(j<=jmax-D))then
              call q1xlayer(U,qq,q1dt,k,i,j) 
           end if
   
           if((j>=1).and.(j<=jmaxpml))then
              call q1corner(U,qq,q1dt,k,i,j)      
           end if
   
           if((j>=jmax-D+1).and.(j<=jmax))then
              call q1corner(U,qq,q1dt,k,i,j)  
           end if
   
        end if
   
        if( (i>=imax-D+1).and.(i<=imax))then
   
           if((j>=jmaxpml+1).and.(j<=jmax-D))then
              call q1xlayer(U,qq,q1dt,k,i,j)
           end if 
   
           if((j>=1).and.(j<=jmaxpml))then
              call q1corner(U,qq,q1dt,k,i,j)    
           end if 
   
           if((j>=jmax-D+1).and.(j<=jmax))then
              call q1corner(U,qq,q1dt,k,i,j)    
           end if 
   
        end if
   
        if( (i>=imaxpml+1).and.(i<=imax-D))then
   
           if((j>=1).and.(j<=jmaxpml))then
              call q1ylayer(U,qq,q1dt,k,i,j)     
           end if 
   
           if((j>=jmax-D+1).and.(j<=jmax))then
              call q1ylayer(U,qq,q1dt,k,i,j) 
           end if 
             
        end if
   
      end do
   
   end do

end do


end subroutine

subroutine  rhsq2(u,qq,q2dt,i,j)

      implicit none
      integer::ii,fi,ij,fj,k,i,j
      real(kind=ip),dimension(4,im,jm):: u,qq,q2dt

     
      !call q2corner(U,qq,q2dt,1,imaxpml,1,jmaxpml)      
      !call q2ylayer(U,qq,q2dt,imaxpml+1,imax-D,1,jmaxpml)     
      !call q2corner(U,qq,q2dt,imax-D+1,imax,1,jmaxpml)    

      !call q2xlayer(U,qq,q2dt,1,imaxpml,jmaxpml+1,jmax-D) 
      !call q2xlayer(U,qq,q2dt,imax-D+1,imax,jmaxpml+1,jmax-D)

      !call q2corner(U,qq,q2dt,1,imaxpml,jmax-D+1,jmax)  
      !call q2ylayer(U,qq,q2dt,imaxpml+1,imax-D,jmax-D+1,jmax) 
      !call q2corner(U,qq,q2dt,imax-D+1,imax,jmax-D+1,jmax)    

end subroutine


subroutine q1xlayer(U,qq,q1dt,k,i,j)   

      implicit none
      integer::ii,fi,ij,fj,k,i,j
      real(kind=ip),dimension(4,i_s:d_s,s_s-g_p:e_s+g_p):: u,qq,q1dt

         q1dt(k,i,j) = U(k,i,j) &
       & - delta*dq1dm(k,i,j)/meshx(i) - sigmax(i)*Beta*delta*qq(k,i,j)

end subroutine


subroutine q1ylayer(U,qq,q1dt,k,i,j)   

      implicit none
      integer::ii,fi,ij,fj,k,i,j
      real(kind=ip),dimension(4,i_s:d_s,s_s-g_p:e_s+g_p):: u,qq,q1dt

      q1dt(k,i,j) = U(k,i,j) &
      & - delta*dq1dn(k,i,j)/meshy(j)  - sigmay(j)*Beta*delta*qq(k,i,j)

end subroutine

subroutine q1corner(U,qq,q1dt,k,i,j)   

      implicit none
      integer::ii,fi,ij,fj,k,i,j
      real(kind=ip),dimension(4,i_s:d_s,s_s-g_p:e_s+g_p):: u,qq,q1dt
      real(kind=ip):: dt2


      q1dt(k,i,j) = U(k,i,j) &
        & - delta*dq1dm(k,i,j)/meshx(i) - sigmax(i)*Beta*delta*qq(k,i,j) &
        & - delta*dq1dn(k,i,j)/meshy(j) - sigmay(j)*Beta*delta*qq(k,i,j) &
        & - sigmax(i)*Beta*delta*qq(k,i,j) - sigmay(j)*Beta*delta*qq(k,i,j)
   
end subroutine

subroutine q2xlayer(U,qq,q2dt,k,i,j)   

      implicit none
      integer::ii,fi,ij,fj,k,i,j
      real(kind=ip),dimension(4,i_s:d_s,s_s-g_p:e_s+g_p):: u,qq,q2dt
      real(kind=ip):: dt2


      q2dt(k,i,j) = dUdn(k,i,j) - dUbdn(k,i,j) - sigmax(i)*qq(k,i,j)   
   
end subroutine

subroutine q2ylayer(U,qq,q2dt,k,i,j)  

      implicit none
      integer::ii,fi,ij,fj,k,i,j
      real(kind=ip),dimension(4,i_s:d_s,s_s-g_p:e_s+g_p):: u,qq,q2dt
      real(kind=ip):: dt2


      q2dt(k,i,j) = dUdn(k,i,j)/meshy(j) - dUbdn(k,i,j)/meshy(j) - qq(k,i,j)*sigmay(j)  
   
end subroutine

subroutine q2corner(U,qq,q2dt,k,i,j)    

      implicit none
      integer::ii,fi,ij,fj,k,i,j
      real(kind=ip),dimension(4,i_s:d_s,s_s-g_p:e_s+g_p):: u,qq,q2dt
      real(kind=ip):: dt2


      q2dt(k,i,j) = dUdn(k,i,j)/meshy(j) - dUbdn(k,i,j)/meshy(j)   - qq(k,i,j)*sigmay(j) -qq(k,i,j)*sigmax(i)
   
end subroutine

end module
