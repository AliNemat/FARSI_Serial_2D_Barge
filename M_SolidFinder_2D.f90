Module M_SolidFinder_2D
use M_Math,                 only:  HFF
use M_General_2D,           only:  x,y,nx,ny,test
use M_Platform_Constant_2D, only:  r2,Len
use M_Mesh_2D,              only:  hx,hy
implicit none
parameter smooth=1
integer         ,dimension(1:nx-1,1:ny-1)        :: ccL       
double precision,dimension(1:nx-1,1:ny-1)        :: cc,CI,CF,CFS
double precision,dimension(8)                    :: mm,mm1,mm2,alpha
!DOUBLE PRECISION                           :: HFF


contains






subroutine ImportSolidPoint1(npoints,point,tag)
  use M_General_2D,           only: YFreeIni
  use M_Platform_Constant_2D, only: Floatx
  implicit none

  double precision,intent(out),dimension (npoints+1,2)     ::   point
  integer         ,intent(in)                              ::   npoints,tag

  double precision                                         ::   IniH 
  
  if (tag.eq.1) then 
    include 'Par_SolidFinder1_2D.txt'
  else if (tag.eq.2) then
    include 'Par_SolidFinder2_2D.txt'
  end if 

  point(npointS+1,1)=point(1,1)
  point(npointS+1,2)=point(1,2)

  
  
return 

end subroutine


subroutine SolidCellFinder(npoints,point,XBarV,CFSOut)
implicit none 
double precision,intent(in),dimension(npoints+1,2)     :: point
double precision,intent(in)                            :: XBarV(2)
double precision,intent(out),dimension(1:nx-1,1:ny-1)  :: CFSOut
integer         ,intent(in)                            :: npoints

  call BorderFinder(npoints,point,XBarV)
  call insideFinder(npoints,point,XBarV)
  call FractionFinder(npoints,point)
  call FractionSmoother()

  CFSOut(:,:)=CFS(:,:)
return

end subroutine



     

subroutine BorderFinder(npoints,point,XBarV)
implicit none

double precision,intent(in),dimension(npoints+1,2)     :: point
double precision,intent(in)                            :: XBarV(2)
integer         ,intent(in)                            :: npoints
integer kk,i,j


  
  do kk=1,npointS

    if (point(kk+1,1).ne.point(kk,1)) then 

      mm(kk)=(point(kk+1,2)-point(kk,2))/(point(kk+1,1)-point(kk,1))

      mm1(kk)=-mm(kk)/sqrt(1+mm(kk)*mm(kk))
      mm2(kk)=  1.0/sqrt(1+mm(kk)*mm(kk))
       
      alpha(kk)=(point(kk,2)-mm(kk)*point(kk,1))/sqrt(1+mm(kk)*mm(kk))

    else

       mm1(kk)=1
       mm2(kk)=0
       alpha(kk)=point(kk,1)

    end if

  end do 
!print*, "mm1",mm1(:)
!print*, "mm2",mm2(:)
!print*, "alpha",alpha(:)

ccl(:,:)=0
  do i=1, nx-1
    if ( abs (x(i)-XBarV(1)).gt.(2*r2) )  then
      cc(i,:)=0
    else
     
      do j=1,ny-1
        if ( abs (y(j)-XBarV(2)).gt.(len) ) then
          cc(i,j)=0
        else 
          do kk=1,npoints

            if (mm2(kk).eq.0) then
          !    print*, " I am here1 " 
              if ( (x(i)-0.5*hx(i)).lt.alpha(kk).AND.(x(i)+0.5*hx(i)).gt.Alpha(kk) ) then
                if ( ( (x(i)-point(kk,1  ))**2 + (y(j)-point(kk,2  ))**2 )  .lt.( (point(kk,1)-point(kk+1,1))**2 + (point(kk,2)-point(kk+1,2))**2 ).AND. &
                     ( (x(i)-point(kk+1,1))**2 + (y(j)-point(kk+1,2))**2 )  .lt.( (point(kk,1)-point(kk+1,1))**2 + (point(kk,2)-point(kk+1,2))**2 ) ) then

                !     print*, " I am here2 " 
                  cc(i,j)=1.0
                  ccl(i,j)=kk
                  exit
                end if  
              end if 
            else
              if (      max(  (Alpha(kk)-mm1(kk)*(x(i)-0.5*hx(i)))/mm2(kk),(Alpha(kk)-mm1(kk)*(x(i)+0.5*hx(i)))/mm2(kk) ).lt. (y(j)-0.5*hy(j)) ) then
                cc(i,j)=0
              else if ( min ( (Alpha(kk)-mm1(kk)*(x(i)-0.5*hx(i)))/mm2(kk),(Alpha(kk)-mm1(kk)*(x(i)+0.5*hx(i)))/mm2(kk) ).gt. (y(j)+0.5*hy(j)) ) then 
                cc(i,j)=0
              else
                 
                if ( ( (x(i)-point(kk,1  ))**2 + (y(j)-point(kk,2  ))**2 )  .lt.( (point(kk,1)-point(kk+1,1))**2 + (point(kk,2)-point(kk+1,2))**2 ).AND. &
                     ( (x(i)-point(kk+1,1))**2 + (y(j)-point(kk+1,2))**2 )  .lt.( (point(kk,1)-point(kk+1,1))**2 + (point(kk,2)-point(kk+1,2))**2 ) ) then

                !     print*, " I am here2 " 
                  cc(i,j)=1.0
                  ccl(i,j)=kk
                  exit
                end if  
              end if 
            end if
 
          end do
       
        end if 
      end do

    end if
  end do

!OPEN(2350,file='BorderFinder.plt')

!write(2350,*) 'zone i=',nx-1,' k=',ny-1
!do j=1,ny-1 ;  Do i=1,nx-1

!write(2350,3500) x(i),y(j),cc(i,j)
!end do ; end do
!3500  format (3(1x,e15.7))
!call flush (2350)


return

end subroutine


subroutine InsideFinder(npoints,point,XBarV)
implicit none

double precision,intent(in)                            :: XBarV(2)
double precision,intent(in),dimension(npoints+1,2)     :: point
integer         ,intent(in)                            :: npoints
double precision                                       :: YCross
integer i,j,kk

  CI(:,:)=0

  do i=1, nx-1
    if ( abs (x(i)-XBarV(1)).gt.(2*r2) )  then
      cI(i,:)=0
    else
     
      do j=1,ny-1
        if ( abs (y(j)-XBarV(2)).gt.(len) ) then
          cI(i,j)=0
        else
          do kk=1,npoints
            
            if (mm2(kk).ne.0) then
                    YCross=(Alpha(kk)-mm1(kk)*x(i))/mm2(kk)
                    if (((YCross-point(kk  ,2))**2 + (X(i)-point(kk  ,1))**2).le.( (point(kk,1)-point(kk+1,1))**2 + (point(kk,2)-point(kk+1,2))**2 ).AND. &
                    &  ((YCross-point(kk+1,2))**2 + (X(i)-point(kk+1,1))**2).le.( (point(kk,1)-point(kk+1,1))**2 + (point(kk,2)-point(kk+1,2))**2 ).AND. &
                    &  ((Alpha(kk)-mm1(kk)*x(i))/mm2(kk)).gt.y(j)                                                                                      ) then
                    
                       cI(i,j)=cI(i,j)+1
                    end if  
            end if
          
          end do
       
          if (cI(i,j).eq.1.AND.cc(i,j).eq.0) then
            CI(i,j)=1.0
          else 
            CI(i,j)=0
          end if 
        end if 
      end do

    end if
  end do
 
         
!OPEN(2351,file='InsideFinder.plt')

!write(2351,*) 'zone i=',nx-1,' k=',ny-1
!do j=1,ny-1 ;  Do i=1,nx-1

!write(2351,3501) x(i),y(j),cI(i,j),cc(i,j)
!end do ; end do
!3501  format (4(1x,e15.7))
!call flush (2351)







 

return
end subroutine 


subroutine FractionFinder(npoints,point)
implicit none
double precision,intent(in),dimension(npoints+1,2)     :: point
double precision,dimension(npoints)                    :: mm1N,mm2N,alphaN
integer         ,intent(in)                            :: npoints
integer kk,i,j


  CF(:,:)=0 
  do i=1,nx-1
    do j=1,ny-1
      if (cc(i,j).ne.0) then
        kk=ccl(i,j) 
        if       ( (point(kk+1,1)-point(kk,1)).gt.0.AND.( point(kk+1,2)-point(kk,2) ).gt.0 ) then
          mm1N(kk)= mm1(kk)
          mm2N(kk)=-mm2(kk)
          AlphaN(kk)=Alpha(kk) -mm1(kk)*(x(i)-0.5*hx(i)) -mm2(kk)*(y(j)+0.5*hy(j))
          !print*, "hx,hy",hx(i),hy(j)
          !print*,"old value",mm1(1),mm2(1), Alpha(1)
          !print*,"new value",mm1N(1),mm2N(1), AlphaN(1)
          !print*, "x(i),y(j)", x(i),y(j) 
 
          !print*,  "I am here 1"
          if (AlphaN(kk).le.0.0) then
            mm1N(kk)=-mm1N(kk)
            mm2N(kk)=-mm2N(kk)
            AlphaN(kk)=-AlphaN(kk)
          end if 
        else if  ( (point(kk+1,1)-point(kk,1)).le.0.AND.( point(kk+1,2)-point(kk,2) ).ge.0 ) then
          mm1N(kk)= mm1(kk)
          mm2N(kk)= mm2(kk)
          AlphaN(kk)=Alpha(kk) -mm1(kk)*(x(i)-0.5*hx(i)) -mm2(kk)*(y(j)-0.5*hy(j))
        ! print *, " I am here 2"
          if (AlphaN(kk).le.0.0) then
            mm1N(kk)=-mm1N(kk)
            mm2N(kk)=-mm2N(kk)
            AlphaN(kk)=-AlphaN(kk)
          end if 
        else if ( (point(kk+1,1)-point(kk,1)).le.0.AND.( point(kk+1,2)-point(kk,2) ).le.0 ) then
        ! print*, "I am here 3 "
          mm1N(kk)=-mm1(kk)
          mm2N(kk)= mm2(kk)
          AlphaN(kk)=Alpha(kk) -mm1(kk)*(x(i)+0.5*hx(i)) -mm2(kk)*(y(j)-0.5*hy(j))
          if (AlphaN(kk).le.0.0) then
            mm1N(kk)=-mm1N(kk)
            mm2N(kk)=-mm2N(kk)
            AlphaN(kk)=-AlphaN(kk)
          end if 
        else
          !print*,  "I am here 4"
          mm1N(kk)=-mm1(kk)
          mm2N(kk)=-mm2(kk)
          AlphaN(kk)=Alpha(kk) -mm1(kk)*(x(i)+0.5*hx(i)) -mm2(kk)*(y(j)+0.5*hy(j))
          if (AlphaN(kk).le.0.0) then
            mm1N(kk)=-mm1N(kk)
            mm2N(kk)=-mm2N(kk)
            AlphaN(kk)=-AlphaN(kk)
          end if 
        end if
        
        if     (mm2N(kk).eq.0) then
            CF(i,j)=(AlphaN(kk)/mm1N(kk))/hx(i)
         !print*," I am here 5 "
        else if(mm1N(kk).eq.0) then
            CF(i,j)=(AlphaN(kk)/mm2N(kk))/hy(j)
        !print *,  "I am here 6" 
        else
          ! print*, "I am here 7"
            CF(i,j)=(AlphaN(kk)**2)/(2*mm1N(kk)*mm2N(kk))/(hx(i)*hy(j))*(1 &
                                                           & -HFF( AlphaN(kk)-mm1N(kk)*hx(i) )*( (AlphaN(kk)-mm1N(kk)*hx(i))/AlphaN(kk) )**2 &
                                                           & -HFF( AlphaN(kk)-mm2N(kk)*hy(j) )*( (AlphaN(kk)-mm2N(kk)*hy(j))/AlphaN(kk) )**2  )
                                                        
          !print*,"CF",CF(i,j) 
        end if
      
      end if
      if (cI(i,j).eq.1.0) then
         CF(i,j)=1.0
      end if 


    end do
  end do 
  

  !if (Test) then   
!write(2352,*) 'zone i=',nx-1,' k=',ny-1
!  do j=1,ny-1 ;  Do i=1,nx-1

!    write(2352,3502) x(i),y(j),cI(i,j),cc(i,j),CF(i,j),dble (ccl(i,j))
!  end do ; end do
!  3502  format (6(1x,e15.7))
!  call flush (2352)

!end if 

return

end subroutine 



Subroutine FractionSmoother()
implicit none
integer i,j,k

 CFS(:,:)=CF(:,:)  !! for boundary points !!


  if (smooth.ne.0) then
 
    do k=1,smooth
      do i=2,nx-2 ; do j=2,ny-2 
        CFS(i,j)= 0.25*(CF(i,j))+ &
        &   0.125*(CF(i+1,j)+CF(i-1,j)+CF(i,j-1)+CF(i,j+1)) + &
        &  0.0625*(CF(i+1,j+1)+CF(i+1,j-1)+CF(i-1,j+1)+CF(i-1,j-1))
      end do 
      end do 
    end do 
  end if 
  
return 
end subroutine   
          
      

end  module 


  

 
         
          
     
        
        
      
          


  

  
  
 


  
 


