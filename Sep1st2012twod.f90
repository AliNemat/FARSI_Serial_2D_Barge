
 program June5th2012twod
 use M_General_2D          
 use M_Platform_Constant_2D
 use M_Mesh_2D          
 use M_SolidFinder_2D 
 use M_output_2D  
    implicit none

parameter iterate=2
parameter npoints1=8
parameter npoints2=8

dOUBLE PRECISION Tx(2:nx-1,1:ny-1),Ty(1:nx-1,2:ny-1)

dOUBLE PRECISION miu,miudrop,miuair,rodrop,roair
double precision dt,div,divmax,beta,maxp,pdif
double precision gx,omegaz,MomentFx,MomentFy


double precision divmax2,Tmass,gap,eps
double precision omegazold,torder,teta,sumphi,sumphi0,yfree
double precision tolen
double precision tetaold,leg,ks,sumfp
double precision tetaplot,sumu,sumv,sumuold,sumvold,Lift,Liftold

double precision,dimension (2)                    ::   XBarV,XBarVOld,UBarV,UBarVOld
double precision,dimension (9)                    ::   YFreeG
double precision,dimension (5)                    ::   drag,dragold
double precision,dimension (1:3)                  ::   lnew,ox,oy
double precision,dimension (npoints1+1,2)         ::   point1,point1old,point2,point2old
double precision,dimension (1:nx,0:ny)            ::   advectu
double precision,dimension (0:nx,1:ny)            ::   advectv
double precision,dimension (0:nx+1,-1:ny+1)       ::   uold
double precision,dimension (-1:nx+1,0:ny+1)       ::   vold

double precision,dimension (-1:nx+1,-1:ny+1)      ::   phiold
double precision,dimension (1:nx-1,1:ny-1)        ::   mom,vort,CFS1,CFS1Old,CFS2,CFS2Old
double precision,dimension (0:nx,0:ny)            ::   ro,rt,p,miuv,miuw,row,div2

real HV
integer i,j,k,kk,tp,count,it,ADVECT,plot,solver,viscose,tstep,number2,tpstar,TCount


!! input parameters !!
call Platform_Constant_Ini() 
include 'par.txt'
!!!!!!!!!!!!!!!!!!!!!!!
Call meshgenerator()

 


 
do  i=1,nx-1
 
  if (x(i).le.floatx.AND.x(i+1).ge.floatx) then 
     gap=0.75*(x(i+1)-x(i))
     exit 
  end if 
end do 

eps=3*( y(int(ny/2))-y(int(ny/2)-1) )  !2*hy(ny/2)
count=0


XbarV(1)=FloatX
XBarV(2)=landa-0.115-0.0055
!xbar=lx/4.0  ; ybar=Ly/2.0

call ImportSolidPoint1(npoints1,point1,1)
call ImportSolidPoint1(npoints2,point2,2)
call SolidCellFinder(npoints1,point1,XBarV,CFS1)
call SolidCellFinder(npoints2,point2,XBarV,CFS2)

call Inicondition(p,yfree)
call boundarycond(0,dt,yfree,yfreeg)

!! end intial values !!

 
OPEN(25,file='solidlocation.plt') 
OPEN(35,file='result.plt')
OPEN(45,file='inside.plt')
!OPEN(75,file='tetherforce.plt')
OPEN(105,file='moment.plt')
OPEN(125,file='neweq.plt')
open(135,file='boundary condition.plt')
open(145,file='boundary condition2.plt')
open(155,file='Pressure contribution.plt')
open(165,file='WaveGage.plt')



!write(25,*) 'variables="t","momx","DTotal","DU","DAdv","DVisc","DP"'
write(25,'(A100)') 'variables="t","xbar","Ybar","teta (deg)","Omega","mass","DT","DUgh","DP"'
write(35,'(A100)') 'variables="x","y","u","v","ro","phi","vort","solid","p","div"'
write(45,'(A100)') 'variables="t","teta","sumphi","yfree","yfreeg(1)"'
!write(75,'(A100)') 'variables="time","ftether1","ftether2"'
write(105,'(A50)') 'variables="tp","phit"'
write(125,'(A75)') 'variables="x","y","uext","vext","phi","ib"'
write(135,'(A75)') 'variables="x","y","uimp","vimp","ib"'
write(145,'(A50)') 'variables="time","sum"'
write(155,'(A75)') 'variables="x","y","indexi","indexj","Ib","Dp"'
write(155,'(A50)') 'zone i=',nx-2,' j=',ny-2
write(165,'(A150)') 'variables="t","Gage1","Gage2","Gage3","Gage4","Gage5","Gage6","Gage7","Gage8","Gage9","Gage10","Gage11","Gage12","Gage13","Gage14"'

!1100  format(A150)
sumu=0 ; sumv=0 ; Drag(:)=0 ; Lift=0 
      !! main loop of the code !!! 
      do tp=1,tstep !time step
      print*,tp
     
      !! parameters for Make solution second order in time !! 
      uold=u  ; vold=v  
      phiold=phi 
      !Icold=Ic 
      
      CFS1old(:,:)=CFS1(:,:)
      Point1Old(:,:)=Point1(:,:)

      CFS2old(:,:)=CFS2(:,:)
      Point2Old(:,:)=Point2(:,:)
      XBarVOld(:)=XBarV(:)       
      UBarVOld(:)=UBarV(:)
      omegazold=omegaz
      
             
      tetaold=teta
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      sumuold=sumu ; sumvold=sumv  
      Dragold(:)=Drag(:) ; Liftold=Lift
      !! Second order time step !! 
      do torder=1,2
        tpstar=tp-1+torder
          
        
    


 


divmax2=1



 
 
  
 !! density and viscosity definition !! 
     
   ! call property(Ib,Ic,ro,rodrop,roair,rosolid,rocon,gap)
   ! call property(Ib,Ic,miuv,miudrop,miuair,miusolid,miucon,gap)
     
     do i=0,nx ;do j=0,ny 
     ro(i,j)=roair+HV( -phi(i,j),eps )*(rodrop-roair)
     ro(i,j)=ro(i,j)+CFS1(i,j)*( rosolid-ro(i,j) )
     ro(i,j)=ro(i,j)+CFS2(i,j)*( rocon )
     end do ;end do 
!     
!   

   
     do i=0,nx ;do j=0,ny 
     miuv(i,j)=miuair+HV( -phi(i,j),eps )*(miudrop-miuair)
     miuv(i,j)=miuv(i,j)+CFS1(i,j)*( miusolid-miuv(i,j) )
     miuv(i,j)=miuv(i,j)+CFS2(i,j)*( miucon )
     end do ; end do 
      
! writing the initial domian parameters!!       
if (tpstar.eq.1.) then     
call Output_vtk(tp,ro,CFS1)
write(35,*) 'zone i=',nx-1,' j=',ny-1
 Do j=1,ny-1 ; Do i=1,nx-1
write(35,122) x(i),y(j),0.5*(u(i,J)+u(i+1,j)),0.5*(v(i,j)+V(i,j+1)),ro(i,j),phi(i,j),0.0,CFS1(i,j),p(i,j),0.0
 122  format (10(1x,e15.7))
end do ; end do  
end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! advection and viscose terms in Navier-Stokes equation!!! 
call advection (advectu,advectv)
call viscosity (ro,miuv,Tx,Ty)

Do i=2,nx-1 ; Do j=1, ny-1
U(i,j)=U(i,j)+advectu(i,j)*dt+Tx(i,j)*dt+gx*dt

end do ; end do
Do j=2,ny-1 ; Do i=1,nx-1
V(i,j)=V(i,j)+advectv(i,j)*dt+Ty(i,j)*dt+gy*dt
end do ; end do 

call boundarycond(tpstar,dt,yfree,yfreeg) !! ? !!


call poisson (ro,dt,pdif,p,beta)        

 !!! start update velocity !!!!!!!!!!!!!!!!!!!
do i=2,nx-1 ; do j=1,ny-1 

U(i,j)=U(i,j)-(  1/(  0.5*(ro(i,j)+ro(i-1,j))  )   )*( p(i,j)-p(i-1,j) )*dt/( x(i)-x(i-1)  ) 
end do ; end do 

do i=1,nx-1 ; do j=2,ny-1 
V(i,j)=V(i,j)-(  1/(  0.5*(ro(i,j)+ro(i,j-1))  )   )*( p(i,j)-p(i,j-1) )*dt/( y(j)-y(j-1)  )
end do ; end do
!! 3 !!

 
!!!!! computing maximum Divergence in the domain before  Solid correction in "Solid Subroutine !!!
divmax2=0
do i=2,nx-1 ;do j=2,ny-1 
div=(U(i+1,j)-U(i,j))/( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )+ &
&   (v(i,j+1)-v(i,j))/( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )

if(abs(div).gt.divmax2)then
divmax2=abs(div)
end if 

end do ; end do 

WRITE(*,*)'DIV before  solid Correct=',DIVMAX2 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Drag(:)=0 ; Lift=0 
MomentFX=0 ; MomentFY=0 


sumfp=0
 do j=2,ny-1 ; do i=2,nx-1  


Drag(1)=Drag(1)+ 0.5*(ro(i,j)+ro(i-1,j))*0.5*(CFS1(i,j)+CFS1(i-1,j))*(x(i)-x(i-1))*( 0.5*(y(j)+y(j+1))-0.5*(y(j)+y(j-1)) )* &
&          ( +advectu(i,j)-(  1/(  0.5*(ro(i,j)+ro(i-1,j))  )   )*( p(i,j)-p(i-1,j) )/( x(i)-x(i-1)  )+Tx(i,j)+gx )

Lift=Lift+ 0.5*(ro(i,j)+ro(i,j-1))*0.5*(CFS1(i,j)+CFS1(i,j-1))*(y(j)-y(j-1))*( 0.5*(x(i)+x(i+1))-0.5*(x(i)+x(i-1)) )* &
&          ( +advectv(i,j)-(  1/(  0.5*(ro(i,j)+ro(i,j-1))  )   )*( p(i,j)-p(i,j-1) )/( y(j)-y(j-1)  )+Ty(i,j)+gy )

MomentFY=MomentFY+(x(i)-XBarV(1))* &
&              0.5*(ro(i,j)+ro(i,j-1))*0.5*(CFS1(i,j)+CFS1(i,j-1))*(y(j)-y(j-1))*( 0.5*(x(i)+x(i+1))-0.5*(x(i)+x(i-1)) )*   &
&              ( +advectv(i,j)-(  1/(  0.5*(ro(i,j)+ro(i,j-1))  )   )*( p(i,j)-p(i,j-1) )/( y(j)-y(j-1)  )+Ty(i,j)+gy ) 

MomentFX=MomentFX-(y(j)-YFreeIni)* &
&             0.5*(ro(i,j)+ro(i-1,j))*0.5*(CFS1(i,j)+CFS1(i-1,j))*(x(i)-x(i-1))*( 0.5*(y(j)+y(j+1))-0.5*(y(j)+y(j-1)) )* &
&             ( +advectu(i,j)-(  1/(  0.5*(ro(i,j)+ro(i-1,j))  )   )*( p(i,j)-p(i-1,j) )/( x(i)-x(i-1)  )+Tx(i,j)+gx )
end do ; end do
 



 !! immersed boundary method for solid correction !!    
      call solid(ro,XBarV,UBarV,dt,omegaz,gy,Tmass,teta,tp,sumu,sumv,CFS1)

      XbarV(:)=XbarV(:)+dt*UbarV(:)
      
      call SolidUpdatePosition(npoints1,point1,UBarV,XBarV,omegaz,dt)
      call SolidUpdatePosition(npoints2,point2,UBarV,XBarV,omegaz,dt)

      call SolidCellFinder(npoints1,point1,XBarV,CFS1)
      call SolidCellFinder(npoints2,point2,XBarV,CFS2)

      call levelset(dt)
      call reinitialize(eps,dt)
 
!
divmax2=0
do i=2,nx-1 ;do j=2,ny-1 
div=(U(i+1,j)-U(i,j))/( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )+ &
&   (v(i,j+1)-v(i,j))/( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )

if(abs(div).gt.divmax2)then
divmax2=abs(div)
end if 

end do ; end do 

WRITE(*,*)'DIVERGENCEiterate=',DIVMAX2 


end do   !! torder !!
!!!!!!!!!!!!!!!!!!! averaging to consequence time step for reaching second order in time for velocity field !!!!!!!!!!!

!! fluid values !!
u=0.5d0*(u+uold)                   ; v=0.5d0*(v+vold)  
!! property  of solid and fluid (density and viscosity)               
phi=0.5d0*(phi+phiold)           
!Ic=0.5d0*(Icold+Ic)  
 
 CFS1(:,:)=0.5d0*(CFS1old(:,:)+CFS1(:,:))
 CFS2(:,:)=0.5d0*(CFS2old(:,:)+CFS2(:,:))

point1(:,:)=0.5d0*(point1Old(:,:)+point1(:,:))  
point2(:,:)=0.5d0*(point2Old(:,:)+point2(:,:))  
!!!!
!! solid and tether positions/velocity/accelaration !! 
xbarV(:)=0.5d0*(xbarVold(:)+xbarV(:))           
ubarV(:)=0.5d0*(ubarVold(:)+ubarV(:))    
omegaz=0.5d0*(omegazold+omegaz)
 

!! just for plotting !!!       
 teta=0.5*(tetaold+teta)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


sumu=0.5*(sumu+sumuold) 
sumv=0.5*(sumv+sumvold)

Drag(:)=0.5*(Drag(:)+Dragold(:))
Lift=0.5*(Lift+Liftold)




!! Post processing start !! 


!!!!!!!!!!!!!!!!!! computing divergence at the end of marching 1 step in time !!!!!!!!!!!!!!!!!!!!!!!!!!!
divmax2=0
div2=0
do i=2,nx-1 ;do j=2,ny-1 
div2(i,j)=(U(i+1,j)-U(i,j))/( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )+ &
      &   (v(i,j+1)-v(i,j))/( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )

!if(abs(div2(i,j)).gt.divmax2)then
!divmax2=abs(div2(i,j))
!else
!end if 

end do ; end do
WRITE(*,*)'DIVERGENCE Second order in time =', divmax2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!! writing  1 !!!!!!!!!!!!!!!!!!!!
tetaplot=180/pi*Atan( ( point1(5,2)-point1(6,2) )/(point1(5,1)-point1(6,1))  )
write (25,'(9(1x,e15.7))') tp*dt,XBarV(1),XBarV(2),tetaplot,omegaz,Drag(1),Lift,MomentFX,MomentFY
 !write (25,*) tp*dt,xbar,ybar,tetaplot,(sumu-0.0)/dt/(0.5*rodrop*(0.375**2)*2*r2),(sumv-0.0)/dt/(0.5*rodrop*(0.375**2)*2*r2)
! write (25,*) tp*dt,xbar,(sumu)/dt/(0.5*rodrop*(0.375**2)*2*r2),(sumv)/dt/(0.5*rodrop*(0.375**2)*2*r2),&
! &            Drag/(0.5*rodrop*(0.375**2)*2*r2),Lift/(0.5*rodrop*(0.375**2)*2*r2)
! write (25,*) tp*dt,(sumu)/dt/(0.5*rodrop*(0.375**2)*2*r2),Drag(1)/(0.5*rodrop*(0.375**2)*2*r2),Drag(2)/(0.5*rodrop*(0.375**2)*2*r2),&
! &            Drag(3)/(0.5*rodrop*(0.375**2)*2*r2),Drag(4)/(0.5*rodrop*(0.375**2)*2*r2),Drag(5)/(0.5*rodrop*(0.375**2)*2*r2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i=1,nx-1 ; do j=1,ny-1 
vort(i,j)= 0.25*(( v(i,j)    -v(i-1,j)   )/(x(i)-x(i-1))   &
        &       +( v(i+1,j)  -v(i,j)     )/(x(i+1)-x(i))   &
        &       +( v(i,j+1)  -v(i-1,j+1) )/(x(i)-x(i-1))   &
        &       +( v(i+1,j+1)-v(i,j+1)   )/(x(i+1)-x(i))  )&
        & -0.25*(( u(i,j)    -u(i,j-1)   )/(y(j)-y(j-1))   &
        &       +( u(i+1,j)  -u(i+1,j-1) )/(y(j)-y(j-1))   &
        &       +( u(i,j+1)  -u(i,j)     )/(y(j+1)-y(j))   &
        &       +( u(i+1,j+1)-u(i+1,j)   )/(y(j+1)-y(j))  ) 
        
end do ; end do  

!!!!!!!!!!!!!!! writing  2 !!!!!!!!!!!!!!!!!!!!
count=count+1
if (count.eq.plot.AND.DIVMAX2.lt.500)then 
  count=0
 call Output_vtk(tp,ro,CFS1)
  print*,"data is wrriten"
  write(35,*) 'zone i=',nx-1,' j=',ny-1
  Do j=1,ny-1 ; Do i=1,nx-1
    write(35,122) x(i),y(j),0.5*(u(i,J)+u(i+1,j)),0.5*(v(i,j)+V(i,j+1)),ro(i,j),phi(i,j),vort(i,j),CFS1(i,j),p(i,j),abs(div2(i,j))
  end do ; end do 

end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!! writing  3 !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! checking mass conservation in Level Set method !!!!!!!!!!!!!!!!!!!
sumphi=0 
do i=1,nx-1 ; do j=1,ny-1 
if (phi(i,j).gt.0) then 
sumphi=sumphi+1
end if 
end do ; end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
write(45,'(5(1x,e15.7))') tp*dt,teta*180/pi,sumphi,yfree,yfreeg(1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!! write for tethers 4 (not now )!!!!!!!!!!!!!!!!!!!!!!!
!write (75,'(3(1x,e15.7))') tp*dt, ks*( max(Lzero(1),Lnew(1))-Lzero(1) ),ks*( max(Lzero(2),Lnew(2))-Lzero(2) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
write(165,1200) tp*dt,YFreeG(1),YFreeG(2),YFreeG(3),YFreeG(4),YFreeG(5),YFreeG(6),YFreeG(7),YFreeG(8),YFreeG(9), &
&                  YFreeG(10),YFreeG(11),YFreeG(12),YFreeG(13),YFreeG(14)

1200  format (15(1x,e15.7))
!! Post processing End  !! 


 

end do !time step !! end of the main loop of the code !! 

write(*,*)'end'
read (*,*)

    end program 




!!!!!!!!!!!!!!!!!!!!!!! subroutines!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




Subroutine Inicondition(p,yfree)
use M_General_2D,  only: nx,ny,y,u,v,phi


implicit none 

double precision p(0:nx,0:ny)

double precision yfree

integer i,j



phi=10000

 do i=-1,nx+1
  do j=-1,ny+1

   phi(i,j)=y(j)-yfree

  end do 
 end do 

v=0 ; U=0 ; p=0 !!initial values for velocity and pressure 



return 
end 



subroutine boundarycond(tpstar,dt,yfree,yfreeg)
use M_General_2D,  only : phi,u,v,nx,ny,x,y,pi,froude,gy,lx,landa
implicit none 
integer i,j,ii,tpstar,imax,jmax,Cb,igage
double precision dt
double precision period,hait,HH,HH1,yfree,kwave,wwave,sumvv,wavegen
double precision,dimension (14) ::  xgage,YFreeG
!! boundary condition for Velocity (Including generating wave ) and Level Set function !! 
!! generated wave from Stokes second order solution  but only velocity in x direction is given velocity in y direction will be 
!! satisifed by mass conservation  
!! inflow out flow on left and top !! 
!! Slip boundary condition on the others !! 





call  boundarycondphi()

include "Par_Wave_2D.txt"

period=landa/(  sqrt(-gy*hait)*sqrt( tanh(2*pi/landa)/(2*pi/landa*hait) )   )
kwave=2*pi/landa
wwave=2*pi/period


do ii=1,size (XGage)
!!!!!!!!!!!!!!!!!! finding the wave height at the  wave gage place !!!!!!!!!!!!!!!!!
igage=1 !! for problems with out free surface 
do  i=1,nx-1
 
  if (x(i).le.xgage(ii).AND.x(i+1).ge.xgage(ii)) then 
     igage=i
     exit 
  end if 
end do 
yfreeg(ii)=0
do j=1,ny-1
  if (0.5*( phi(igage+1,j)+phi(igage,j) ).le.0.AND.0.5*( phi(igage+1,j+1)+phi(igage,j+1) ).ge.0 ) then             !! may be more precsiness can be done !! 
    yfreeg(ii)=y(j)+ ( y(j+1)-y(j) )/ ( 0.5*( phi(igage+1,j+1)+phi(igage,j+1) )-0.5*( phi(igage+1,j)+phi(igage,j) ) ) & 
    &    *( -0.5*( phi(igage+1,j)+phi(igage,j) ) ) 
    exit 
  end if 
end do 

if (yfreeg(ii).eq.0)then 
print*,"errorgage(ii)" 
end if
end do 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 
!!!!!!!!! indicator for error !!!!!!!

 


!!!!!!!!!!!!!!! tracking the free surface place in left hand side of the domain for using in wave generator !!!!!!!!!!!!!
yfree=0
do j=1,ny-1
  if (0.5*( phi(1,j)+phi(0,j) ).le.0.AND.0.5*( phi(1,j+1)+phi(0,j+1) ).ge.0 ) then             !! may be more precsiness can be done !! 
    yfree=y(j)+ ( y(j+1)-y(j) )/ ( 0.5*( phi(1,j+1)+phi(0,j+1) )-0.5*( phi(1,j)+phi(0,j) ) ) *( -0.5*( phi(1,j)+phi(0,j) ) ) 
    exit 
  end if 
end do  
if (yfree.eq.0)then 
print*,"error" 
end if





     

!!!!!! analytical solution for second order wave, the velocity for the air is fake (only for continuty )

!if ( (tpstar*dt).lt.(1*period) ) then 
!  HH=0

!else
 if ( (tpstar*dt).lt.(1*period) ) then 
  !HH=(tpstar*dt-period)/(1*period)*HH1

  HH=(tpstar*dt)/(1*period)*HH1
else 
HH=HH1
end if

print*,"time/T=",tpstar*dt/period,"Waveheight percentage=",HH/HH1 
sumvv=0
if (wavegen.eq.1 ) then 

  do j=1 ,ny-1 

    if (phi(0,j).lt.0 ) then

    u(1,j)=HH/2* wwave *cos(-wwave*(tpstar*dt))*cosH(kwave* ( (y(j)-yfree)+hait) )/sinH(kwave*hait)+ &
    & 3/4*((pi*HH)**2)/(period*landa)*cos(-2*wwave*(tpstar*dt))*cosH( 2*kwave* ((y(j)-yfree)+hait) )/( (sinH(kwave*hait))**4 )
 
    u(0,j)=HH/2* wwave*cos(-wwave*(tpstar*dt))*cosH(kwave* ( (y(j)-yfree)+hait) )/sinH(kwave *hait)+ &
    & 3/4*((pi*HH)**2)/(period*landa)*cos(-2*wwave*(tpstar*dt))*cosH( 2*kwave* ((y(j)-yfree)+hait) )/( (sinH(kwave*hait))**4 )
                        
    else

    u(1,j)=HH/2* wwave *cos(-wwave*(tpstar*dt))*cosH(kwave* ( -(y(j)-yfree)+hait) )/sinH(kwave*hait)+ &
    & 3/4*((pi*HH)**2)/(period*landa)*cos(-2*wwave*(tpstar*dt))*cosH( 2*kwave* (-(y(j)-yfree)+hait) )/( (sinH(kwave*hait))**4 )
 
    u(0,j)=HH/2* wwave*cos(-wwave*(tpstar*dt))*cosH(kwave* ( -(y(j)-yfree)+hait) )/sinH(kwave *hait)+ &
    & 3/4*((pi*HH)**2)/(period*landa)*cos(-2*wwave*(tpstar*dt))*cosH( 2*kwave* (-(y(j)-yfree)+hait) )/( (sinH(kwave*hait))**4 )
    
    end if 
 sumvv=sumvv+u(1,j)*0.5*(y(j+1)-y(j-1))*1.0d0 

  !U(1,j)=0      
  !U(0,j)=-U(2,j)
   u(nx,j)=0 
   u(nx+1,j)=-u(nx-1,j) 
  end do 

else if (wavegen.eq.2) then 

   do j=1 ,ny-1 

     if (phi(0,j).gt.0 ) then
     
     u(1,j)=HH*y(j)*sin(kwave*0-wwave*tpstar*dt)      
     u(0,j)=HH*y(j)*sin(kwave*-2*x(1)-wwave*tpstar*dt)
     
     else 
     
     u(1,j)=HH*( 2*yfree-y(j)   )*sin(kwave*0-wwave*tpstar*dt)     
     u(0,j)=HH*( 2*yfree-y(j)   )*sin(kwave*-2*x(1)-wwave*tpstar*dt)
     
     end if 
 !  sumvv=sumvv+u(1,j)*hy(j)*1.0d0
    sumvv=sumvv+u(1,j)*0.5*(y(j+1)-y(j-1))*1.0d0 
   !U(1,j)=0      
   !U(0,j)=-U(2,j)
   u(nx,j)=0 
   u(nx+1,j)=-u(nx-1,j) 
  end do
else if (wavegen.eq.3) then 
 
  do j=1 ,ny-1 
   U(1,j)=0      
   U(0,j)=-U(2,j)
   u(nx,j)=0 
   u(nx+1,j)=-u(nx-1,j) 
  end do
  
 else  !! unifrom inflow and not gradient outflow !!
 do j=1 ,ny-1 
   u(1,j)   =0.375 !u(nx-1,j)      
   u(0,j)   =0.375 !u(nx-2,j)
   u(nx,j)  =u(nx-1,j) 
   u(nx+1,j)=u(nx-2,j) 
  end do
  
  

end if 

!! at first smaller and then use them for boundary of bigger !!!


!!!slip!! 
Do i=0, nx+1
U(i,0)=U(i,1)
U(i,-1)=U(i,2)
U(i,ny)=U(i,ny-1)
U(i,ny+1)=U(i,ny-2)
end do

!no-slip!! cavity 
!Do i=0, nx+1
!U(i,0)=-U(i,1)
!U(i,-1)=-U(i,2)
!U(i,ny)=-U(i,ny-1)  !!1.0 
!U(i,ny+1)=-U(i,ny-2) !!1.0
!end do



if (wavegen.eq.1.OR.wavegen.eq.2) then
   Do i=1, nx-1 
   V(i,1)=0
   V(i,0)=-V(i,2)
   V(i,ny)=sumvv/( real(nx-1) *0.5*(x(i+1)-x(i-1)) )         !!!!!!!!!!!!!inflow out flow  boundary condition !!!!!!!!!!!!!!
   V(i,ny+1)=V(i,ny)
   end do
else
   Do i=1, nx-1 
   V(i,1)=0
   V(i,0)=-V(i,2)
   V(i,ny)=0
   V(i,ny+1)=-V(i,ny-1)
   end do
end if 



if (wavegen.eq.5) then
 
   do j=0,ny+1
   
   if (phi(0,j).lt.0 ) then
   v(0,j)=HH/2*wwave *sin(-wwave*(dble(tpstar)*dt))*sinH(kwave* ( (y(j)-yfree)+hait) )/sinH(kwave *hait)+ &
   & 3/4*((pi*HH)**2)/(period*landa)*sin(-2*wwave*(tpstar*dt))*sinH( 2*kwave* ((y(j)-yfree)+hait) )/( (sinH(kwave*hait))**4 )

   v(-1,j)=HH/2* wwave *sin(-wwave*(dble(tpstar)*dt))*sinH(kwave* ( (y(j)-yfree)+hait) )/sinH(kwave *hait)+ &
   & 3/4*((pi*HH)**2)/(period*landa)*sin(-2*wwave*(tpstar*dt))*sinH( 2*kwave* ((y(j)-yfree)+hait) )/( (sinH(kwave*hait))**4 )
   else 
   v(0,j)=HH/2*wwave *sin(-wwave*(dble(tpstar)*dt))*sinH(kwave* ( -(y(j)-yfree)+hait) )/sinH(kwave *hait)+ &
   & 3/4*((pi*HH)**2)/(period*landa)*sin(-2*wwave*(tpstar*dt))*sinH( 2*kwave* (-(y(j)-yfree)+hait) )/( (sinH(kwave*hait))**4 )

   v(-1,j)=HH/2* wwave *sin(-wwave*(dble(tpstar)*dt))*sinH(kwave* ( -(y(j)-yfree)+hait) )/sinH(kwave *hait)+ &
   & 3/4*((pi*HH)**2)/(period*landa)*sin(-2*wwave*(tpstar*dt))*sinH( 2*kwave* (-(y(j)-yfree)+hait) )/( (sinH(kwave*hait))**4 )
   end if
   
   V(nx,j)=V(nx-1,j)
   V(nx+1,j)=V(nx-2,j)
   
   end do 
else

!! slip !!
   do j=0,ny+1 
   
   V(0,j)=V(1,j)
   V(-1,j)=V(2,j)
   V(nx,j)=V(nx-1,j)
   V(nx+1,j)=V(nx-2,j)
   
   end do

!! No-Slip cavity !!
!do j=0,ny+1 
!   
!   V(0,j)=-V(1,j)
!   V(-1,j)=-V(2,j)
!   V(nx,j)=-V(nx-1,j)
!   V(nx+1,j)=-V(nx-2,j)
!   
! end do

! periodic !!!!!!!!!
!do j=0,ny+1 
!   
!   V(nx,j)=V(1,j)
!   V(nx+1,j)=V(2,j)
!   V(0,j)=V(nx-1,j)
!   V(-1,j)=V(nx-2,j)
!   
!   end do

!! out flow 
!do j=0,ny+1 
   
!   V(0,j)=V(1,j)
!   V(-1,j)=V(2,j)
!   V(nx,j)=V(nx-1,j)
!   V(nx+1,j)=V(nx-2,j)
   
!   end do
end if    
 


return 
end 

subroutine levelset(dt)
use M_General_2D,  only:phi,u,v,x,y,nx,ny 
implicit none

 
Double precision phioldc(-1:nx+1,-1:ny+1)
double precision, DIMENSION (0:nx,0:ny)       ::Dpphix,Dmphix,Dpphiy,Dmphiy
double precision, DIMENSION (1:nx-1,1:ny-1)   ::phix,phiy,Lphin
Double precision Lphis,dt
integer i,j,kk

!if (numimp.eq.2) then 
phioldc(1:nx-1,1:ny-1)=phi(1:nx-1,1:ny-1)
!end if 


!! Second ordeer ENO convective terms for Advection terms of Level set equation !! 
!! Second order in time !! 

do kk=1,2  !!prediction correction method!!


do i=0,nx ; do j=0,ny                                
Dpphix(i,j)=( phi(i+1,j)-phi(i,j)   )/(  x(i+1)-x(i)  )
Dmphix(i,j)=( phi(i,j)  -phi(i-1,j) )/(  x(i)-x(i-1)  )
Dpphiy(i,j)=( phi(i,j+1)-phi(i,j)   )/(  y(j+1)-y(j)  )
Dmphiy(i,j)=( phi(i,j)  -phi(i,j-1) )/(  y(j)-y(j-1)  )

end do ;end do 

do i=1,nx-1 ; do j=1,ny-1    !!A!!

if (0.5*( u(i,j)+u(i+1,j) ).gt.0.d0) then


  if (  abs(  Dmphix(i,j)-Dmphix(i-1,j) ).lt.abs(  Dpphix(i,j)-Dpphix(i-1,j) )   ) then
phix(i,j)=Dmphix(i,j)+  0.5*(  Dmphix(i,j)-Dmphix(i-1,j)  ) 
   else
phix(i,j)=Dmphix(i,j)+  0.5*(  Dpphix(i,j)-Dpphix(i-1,j)  )
   end if 
    
else

   if (  abs(  Dmphix(i+1,j)-Dmphix(i,j) ).lt.abs(  Dpphix(i+1,j)-Dpphix(i,j) )   ) then
phix(i,j)=Dpphix(i,j)-  0.5*(  Dmphix(i+1,j)-Dmphix(i,j)  ) 
   else
phix(i,j)=Dpphix(i,j)-  0.5*(  Dpphix(i+1,j)-Dpphix(i,j)  )
   end if 

end if 



if (0.5*( V(i,j)+V(i,j+1) ).gt.0.d0) then


   if (  abs(  DmphiY(i,j)-DmphiY(i,j-1) ).lt.abs(  DpphiY(i,j)-DpphiY(i,j-1) )   ) then
phiY(i,j)=DmphiY(i,j)+  0.5*(  DmphiY(i,j)-DmphiY(i,j-1)  ) 
   else
phiY(i,j)=DmphiY(i,j)+  0.5*(  DpphiY(i,j)-DpphiY(i,j-1)  )
   end if 
    
else

  if (  abs(  DmphiY(i,j+1)-DmphiY(i,j) ).lt.abs(  DpphiY(i,j+1)-DpphiY(i,j) )   ) then
phiY(i,j)=DpphiY(i,j)-  0.5*(  DmphiY(i,j+1)-DmphiY(i,j)  ) 
   else
phiY(i,j)=DpphiY(i,j)-  0.5*(  DpphiY(i,j+1)-DpphiY(i,j)  )
   end if 


end if 


     
end do ;end do  !!A!!

if (kk.eq.1) then

do i=1,nx-1 ; do j=1,ny-1 
Lphin(i,j)=  -0.5*( u(i,j)+u(i+1,j) )* phix(i,j)+&
         &   -0.5*( V(i,j)+V(i,j+1) )* phiy(i,j) 
           
phi(i,j)=phioldc(i,j)+dt*( -0.5*( u(i,j)+u(i+1,j) )* phix(i,j)+&
                     & -0.5*( V(i,j)+V(i,j+1) )* phiy(i,j) )
                
end do ;end do 

                  
end if                        

   if (kk.eq.2) then
   
   do i=1,nx-1 ; do j=1,ny-1 
 
Lphis=   -0.5*( u(i,j)+u(i+1,j) )* phix(i,j)&
       & -0.5*( V(i,j)+V(i,j+1) )* phiy(i,j) 
            
phi(i,j)=phioldc(i,j)+0.5*dt*(Lphis+Lphin(i,j) )                       
   end do ;end do 

   end if 


end do  !!prediction corection method!! 
 
call boundarycondphi()                    

RETURN
END 


subroutine reinitialize(eps,dt)
use M_General_2D,  only:phi,u,v,x,y,nx,ny,pi 
implicit none
Double precision eps
Double precision phixm,phiym,phixp,phiyp,phixr,phiyr,sphi,s,dtau,AA
double precision, Dimension (0:nx,0:ny)       ::  Dpphix,Dmphix,Dpphiy,Dmphiy
double precision, Dimension (1:nx-1,1:ny-1)   ::  phioldCR,phin,lphin
real sgn
Double precision lphis,bb,avgh,dt
integer i,j,kk,m

! Ref; Journal of computational physics vol. 152 pp-493-516 (1999)

phioldCR(1:nx-1,1:ny-1)=phi(1:nx-1,1:ny-1)
avgh=1000 !! large number!! 
do i=-1,nx ; do j=-1, ny 
bb=min (x(i+1)-x(i),y(j+1)-y(j))
if (bb.lt.avgh) then
avgh=bb
end if 
end do ; end do 

dtau=0.5*avgh  !! fictious time step !


do kk=1,3

do m=1,2   !!prediction correction method for time step the same as Level set for advection terms  !! 


do i=0,nx ; do j=0,ny 
Dpphix(i,j)=( phi(i+1,j)-phi(i,j)   )/(  x(i+1)-x(i)  )
Dmphix(i,j)=( phi(i,j)  -phi(i-1,j) )/(  x(i)-x(i-1)  )
Dpphiy(i,j)=( phi(i,j+1)-phi(i,j)   )/(  y(j+1)-y(j)  )
Dmphiy(i,j)=( phi(i,j)  -phi(i,j-1) )/(  y(j)-y(j-1)  )

end do ;end do 

do i=1,nx-1 ; do j=1,ny-1 !!A!!  !!A!!


!sphi=phiold(i,j)/(  dsqrt(phiold(i,j)*phiold(i,j)+h*h)  )
!sphi=sgn(phiold)



if (phioldCR(i,j).gt.eps) then
sphi=1.d0
else if (phioldCR(i,j).lt.-eps) then
sphi=-1.d0
else
sphi=phioldCR(i,j)/eps -(1/pi)*dsin(pi*phioldCR(i,j)/eps) 
end if 

 

   if (  abs(  Dmphix(i,j)-Dmphix(i-1,j) ).lt.abs(  Dpphix(i,j)-Dpphix(i-1,j) )   ) then
phixm=Dmphix(i,j)+  0.5*(  Dmphix(i,j)-Dmphix(i-1,j)  ) 
   else
phixm=Dmphix(i,j)+  0.5*(  Dpphix(i,j)-Dpphix(i-1,j)  )
   end if 
    

    if (  abs(  Dmphix(i+1,j)-Dmphix(i,j) ).lt.abs(  Dpphix(i+1,j)-Dpphix(i,j) )   ) then
phixp=Dpphix(i,j)-  0.5*(  Dmphix(i+1,j)-Dmphix(i,j)  ) 
   else
phixp=Dpphix(i,j)-  0.5*(  Dpphix(i+1,j)-Dpphix(i,j)  )
   end if 



    if (  abs(  DmphiY(i,j)-DmphiY(i,j-1) ).lt.abs(  DpphiY(i,j)-DpphiY(i,j-1) )   ) then
phiYm=DmphiY(i,j)+  0.5*(  DmphiY(i,j)-DmphiY(i,j-1)  ) 
   else
phiYm=DmphiY(i,j)+  0.5*(  DpphiY(i,j)-DpphiY(i,j-1)  )
   end if 
    


   if (  abs(  DmphiY(i,j+1)-DmphiY(i,j) ).lt.abs(  DpphiY(i,j+1)-DpphiY(i,j) )   ) then
phiYp=DpphiY(i,j)-  0.5*(  DmphiY(i,j+1)-DmphiY(i,j)  ) 
   else
phiYp=DpphiY(i,j)-  0.5*(  DpphiY(i,j+1)-DpphiY(i,j)  )
   end if
  
  
   

   
        if (sphi*phixp.ge.0.d0.AND.sphi*phixm.ge.0.d0) then
   phixr=phixm
   else if (sphi*phixp.le.0.d0.AND.sphi*phixm.le.0.d0) then 
   phixr=phixp
   else if (sphi*phixp.gt.0.d0.AND.sphi*phixm.lt.0.d0) then
   phixr=0.0
   else if (sphi*phixp.lt.0.d0.AND.sphi*phixm.gt.0.d0) then
   s=sphi*( abs(phixp)-abs(phixm) )/(phixp-phixm)
                      if (s.gt.0.0) then
                      phixr=phixm
                   else
                      phixr=phixp 
                   end if
   end if
   
   
        if (sphi*phiyp.ge.0.d0.AND.sphi*phiym.ge.0.d0) then
   phiyr=phiym
   else if (sphi*phiyp.le.0.d0.AND.sphi*phiym.le.0.d0) then 
   phiyr=phiyp
   else if (sphi*phiyp.gt.0.d0.AND.sphi*phiym.lt.0.d0) then
   phiyr=0.0
   else if (sphi*phiyp.lt.0.d0.AND.sphi*phiym.gt.0.d0) then
   s=sphi*( abs(phiyp)-abs(phiym) )/(phiyp-phiym)
                      if (s.gt.0.0) then
                      phiyr=phiym
                   else
                      phiyr=phiyp 
                  end if
    end if
    
        
   
   
   if (m.eq.1) then
   lphin(i,j)=sphi*(  1.d0-dsqrt(phiyr*phiyr+phixr*phixr)  )
   phin(i,j)=phi(i,j)
   phi(i,j)=phi(i,j)+dtau*sphi*(  1.d0-dsqrt(phiyr*phiyr+phixr*phixr)  )
   
   end if
   
   
   if (m.eq.2) then
   lphis   =sphi*(  1.d0-dsqrt(phiyr*phiyr+phixr*phixr)  )
   phi(i,j)=phin(i,j)+0.5*dtau*(  lphis+lphin(i,j)  )                   
   end if 

end do ; end do 

call boundarycondphi() 


end do 



end do !!fictious time step!!

return 
END

subroutine solid(ro,XBarV,ubarV,dt,omegaz,gy,Tmass,teta,tp,sumu,sumv,CFS1)
use M_Math,                  only:CROSS12
use M_General_2D          ,  only:u,v,x,y,nx,ny
use M_Platform_Constant_2D,  only:r2,len,floatx
       
implicit none
integer nz,tp
double precision,intent (in)                              :: gy,dt
double precision,intent (in)   ,dimension(0:nx,0:ny)      :: ro
double precision,intent (in)   ,dimension(1:nx-1,1:ny-1)  :: CFS1
double precision,intent (out)                             :: Tmass
double precision,intent (inout)                           :: UbarV(2),teta,omegaz,xbarV(2),sumu,sumv


!! local variables !!

double precision dm,Izzb,sumIzb,AAA,BBB,CCC,mom(1:nx-1,1:ny-1),Txm,Tym
integer i,j,k
!!!!!!!!!!!!!!!!!!!!!!!

 
  


  Tmass=0 
  Txm=0 
  Tym=0
  do i=1,nx-1 ; do j=1,ny-1 

    dm=0.25*(x(i+1)-x(i-1))*(y(j+1)-y(j-1))*ro(i,j)*CFS1(i,j)
    Tmass=Tmass+dm
    Txm=Txm+x(i)*dm
    Tym=Tym+y(j)*dm
  end do ; end do

  XBarV(1)=Txm/Tmass 
  XBarV(2)=Tym/Tmass

  Tmass=0 
  sumu =0 
  sumv =0 
  sumIzb=0     !!!sumIx=Hx , SumIy=Hy, SumIz=Hz !!
  Izzb=0   
  mom(:,:)=0 
  do i=1,nx-1 ; do j=1,ny-1 

    dm=0.25*(x(i+1)-x(i-1))*(y(j+1)-y(j-1))*ro(i,j)*CFS1(i,j)
    Tmass=Tmass+dm
    Izzb=Izzb+ ( (x(i)-xbarV(1))**2+ (y(j)-XbarV(2))**2 )*dm   

    AAA=0.5*(u(i,j)+u(i+1,j))*dm 
    BBB=0.5*(v(i,j)+v(i,j+1))*dm   

    sumu=sumu+AAA 
    sumv=sumv+BBB 
    mom(i,j)=(x(i)-XBarV(1))*BBB -(y(j)-XBarV(2))*AAA 

    sumIzb=sumIzb+ mom(i,j)  
  
  end do ; end do 

print*,"distributed XBar=",Txm/Tmass, "distributed YBar=",Tym/Tmass
print*,"distributed mass=",Tmass,"Dsitributed mass moment of inertia",Izzb

!! calculating linear and angular velocity of the solid !! 

omegaz=   (sumIzb)/(Izzb)
UBarV(1) =(sumu)/Tmass +dt/Tmass*(-197.58*(xbarV(1)-floatx)-19.8*(UBarV(1)))
UBarV(2)  = (sumv)/Tmass


!! Correcting Navier stokes prediction for velocity field !!                                                                
  do i=2,nx-1 ; do j=2,ny-1 

    u(i,j)=u(i,j)+0.5d0*( CFS1(i-1,j)+CFS1(i,j) )*( ( UBarV(1)  -omegaz*(y(j)-XbarV(2)) )-u(i,j) ) !! mistake solved!!
    v(i,j)=v(i,j)+0.5d0*( CFS1(i,j-1)+CFS1(i,j) )*( ( UBarV(2)  +omegaz*(x(i)-XBarV(1)) )-v(i,j) )

  end do ; end do !! be careful if the solid is near boundry then boundry condition of velocity should change !!! 

  print*,"omega=",omegaz,"mass moment of inertia respect to Cg",Izzb

  
  

  

return
 End
 

!!!!!!!!!!!!! smoothed heaviside function !!!!!!!!!!!!!
function HV(phi,eps)

double precision phi,pi,eps


pi=3.141592654D0
if (phi.gt.eps) then
HV=1.d0
else if (phi.lt.-eps) then
HV=0.d0
else
HV=0.5d0*(  1.d0+ phi/eps +(1/pi)*dsin(pi*phi/eps)  )
end if 

return
end 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!subroutine for finding weather a point is inside the solid or not !!!



!!!!!!!!!!!!!!!!!! subroutine for boundary condition of  level set function !!!!!!!!!!!!!!!!!!!!
subroutine boundarycondphi()
use M_General_2D,  only:phi,nx,ny 
implicit none 
integer i,j



do j=1,ny-1 
phi(0,j)=2*phi(1,j)-phi(2,j)
phi(-1,j)=3*phi(1,j)-2*phi(2,j)

phi(nx,j)=2*phi(nx-1,j)-phi(nx-2,j)
phi(nx+1,j)=3*phi(nx-1,j)-2*phi(nx-2,j)
end do 

do i=-1,nx+1 
phi(i,0)=2*phi(i,1)-phi(i,2)
phi(i,-1)=3*phi(i,1)-2*phi(i,2)

phi(i,ny)=2*phi(i,ny-1)-phi(i,ny-2)
phi(i,ny+1)=3*phi(i,ny-1)-2*phi(i,ny-2)
end do 

return 
end 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!! finding tether forces ( not in the code now) !!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!! subroutine for adding tower and nacelle mass (not in the code now ) !!!!!!!!!

 



 
 
 

   
 subroutine advection (advectu,advectv)  

use M_General_2D,  only:u,v,x,y,nx,ny
 implicit none 
 double precision ux,uy,vx,vy
 double precision,dimension (1:nx,0:ny)       ::   dpux,dmux,dpuy,dmuy,advectu
 double precision,dimension (0:nx,1:ny)        ::   dpvx,dmvx,dpvy,dmvy,advectv

 
 integer i,j,k
 !!!!!!!!!!!! second order ENO method for advection terms !!!!!!!!!!!!
 
 do i=1,nx ; do j=0,ny 
 Dpux(i,j)=( u(i+1,j)-u(i,j)   )/( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )
 Dmux(i,j)=( u(i,j)  -u(i-1,j) )/( 0.5*(x(i)  +x(i-1))-0.5*(x(i-1)+x(i-2)) )
 Dpuy(i,j)=( u(i,j+1)-u(i,j)   )/( y(j+1)-y(j) )
 Dmuy(i,j)=( u(i,j)  -u(i,j-1) )/( y(j)-y(j-1) )
 
 end do ;end do 

 do i=2,nx-1 ; do j=1,ny-1    !!A!!
 
 
 !!U1    
 if (u(i,j).gt.0.d0) then

    if (  abs(  Dmux(i,j)-Dmux(i-1,j) ).lt.abs(  Dpux(i,j)-Dpux(i-1,j) )   ) then
 ux=Dmux(i,j)+  0.5*(  Dmux(i,j)-Dmux(i-1,j)  ) 
    else
 ux=Dmux(i,j)+  0.5*(  Dpux(i,j)-Dpux(i-1,j)  )
    end if 
    
 else

    if (  abs(  Dmux(i+1,j)-Dmux(i,j) ).lt.abs(  Dpux(i+1,j)-Dpux(i,j) )   ) then
 ux=Dpux(i,j)-  0.5*(  Dmux(i+1,j)-Dmux(i,j)  ) 
    else
 ux=Dpux(i,j)-  0.5*(  Dpux(i+1,j)-Dpux(i,j)  )
    end if 

 end if 


!!U2
 if (  (0.25*( V(i,j)+V(i,j+1)+v(i-1,j)+v(i-1,j+1) )).gt.0.d0) then

    if (  abs(  DmuY(i,j)-DmuY(i,j-1) ).lt.abs(  DpuY(i,j)-DpuY(i,j-1) )   ) then
 uY=DmuY(i,j)+  0.5*(  DmuY(i,j)-DmuY(i,j-1)  ) 
    else
 uY=DmuY(i,j)+  0.5*(  DpuY(i,j)-DpuY(i,j-1)  )
    end if 
    
 else

    if (  abs(  DmuY(i,j+1)-DmuY(i,j) ).lt.abs(  DpuY(i,j+1)-DpuY(i,j) )   ) then
 uY=DpuY(i,j)-  0.5*(  DmuY(i,j+1)-DmuY(i,j)  ) 
    else
 uY=DpuY(i,j)-  0.5*(  DpuY(i,j+1)-DpuY(i,j)  )
    end if 

 end if 
 
 
 
 advectu(i,j)=-(  u(i,j)*uX+0.25*( V(i,j)+V(i,j+1)+v(i-1,j)+v(i-1,j+1) )*uY  )

 end do ; end do
  
 do i=0,nx ; do j=1,ny 
 Dpvx(i,j)=( v(i+1,j)-v(i,j)   )/( x(i+1)-x(i) )
 Dmvx(i,j)=( v(i,j)  -v(i-1,j) )/( x(i)-x(i-1) )
 Dpvy(i,j)=( v(i,j+1)-v(i,j)   )/( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )
 Dmvy(i,j)=( v(i,j)  -v(i,j-1) )/( 0.5*(y(j)  +y(j-1))-0.5*(y(j-1)+y(j-2)) )
 
 end do ;end do 

 do i=1,nx-1 ; do j=2,ny-1  !!A!!
 
 
!!V1
 if (  (0.25*( u(i,j)+u(i+1,j)+u(i,j-1)+u(i+1,j-1) )).gt.0.d0) then

    if (  abs(  Dmvx(i,j)-Dmvx(i-1,j) ).lt.abs(  Dpvx(i,j)-Dpvx(i-1,j) )   ) then
 vx=Dmvx(i,j)+  0.5*(  Dmvx(i,j)-Dmvx(i-1,j)  ) 
    else
 vx=Dmvx(i,j)+  0.5*(  Dpvx(i,j)-Dpvx(i-1,j)  )
    end if 
    
 else

    if (  abs(  Dmvx(i+1,j)-Dmvx(i,j) ).lt.abs(  Dpvx(i+1,j)-Dpvx(i,j) )   ) then
 vx=Dpvx(i,j)-  0.5*(  Dmvx(i+1,j)-Dmvx(i,j)  ) 
    else
 vx=Dpvx(i,j)-  0.5*(  Dpvx(i+1,j)-Dpvx(i,j)  )
    end if 

 end if 



!!V2
 if ( V(i,j).gt.0.d0) then

    if (  abs(  DmvY(i,j)-DmvY(i,j-1) ).lt.abs(  DpvY(i,j)-DpvY(i,j-1) )   ) then
 vY=DmvY(i,j)+  0.5*(  DmvY(i,j)-DmvY(i,j-1)  ) 
    else
 vY=DmvY(i,j)+  0.5*(  DpvY(i,j)-DpvY(i,j-1)  )
    end if 
    
 else

    if (  abs(  DmvY(i,j+1)-DmvY(i,j) ).lt.abs(  DpvY(i,j+1)-DpvY(i,j) )   ) then
 vY=DpvY(i,j)-  0.5*(  DmvY(i,j+1)-DmvY(i,j)  ) 
    else
 vY=DpvY(i,j)-  0.5*(  DpvY(i,j+1)-DpvY(i,j)  )
    end if 

end if 
 
 
 advectv(i,j)=-(0.25*(  u(i,j)+u(i+1,j)+u(i,j-1)+u(i+1,j-1)  )*vX+ V(i,j)*vY)

 end do ; end do 
 
 
 return 
 
  end 


  



Subroutine viscosity (ro,miuv,Tx,Ty)
use M_General_2D,  only:u,v,x,y,nx,ny  
implicit none 
dOUBLE PRECISION  Txxr,Txxl,Tyxd,Tyxu,Tyyu,Tyyd,Txyr,Txyl
double precision,dimension (0:nx,0:ny)            ::   ro,miuv
double precision  Tx(2:nx-1,1:ny-1),Ty(1:nx-1,2:ny-1)
integer i,j

!!!!!!!!!!!! simple centeral difference for viscose terms !!!!!!!!!!!!!!!!!!

Do i=2,nx-1 ; Do j=1, ny-1 

 Txxr=2*miuv(i,j)  *( U(i+1,j)-U(i,j)   )/( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )
 
 Txxl=2*miuv(i-1,j)*( U(i,j)  -U(i-1,j) )/( 0.5*(x(i)  +x(i-1))-0.5*(x(i-1)+x(i-2)) )

 Tyxu=0.25*( miuv(i,j)+miuv(i,j+1)+miuv(i-1,j)+miuv(i-1,j+1) )*&
   &(      ( U(i,j+1)-U(i,j)     )/( Y(j+1)-y(j) )+ ( V(i,j+1)-V(i-1,j+1) )/( x(i)-x(i-1) )     )
                 
 Tyxd=0.25*( miuv(i,j-1)+miuv(i,j)+miuv(i-1,j-1)+miuv(i-1,j) )*&
   &(      ( U(i,j)-U(i,j-1)     )/( y(j)-y(j-1) )+ ( V(i,j)  -V(i-1,j)   )/( x(i)-x(i-1) )    )

  

 Tx(i,j)= (  (Txxr-Txxl)/(x(i)-x(i-1)) + (Tyxu-Tyxd)/( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )  )/(  0.5*(ro(i,j)+ro(i-1,j))  )
 
   
end do ; end do   




Do j=2,ny-1 ; Do i=1,nx-1 


Tyyu=2*miuv(i,j)   *( V(i,j+1)-V(i,j)   )/( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )
 
 Tyyd=2*miuv(i,j-1)*( V(i,j)  -V(i,j-1) )/( 0.5*(y(j)  +y(j-1))-0.5*(y(j-1)+y(j-2)) )
 
 
 Txyr=0.25*( miuv(i,j)+miuv(i,j-1)+miuv(i+1,j)+miuv(i+1,j-1) )*&
  & (     ( V(i+1,j)-V(i,j)     )/( x(i+1)-x(i) )+ ( U(i+1,j)  -U(i+1,j-1) )/( y(j)-y(j-1)  )   )
  
 Txyl=0.25*( miuv(i-1,j)+miuv(i-1,j-1)+miuv(i,j)+miuv(i,j-1) )*&
  & (     ( V(i,j)-V(i-1,j)     )/( x(i)-x(i-1) )+ ( U(i,j)    -U(i,j-1)   )/( y(j)-y(j-1)  )   )
  
 

 Ty(i,j)= (  (Tyyu-Tyyd)/( y(j)-y(j-1) ) + (Txyr-Txyl)/( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) ) )/(  0.5*(ro(i,j)+ro(i,j-1))  )

 
end do ; end do 
 

return 
end 


!!!!!!!!!!!!!! simple iterative method for solving poisson equation !!!!!!!!!!!!
subroutine poisson (ro,dt,pdif,p,beta)
use M_General_2D,  only:x,y,u,v,nx,ny 
implicit none 

double precision dt,beta,maxp,pdif
double precision,dimension (0:nx,0:ny)            ::   ro,miuv,p,pold
double precision,dimension (1:nx-1,1:ny-1)        ::   Apx,Amx,Apy,Amy,AP,Q
integer i,j,number2
 


do i=1,nx-1 ; do j=1,ny-1 
                         
 
 if (i==nx-1) then
 Apx(i,j)=0
 else
 

 Apx(i,j)=1/( (x(i+1)-x(i))*( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) ) )/(  0.5*(ro(i,j)+ro(i+1,j))  )
 !Apx(i,j)=1/( 0.5d0*(hx(i)+hx(i+1))*hx(i) )/(  0.5*(ro(i,j)+ro(i+1,j))  )
 end if 
 
 if (i==1)then 
 Amx(i,j)=0       
 else 
 Amx(i,j)=1/( (x(i)-x(i-1))*( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) ) )/(  0.5*(ro(i,j)+ro(i-1,j))  )
 !Amx(i,j)=1/( 0.5d0*(hx(i)+hx(i-1))*hx(i) )/(  0.5*(ro(i,j)+ro(i-1,j))  )
 end if 
 
 if (j==ny-1) then
 Apy(i,j)=0
 else
 Apy(i,j)=1/( (y(j+1)-y(j))*( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) ) )/(  0.5*(ro(i,j)+ro(i,j+1))  )
! Apy(i,j)=1/( 0.5d0*(hy(j)+hy(j+1))*hy(j) )/(  0.5*(ro(i,j)+ro(i,j+1))  )

 end if 
 if (j==1) then
 Amy(i,j)=0
 else
 Amy(i,j)=1/( (y(j)-y(j-1))*( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) ) )/(  0.5*(ro(i,j)+ro(i,j-1))  )
 !Amy(i,j)=1/( 0.5d0*(hy(j)+hy(j-1))*hy(j) )/(  0.5*(ro(i,j)+ro(i,j-1))  )
 
 end if 
 
 
 
 AP(i,j)=-( Apx(i,j)+Amx(i,j)+Apy(i,j)+Amy(i,j) )   
 
 !Q(I,J)=(   (U(I+1,J)-U(I,J))/hx(i)+(V(I,J+1)-V(I,J))/hy(j) )/dt
 Q(I,J)=(   (U(I+1,J)-U(I,J))/( 0.5*(x(i+1)+x(i)  )-0.5*(x(i)  +x(i-1)) )+&
          & (V(I,J+1)-V(I,J))/( 0.5*(y(j+1)+y(j)  )-0.5*(y(j)  +y(j-1)) )   )/dt 


 

 end do ; end do 
 
 
 number2=0
 maxp=1.5
           do while (maxp.gt.pdif.AND.number2.lt.5000)

 maxp=0 
 number2=number2+1
 
 do i=1,nx-1 ; do j=1,ny-1 

 pold(i,j)=p(i,j)
 p(i,j)=beta*( ( Apx(i,j)*P(I+1,j)+Amx(i,j)*P(I-1,j)+Apy(i,j)*P(I,j+1)+Amy(i,j)*P(I,j-1)-Q(I,j)  )/(-AP(i,j))  ) &
         & +(1-beta)*p(i,j)
 if (abs( pold(i,j)-p(i,j) ).gt.maxp) then 
 maxp=abs( pold(i,j)-p(i,j) )
 end if 
 
 end do ; end do
 
            end do  !!while !
                        
 
   print *,"pdif=",maxp,"number it=",number2
   
return 
end     
 
 
 
subroutine SolidUpdatePosition(npoints,point,UBarV,XBarV,omegaz,dt) 



implicit none
double precision,intent(inout), dimension(npoints+1,2)    :: point
double precision,intent(in)   , dimension(2)              :: XBarV,UBarV
double precision,intent(in)                               :: dt,omegaz
integer         ,intent(in)                               :: npoints


double precision,dimension(npoints+1,2)    :: pointTmp
double precision,dimension(2)              :: RTmp
integer i


  pointTmp(:,:)=point(:,:)
  do i=1,NPoints+1

    RTmp(:)=point(i,:)-XBarV(:)
    point(i,1)=pointTmp(i,1)+dt*(UbarV(1)-omegaz*(pointTmp(i,2)-xbarV(2)))
    point(i,2)=pointTmp(i,2)+dt*(UbarV(2)+omegaz*(pointTmp(i,1)-XBarV(1)))
    !point(i,:)=point(i,:)+dt*(UbarV(:)+CROSS12(omegaz,RTmp))
    
    ! point(i,:)=point(i,:)+dt*UbarV(:)
   end do 

return

end
 







 
  















