module M_General_2D

     implicit none 
     save
integer, parameter :: nx=1201, ny=382 ! nx=601, ny=181 !nx =400, ny=120 !nx=460, ny=186 !nx=820, ny=373 !nx=512, ny=224 !nx=390,ny=224  !nx=355, ny=204 !nx=710, ny=408  !nx=374, ny=181 !nx=284  , ny=181  !ny=261    !241    !nx=35, ny=35 !nx=17 ,ny=17 !    
   
double precision, parameter  :: pi=3.141592653, gy=-9.81, froude=1.0, landa=6.00 
double precision, parameter  :: lx=22.0*Landa,  Ly=2.0*Landa, YFreeIni=Landa 
double precision,dimension (0:nx+1,-1:ny+1)       ::    u
double precision,dimension (-1:nx+1,0:ny+1)       ::    v
double precision,dimension (-1:nx+1,-1:ny+1)      ::   phi
real(8) x(-1:nx+1),y(-1:ny+1)
logical,parameter                                 :: test=.false.
       




end module
