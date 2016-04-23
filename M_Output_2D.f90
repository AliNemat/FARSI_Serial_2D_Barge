Module M_output_2D
 use M_General_2D, only: u,v,nx,ny
 use M_Mesh_2D,    only: x,y

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Output_vtk_ldos - Plot Density of State (2d plot)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains 

subroutine Output_vtk(tp,ro,CFS1) !pass A and G, c is a counter

double precision,dimension (0:nx,0:ny),intent(in)            ::   ro
double precision,dimension (1:nx-1,1:ny-1),intent(in)        ::   CFS1
double precision,dimension (2)                               ::   z 
integer,intent (in)                                          ::   tp

 integer(kind=4)   :: k,j,i,nz
 !CHARACTER(len=5)  :: NNUMBER 
 CHARACTER(len=1024) :: FileName

nz=2

z(1)=1
!z(2)=1  
    !real(kind=8) :: w, l, Em
    !character*100 :: var
    !character*80 :: Fname

    !integer(kind=4), intent(in) :: c
    !real(kind=8), dimension(:,:), intent(in) :: A,G

    103 format (A)
    104 format (ES12.5)
    105 format (A,1X,I5,I5,I5)
    106 format (A,I5,A)
    107 format (A,I12)

    !this routine just make a filename based on the counter number (c)
    !Fname = Utility_MakeFileName('ldos',c,'vtk')
    !write(var,*) Fname
    !call Errors_Fileopen(err,var)

    write (FileName, fmt='(A6,I5,A4)') "result",tp+10000,'.vtk' 
    OPEN(unit=tp+10000,file=FileName)
!    OPEN(unit=10000,file='results.vtk')
    !OPEN(unit=tp+10000,file=FileName)
!    print*,FileName,"FileName"
!     if (1.eq.2) then 
    !Output header to file
    write(tp+10000,103) '# vtk DataFile Version 2.0' !don't mode
    write(tp+10000,103) 'Result for paraview 2d code' !Change to the title you want
    write(tp+10000,103) 'ASCII' !tell vtk that the points are ascii not binary
    write(tp+10000,103) 'DATASET RECTILINEAR_GRID' !type of element
    write(tp+10000,105) 'DIMENSIONS',Nx-1,Ny-1,Nz-1 !Number of points in each of the 3 directions

  !  l = 0; w = 0
     
    !because we are using a rectilinear grid only have to specify grid discretization
    ! along boarders of the domain
    write(tp+10000,106) 'X_COORDINATES ',Nx-1,' float'
    do i=1,Nx-1
      write(tp+10000,104) x(i) !Points
      !l = l + Gdx(n)
    end do

    write(tp+10000,106) 'Y_COORDINATES ',Ny-1,' float'
    do j=1,Ny-1
      write(tp+10000,104) y(j) !Points
      !w = w + Gdy(j)
    end do

    !Em = maxval(E)
    write(tp+10000,106) 'Z_COORDINATES ',Nz-1,' float'
    do k=1,Nz-1
      write(tp+10000,104) z(k) !E(r)/Em*l !Points
    end do

    ! Writing data about cells to the vtk file
    ! writing 3 different sets of data
    write(tp+10000,107) 'POINT_DATA', (Nx-1)*(Ny-1)*(Nz-1) !Number of elements


    write(tp+10000,103) 'SCALARS X_Velocity float 1' !change variable name
    write(tp+10000,103) 'LOOKUP_TABLE default'
    do k=1,Nz-1
       do j=1,Ny-1
         do i=1,Nx-1
         !j = (k-1)*Ndx + n
           write(tp+10000,104) 0.5*(u(i,j)+u(i+1,j)) !A(r,j)/(2.0*pi)!LDOS Contour Plot
         end do
       end do
    end do

    write(tp+10000,103) 'SCALARS V_Velocity float 1' 
    write(tp+10000,103) 'LOOKUP_TABLE default'
    do k=1,Nz-1
       do j=1,Ny-1
         do i=1,Nx-1
         !j = (k-1)*Ndx + n
           write(tp+10000,104) 0.5*(v(i,j)+v(i,j+1)) !G(r,j)/(2.0*pi) !LDOS Contour Plot
         end do
       end do
    end do

    write(tp+10000,103) 'SCALARS Density float 1'
    write(tp+10000,103) 'LOOKUP_TABLE default'
    do k=1,Nz-1
       do j=1,Ny-1
         do i=1,Nx-1
         !j = (k-1)*Ndx + n
           write(tp+10000,104) ro(i,j)  !(A(r,j)-G(r,j))/(2.0*pi) !LDOS Contour Plot
         end do
       end do
    end do

    write(tp+10000,103) 'SCALARS Marker_F float 1'
    write(tp+10000,103) 'LOOKUP_TABLE default'
    do k=1,Nz-1
       do j=1,Ny-1
         do i=1,Nx-1
         !j = (k-1)*Ndx + n
           write(tp+10000,104) CFS1(i,j)  !(A(r,j)-G(r,j))/(2.0*pi) !LDOS Contour Plot
         end do
       end do
    end do
!end if 
!   close(unit=tp+10000)! close file
end subroutine Output_vtk

end module M_output_2D


