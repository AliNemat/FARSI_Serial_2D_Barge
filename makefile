

NICE    = ifort
Gfort   = gfortran44
FCOMP   = mpif90
F77     = gfortran
CCOMP   = mpicc
CXXCOMP = mpicxx

#HYPRE_DIR = /usr/local/lib #HYPRE_DIR = /afs/crc.nd.edu/x86_64_linux/scilib/hypre/2.0.0/intel
#HYPRE_LIBS =  -L$(HYPRE_DIR)-lHYPRE #-lHYPRE_struct_ls -lHYPRE_struct_mv -lHYPRE_krylov -lHYPRE_utilities -lHYPRE_DistributedMatrix -lHYPRE_IJ_mv


#-----------------------------------------------------------------------------------------------

EXE = July16thL6VVF

OBJ =	      M_Math.o    M_General_2D.o       M_Platform_Constant_2D.o         M_Mesh_2D.o    M_SolidFinder_2D.o    M_Output_2D.o    Sep1st2012twod.o 
INC =	      M_Math.f90  M_General_2D.f90     M_Platform_Constant_2D.f90       M_Mesh_2D.f90  M_SolidFinder_2D.f90  M_Output_2D.f90  Sep1st2012twod.f90  

$(EXE): $(OBJ)
	$(NICE) $(OBJ)  -o $(EXE) -O3

$(OBJ): $(INC)
        
%.o:	%.f90 
	$(NICE) $(FOPTS) -c $<

clean:
	rm -f *.o *.mod *__genmod.f90

