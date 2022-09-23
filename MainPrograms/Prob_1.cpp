#include<iostream>
#include<vector>
#include<cmath>
#include<cstring>
#include<fstream>
#include <ctime>
#include <omp.h>
#include "InputData.h"
#include "Iterative_Solvers.h"
#include "Sparse_Matrix.h"
#include "FE_1D.h"
#include "Direct_Solvers.h"

//#include "solvers.h"

#ifdef _DEBUG
    #define DEBUG = _DEBUG
#endif

//using namespace std;

enum BOUNDARY_CONDITION {DIRICHLET,NEUMANN};

 

int main(int argc , char** argv)
{

    omp_set_num_threads(InputData::NUM_THREADS);
    if(argc < 1){
        std::cerr<<" Insufficient input arguments " << std::endl;
        std::cout << " Args : 1 = DAT FILE " <<std::endl;
    }

   //std::string filename = "file.txt";
    InputData* a = new InputData();
    a->get_input_values(argv[1]);
    // Set up Object for Sparse Matrix
    Sparse_Matrix* matrix = new Sparse_Matrix;

    // generate Mesh and create object for FE System
    FE_1D FE(InputData::FE_ORDER,InputData::N_CELLS,matrix);
    
    // Setup FE System 
    FE.setup_FEsystem();

    FE.Assemble1D(InputData::ASSEMBLY_TYPE);

    // Apply Boundary Condition
    FE.apply_boundary_condition(DIRICHLET,0,DIRICHLET,1,InputData::SYMM_BOUND_COND);
    
  
 #ifdef DEBUG
    for( int i =0 ; i<FE.N_DOF;i++)
    {
        for(int j = 0 ; j <FE.N_DOF; j++)
            std::cout<<matrix->getValues(i,j)<<"\t";
        std::cout<<std::endl;
    }
    for( int i=0; i<FE.N_DOF; i++)
        std::cout<<FE.FGlobal[i]<<std::endl;
  #endif
   
    // Call the Solver Routine
   // Solver_Iterative(Solver_type, matrix, &(FE.FGlobal[0]), &(FE.Solution[0]), 1e-4, 600000);
    if(InputData::SOLVER_TYPE == 1)
        Solver_Direct(matrix,FE.FGlobal.data(),FE.Solution.data());
    else
        Solver_Iterative(InputData::SOLVER_TYPE, matrix, &(FE.FGlobal[0]), &(FE.Solution[0]), InputData::TOLERANCE,InputData::MAX_ITER);

    
    // Write Solution into a VTK File
    //FE.Write_VTK("sol");



    // Plot Solution in GNU PLOT>
    if(InputData::PLOT_GNU)    FE.Plot_GNU();


    //Free memory
    matrix->deleteMatrix();
    delete matrix;
    return 0;
}

