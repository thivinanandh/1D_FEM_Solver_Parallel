#include <omp.h>
#include<iostream>
#include<fstream>
#include<cstring>
#include<string>

//enum BOUNDARY_CONDITION {DIRICHLET,NEUMANN};

#ifndef __INPUT_DATA__
#define __INPUT_DATA__

class InputData
{
    private:
       

    public:
    // DECLARE FILE OBJECT
    std::ifstream file;
    static int SOLVER_TYPE;    
    static int FE_ORDER;
    static int ASSEMBLY_TYPE;
    static int N_CELLS;
    static int MAX_ITER;
    static double TOLERANCE;
    static double BOUND_VAL1;
    static double BOUND_VAL2;
    static bool BOUND_COND1;
    static bool BOUND_COND2;
    static double RELAXATION_JACOBI;
    static double RELAXATION_SOR;
    static bool SYMM_BOUND_COND;
    static int CG_PRECONDITIONER;
    static bool PRINT_RESIDUAL;
    static bool PLOT_GNU;
    static std::string  MATFILE;
    static std::string  BFILE;
    static int RESIDUAL_DISPLAY;
    static int NUM_THREADS;
    static int RESTART_PARAMETER_FOM;
    static int MAX_ITERATION_FOM;

    static double CONVECTION_COEFF;
    static double ADVECTION_COEFF;
    static double SOURCE_COEFF;
    static double FORCE_VAL;

    

  
    
    
      // Constructor 
    InputData()
    {
       
    }

    void get_input_values(std::string filename);
    
};
#endif