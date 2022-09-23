
#include "InputData.h"

int InputData::SOLVER_TYPE =1;    
int InputData::FE_ORDER = 1;
int InputData::ASSEMBLY_TYPE= 1;
int InputData::N_CELLS= 10;
int InputData::MAX_ITER= 10000;
double InputData::TOLERANCE= 1e-6;
double InputData::BOUND_VAL1= 0;
double InputData::BOUND_VAL2= 1;
double InputData::RELAXATION_JACOBI= 1.98;
double InputData::RELAXATION_SOR= 1.98;
bool InputData::SYMM_BOUND_COND = 0;
int InputData::CG_PRECONDITIONER = 0;
bool InputData::PRINT_RESIDUAL = 0;
bool InputData::PLOT_GNU = 0;
std::string InputData::MATFILE;
std::string InputData::BFILE;
int InputData::RESIDUAL_DISPLAY;
int InputData::NUM_THREADS;
int InputData::RESTART_PARAMETER_FOM;
int InputData::MAX_ITERATION_FOM;


void InputData::get_input_values(std::string filename)
{
        file.open(filename);

        std::string line;
        std::string temp;
        std::string temp1;
        while(!file.eof())
        {   
            file >> temp;
            if(temp == "SOLVER_TYPE:"){
                file >> temp1;
                SOLVER_TYPE = stoi(temp1);
            }
            else if(temp == "FE_ORDER:"){
                file >> temp1;
                FE_ORDER = stoi(temp1);
            }
            else if(temp == "ASSEMBLY_TYPE:"){
                file >> temp1;
                ASSEMBLY_TYPE = stoi(temp1);
            }
            else if(temp == "N_CELLS:"){
                file >> temp1;
                N_CELLS = stoi(temp1);
            }    
            else if(temp == "MAX_ITER:"){
                file >> temp1;
                MAX_ITER = stod(temp1);
            }    
            else if(temp == "TOLERANCE:"){
                file >> temp1;
                TOLERANCE = stod(temp1);
            }
            else if(temp == "BOUND_VAL1:"){

                file >> temp1;
                BOUND_VAL1 = stod(temp1);
            }
            else if(temp == "BOUND_VAL2:"){
                file >> temp1;
                
                BOUND_VAL2 = stod(temp1);
            }
            else if(temp == "BOUND_COND1:"){
                file >> temp1;
                BOUND_VAL2 = stoi(temp1);
            }
            else if(temp == "BOUND_COND2:"){
                file >> temp1;
                BOUND_VAL2 = stoi(temp1);
            }
            else if(temp == "RELAXATION_JACOBI:"){
                file >> temp1;
                RELAXATION_JACOBI = stod(temp1);
            }
            else if(temp == "RELAXATION_SOR:"){
                file >> temp1;
                RELAXATION_SOR = stod(temp1);
            }

            else if(temp == "SYMM_BOUND_COND:"){
                file >> temp1;
                SYMM_BOUND_COND = stoi(temp1);
            }

            else if(temp == "CG_PRECONDITIONER:"){
                file >> temp1;
                CG_PRECONDITIONER = stoi(temp1);
            }

            else if(temp == "PRINT_RESIDUAL:"){
                file >> temp1;
                PRINT_RESIDUAL = stoi(temp1);
            }

            else if(temp == "PLOT_GNU:"){
                file >> temp1;
                PLOT_GNU = stoi(temp1);
            }

            else if(temp == "MATFILE:"){
                file >> temp1;
                MATFILE = temp1;
            }

            else if(temp == "BFILE:"){
                file >> temp1;
                BFILE = temp1;
            }

            else if(temp == "RESIDUAL_DISPLAY:"){
                file >> temp1;
                RESIDUAL_DISPLAY = stoi(temp1);
            }

            else if(temp == "NUM_THREADS:"){
                file >> temp1;
                NUM_THREADS = stoi(temp1);
            }
            else if(temp == "RESTART_PARAMETER_FOM:"){
                file >> temp1;
                RESTART_PARAMETER_FOM = stoi(temp1);
            }
            
            else if(temp == "MAX_ITERATION_FOM:"){
                file >> temp1;
                MAX_ITERATION_FOM = stoi(temp1);
            }

            else
            {
            }
            
                
        }

    }
