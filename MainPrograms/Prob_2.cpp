#include<iostream>
#include<cstring>
#include<vector>
#include<cmath>
#include<fstream>
#include <ctime>
#include <cstdlib>
#include "InputData.h"
#include "Iterative_Solvers.h"
#include "Sparse_Matrix.h"
#include "FE_1D.h"
#include "Direct_Solvers.h"

void getBValfromFile(std::string filename, int col , double* b)
{
    std::string tmp;
    int size;
    std::ifstream  file;
    file.open(filename);
    file >> tmp;
    size = stoi(tmp);
    if(size !=  col){
        std::cout << " The Matrix Dimension and the B vector dimension dows not match" <<std::endl;
    }

    for ( int i= 0 ; i < size ; i++){
        file >> tmp;
        b[i] = stod(tmp);
    }


    file.close();
}

int main(int argc , char** argv)
{
    if(argc < 1){
        std::cerr<<" Insufficient input arguments " << std::endl;
        std::cout << " Args : 1 = DAT FILE " <<std::endl;
    }

   //Get the Input Parameter values from DAT FILE 
    InputData* a = new InputData();
    a->get_input_values(argv[1]);

    Sparse_Matrix* matrix = new Sparse_Matrix;

    matrix->importValuesfromFile(InputData::MATFILE);

    std::cout << " Matrix Dimensions :  " << matrix->row <<std::endl;
    std::cout << " NNZ :  " << matrix->NNZSize <<std::endl;
    // setup Solution and B Array 
    std::vector<double> b(matrix->row,0);
    std::vector<double> sol(matrix->row,0);  

    getBValfromFile(InputData::BFILE,matrix->col , b.data());

       // Call the Solver Routine
   // Solver_Iterative(Solver_type, matrix, &(FE.FGlobal[0]), &(FE.Solution[0]), 1e-4, 600000);
    if(InputData::SOLVER_TYPE == 1)
        Solver_Direct(matrix,b.data(),sol.data());
    else
        Solver_Iterative(InputData::SOLVER_TYPE, matrix,b.data(), sol.data(), InputData::TOLERANCE,InputData::MAX_ITER);
    

    return 0;
}