#include<vector>
#include "mkl.h"
#include "mkl_spblas.h"
#include "mkl_types.h"
#include "InputData.h"
#include "Sparse_Matrix.h"

#ifndef __FE_1D__
#define __FE_1D__

// Code for 1st Order FE Method
class FE_1D
{
private:
    /* data */
    int quadRule;                       // Quadrature Rule
    std::vector<double> quadPoints;     // Quadrature Points
    std::vector<double> quadWeight;     // Quadrature Weights
    std::vector<double> xiAtNode;       // xi at node
    std::vector<double> Nx;   // gradient of basis func evaluated at quad pts
    std::vector<double> N;  // Basis Function evaluated at quadrature Pts
    std::vector<double> detJ;


public:
    // Order of the Finite element Space 
    int FE_order;
    int N_Cells;
    int N_Nodes_Element;
    int N_DOF;
    std::vector<int> local2Global;
    //Sparse_Matrix* matrix;
    Sparse_Matrix* matrix;
    double h;
    std::vector<double> FGlobal;
    std::vector<double> Solution;

    // array for cells
    // Constructor
    FE_1D(int _order, int _N_cells, Sparse_Matrix* _matrix)
    {
        FE_order = _order;
        N_Nodes_Element = _order + 1; 
        N_Cells = _N_cells;
        N_DOF = (N_Nodes_Element-1)*N_Cells + 1;
        h = 1./double(N_DOF-1) ;
        matrix =  _matrix;
        // Initialise Nodal Co - ordinate Values 

      }
    

    // Function Declarations
    // To get the values of functions at the nodal points
    double basisFunctionValues(unsigned int , double);

    // To get the values of Gradient of basis functions at the nodal points
    double basisFunctionGradients(unsigned int , double);

    void  Initialise_xi_at_node();

    void Initialise_Local2Global();

    // Assembly Function 
    void assemble();

    // Assembly Function for symmetric Matrices
    void Smart_Assemble();

    //Quadrature Formulas 
    void Initialise_quadRules();

    //initialise matrix
    void Initialise_SparseMatrix(Sparse_Matrix* matrix);  

    //Initialise values of basis functions and basis gradients at the Quadrature Points
    void Initialise_values_at_Quad();

    // Initialise the determinant oof Jacobian 
    void Initialise_det_jacobian();

    // Setup all the initialisation functions
    void setup_FEsystem();

    // INitialise the Solution Vector and RHS Vector 
    void Initialise_sol_rhs();

    //Write the solution into an VTK File
    void Write_VTK(std::string name );

    // Write the Solution for Plotting the Solution
    void Plot_GNU();
    

        // Apply the boundary Condition 
    void apply_boundary_condition_normal(int a,double val1,int b,double val2);

    // Apply Boundary Condition while preserving Symmetricity
    void apply_boundary_condition_symmetric(int a , double val1,int b , double val2);

    void Assemble1D(int i){
        
        std::cout << "---------- Assembly Started --------------- " << std::endl;
        
        if(i ==1){
            std::cout << " Assembly Type : Normal" << std::endl;
            assemble();
        }
        else if (i ==2){
                std::cout << " Assembly Type : jugad" << std::endl;
            Smart_Assemble();
        }
        else
            std::cerr<<" Solver Type not implemented " << std::endl;
        
    }

    void apply_boundary_condition(int a,double val1,int b,double val2,bool i){
        if(i !=1 ){
            std::cout << " Boubdary Condition Type : Normal" << std::endl;
            apply_boundary_condition_normal(a,val1,b,val2);
        }
        else if (i == 1){
                std::cout << " Boubdary Condition Type : Symmetric " << std::endl;
            apply_boundary_condition_symmetric(a,val1,b,val2);
        }
        else
            std::cerr<<" Boundary Implementation Type not implemented " << std::endl;
    }


};

#endif

