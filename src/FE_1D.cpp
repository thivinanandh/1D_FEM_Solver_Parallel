#include "FE_1D.h"
#include "InputData.h"
#include<iostream>
#include<fstream>

void FE_1D::setup_FEsystem()
{
    std::cout << "-------------- Meshing ----------------" << std::endl;
    Initialise_xi_at_node();
    Initialise_Local2Global();       // Initialise the Local to global Mapping vector
    std::cout<<"No of Nodes : " << N_Nodes_Element<<std::endl;
    std::cout<<"FE order    : " << FE_order<<std::endl;
    std::cout<<"No of Cells : " << N_Cells<<std::endl;
    std::cout << "-------------- Meshing [Completed] ----------------" << std::endl; 
    std::cout << "-------------- Setting Up FE System ----------------" << std::endl;
    Initialise_quadRules();          // Initialise quadrature Rules
    Initialise_SparseMatrix(matrix); // Initialise object for the matrix
    Initialise_values_at_Quad();
    Initialise_det_jacobian();
    Initialise_sol_rhs();
    std::cout << "-------------- Setting Up FE System [ COMPLETED ]----------------" << std::endl;
    std::cout<<" N_DOF    : " << N_DOF<<std::endl;
    std::cout<<" NNZ      : " << matrix->NNZSize<<std::endl;
    std::cout<<" Sparsity : " <<double( 1 - (double(matrix->NNZSize)/double(N_DOF*N_DOF)))<<std::endl;

}

void FE_1D::Initialise_xi_at_node()
{
    xiAtNode.resize(N_Nodes_Element,0);

    xiAtNode[0] = -1;
    xiAtNode[1] = 1;

    if(N_Nodes_Element> 2)
        for ( int i = 2; i<N_Nodes_Element ; i++){
            xiAtNode[i] = -1 + (2./(FE_order))*(i-1);
        }
}

// Provides the Values of Basis Function Values
double FE_1D::basisFunctionValues(unsigned int node, double xi)
{
    double value = 1.;

    for ( int i = 0 ; i< N_Nodes_Element ; i++)
        if ( i != node)
            value *= (xi - xiAtNode[i] ) / (xiAtNode[node] - xiAtNode[i]) ; 
    return value;
}

//Provides the Values of Basis Function Gradients
double FE_1D::basisFunctionGradients(unsigned int j, double xi)
{
    // double value = 0;

    // for ( int i = 0 ; i< N_Nodes_Element ; i++)
    //     if ( i != node)
    //         value += 1. / (xi - xiAtNode[i]) ; 
    
    // value *= basisFunctionValues(node,xi);

    double value = 0.0;
    for ( int i = 0 ; i< N_Nodes_Element ; i++){
        if(i == j) continue;
        double dummy = 1.0/(xiAtNode[j] - xiAtNode[i]);
        for(int m = 0; m < N_Nodes_Element; m++){
            if(m == i || m == j) continue;
            dummy *= (xi - xiAtNode[m])/(xiAtNode[j] - xiAtNode[m]);
        }
        value += dummy;
    }

    return value;
}


// Creates and sets the GLOBAL to local DOF Mappings
void FE_1D::Initialise_Local2Global()
{
    int size =  N_Cells*N_Nodes_Element;
    int cell;
    local2Global.resize(size,0);
    for ( int i = 0  ; i< size ; i++){
        cell = i/N_Nodes_Element;
        if( i % N_Nodes_Element == 0) // 0th local node in each element - corresponds to start Global DOF
            local2Global[i] = cell*FE_order;
        else if (i % N_Nodes_Element == 1) // 1st local node in each element - corresponds to end globl DOF
            local2Global[i] = cell*FE_order + FE_order; 
        else
            local2Global[i] = cell*FE_order + (i%N_Nodes_Element)-1;  
    }
    //std::cout<<" Meshing Completed" << std::endl;

}

void FE_1D::Initialise_quadRules()
{
    quadRule = N_Nodes_Element;

   if(quadRule == 2)
   {    
       quadWeight.resize(2,0);
       quadWeight[0] = 1.0000000000000000;
       quadWeight[1] = 1.0000000000000000;
       
       quadPoints.resize(2,0);
       quadPoints[0] = -0.5773502691896257;
       quadPoints[1] = 0.5773502691896257;
   }
   
   else if(quadRule == 3)
    {
        quadWeight.resize(3,0);
        quadWeight[0] = 0.8888888888888888;
        quadWeight[1] = 0.5555555555555556;
        quadWeight[2] = 0.5555555555555556;

        quadPoints.resize(3,0);
        quadPoints[0] = 0.0000000000000000;
        quadPoints[1] = -0.7745966692414834;
        quadPoints[2] = 0.7745966692414834;
    }

    else if(quadRule==4)
    {
        quadWeight.resize(4,0);
        quadWeight[0] = 0.6521451548625461;
        quadWeight[1] = 0.6521451548625461;
        quadWeight[2] = 0.3478548451374538;
        quadWeight[3] = 0.3478548451374538;

        quadPoints.resize(4,0);
        quadPoints[0] = -0.3399810435848563;
        quadPoints[1] = 0.3399810435848563;
        quadPoints[2] = -0.8611363115940526;
        quadPoints[3] = 0.8611363115940526;
    }

    else if(quadRule = 5)
    {
        quadWeight.resize(5,0);
        quadWeight[0] = 0.5688888888888889;
        quadWeight[1] = 0.4786286704993665;
        quadWeight[2] = 0.4786286704993665;
        quadWeight[3] = 0.2369268850561891;
        quadWeight[4] = 0.2369268850561891;

        quadPoints.resize(5,0);
        quadPoints[0] = 0.0000000000000000;
        quadPoints[1] = -0.5384693101056831;
        quadPoints[2] = 0.5384693101056831;
        quadPoints[3] = -0.9061798459386640;
        quadPoints[4] = 0.9061798459386640;
    }

    else
    {
        std::cout<<" The QUADRATURE FORMULA IS NOT  SET FOR FE_ORDER " << FE_order << std::endl;
        exit(0);     
    }
    std::cout<<" The Quadrature Rule Selected is " << quadRule << std::endl;

}


void FE_1D::Initialise_SparseMatrix(Sparse_Matrix* matrix)
{
    matrix->setSparsityPattern(N_Nodes_Element,N_Cells,local2Global);
    std::cout<<" Matrix Initialisation Completed "<< std::endl;
}

void FE_1D::Initialise_values_at_Quad()
{
    // Calculate values of shape functions and shape functions gradient
    N.resize(quadRule*N_Nodes_Element,0);
    Nx.resize(quadRule*N_Nodes_Element,0);
   for(int i = 0 ; i<N_Nodes_Element ;i++)
        for ( int q = 0 ; q<quadRule ; q++) {
            N[i*quadRule + q]   = basisFunctionValues(i,quadPoints[q]);
            Nx[i*quadRule + q]  = basisFunctionGradients(i,quadPoints[q]);  
    }
    
    std::cout << " Values initiated at quadrature points"<<std::endl;
  
}

void FE_1D::Initialise_det_jacobian()
{
    detJ.resize(quadRule,0);
    for (int i = 0 ; i < quadRule ; i++){
        detJ[i] = quadWeight[i] * 2/(h);     // where h * 0.5 is the det of Jacobian matrix for an 1D Element 
    }
    std::cout << " Values initiated for Det J"<<std::endl;
}
   
void FE_1D::Initialise_sol_rhs()
{
    FGlobal.resize(N_DOF,0);
    Solution.resize(N_DOF,0);
}



void FE_1D::apply_boundary_condition_normal(int a,double Boundval1 , int b, double Boundval2)
{
     if(a + b > 1)
    {
        std::cout<< "Cannot provide boundary condition as Neumann on both the edges for a stationary problem " << std::endl;
        exit(0);
    }
   
    if(a == 0) 
    {
        int end = matrix->rowPtr[1];
        for( int i = 0 ; i < end; i++){
            int j = matrix->colPtr[i];
            matrix->getValues(0,j) = 0.;
            matrix->getValues(0,0) = 1.;
        }
        FGlobal[0] = Boundval1;
    }
    
    if(b == 0){
        int end = matrix->rowPtr[N_DOF ];
        int start = matrix->rowPtr[N_DOF -1];
        for ( int i = start ; i<end; i++)
        {
            int j = matrix->colPtr[i];
            matrix->getValues(N_DOF-1,j) = 0.;
            matrix->getValues(N_DOF-1,N_DOF-1) = 1.;
        }
        FGlobal[N_DOF-1] = Boundval2;
    }
    
    if(a == 1)
        FGlobal[0] += Boundval1;
    
    if(b == 1)
        FGlobal[N_DOF-1] +=Boundval2;
}

void FE_1D::Smart_Assemble()
{   
   	mkl_set_num_threads(InputData::NUM_THREADS);
	mkl_set_dynamic(false);
   
    // Declare a 2D Array - for Local Stiffness
    std::vector<double> flocal(N_Nodes_Element);
    std::vector<std::vector<double>> klocal(N_Nodes_Element, std::vector<double>(N_Nodes_Element,0));

    // Initialise Global F Vector and the Global Solution Vector
    
    std::fill(FGlobal.begin(), FGlobal.end(), 0);
    
    std::fill(Solution.begin(), Solution.end(), 0);
    double val = 0;

    auto start = clock();

    for ( int i = 0 ; i < N_Nodes_Element ; i++)
            for(int j=0 ; j < N_Nodes_Element; j++)
                for(int q = 0 ; q < quadRule ; q++)
                    klocal[i][j] += Nx[i*quadRule + q]*Nx[j*quadRule + q]*detJ[q];// Fill the Values here; 

    
    // Loop Over cells
    for(int cell = 0; cell<N_Cells;cell++)
    {   
        int cell_index = cell*N_Nodes_Element;
        // Loop over all Quadrature Points 
            for ( int i = 0 ; i < N_Nodes_Element ; i++)
            {   
                flocal[i] = 0;
                FGlobal[local2Global[cell_index + i]] += flocal[i];
                for(int j=0 ; j < N_Nodes_Element; j++)
                    matrix->getValues(local2Global[cell_index+i],local2Global[cell_index+j]) += klocal[i][j];
                flocal[i] = 0;
            }
    }
    
    auto stop = clock(); 
    auto duration = (stop -start)/double(CLOCKS_PER_SEC)*1000; 
    std::cout << " Time Taken for assembly : " << duration <<" ms" << std::endl;
    std::cout << " ----------------- Assembly [Completed] -------------------- " << std::endl;

    // std::cout<<std::endl;
    // for( int i =0 ; i<N_DOF;i++)
    // {
    //     for(int j = 0 ; j <N_DOF; j++)
    //         std::cout<<matrix->getValues(i,j)<<"\t";
    //     std::cout<<std::endl;
    // }

    
}



void FE_1D::assemble()
{   
     mkl_set_num_threads(InputData::NUM_THREADS);
	mkl_set_dynamic(false);
      auto start = clock(); 
    // Declare a 2D Array - for Local Stiffness
    std::vector<double> flocal(N_Nodes_Element);

    // Initialise Global F Vector and the Global Solution Vector
     std::fill(FGlobal.begin(), FGlobal.end(), 0);
    
    std::fill(Solution.begin(), Solution.end(), 0);
    double val = 0;
    double fval = 0;
    int global_i;
    int global_j;

    
    // Loop Over cells
    // #pragma omp parallel for
    for(int cell = 0; cell<N_Cells;cell++)
    {   
        int cell_index = cell*N_Nodes_Element;
        // Loop over all Nodal Points - i // for test

        for(int q = 0 ; q < quadRule ; q++)
        {
            for ( int i = 0 ; i < N_Nodes_Element ; i++)
            {
                global_i = local2Global[cell_index + i];
                
                FGlobal[global_i] += InputData::FORCE_VAL * N[i*quadRule + q]*detJ[q];
                // FGlobal[global_i] = 1;
                for(int j=0 ; j < N_Nodes_Element; j++)
                {   
                    global_j = local2Global[cell_index + j];
                    matrix->getValues(global_i,global_j) += ( (InputData::CONVECTION_COEFF* Nx[i*quadRule + q]*Nx[j*quadRule + q] ) + (InputData::ADVECTION_COEFF * Nx[j*quadRule + q] * N[i*quadRule * q] ) ) * detJ[q];
                }
            }
        }


        
    }
 
    auto stop = clock(); 
    auto duration = (stop -start)/double(CLOCKS_PER_SEC)*1000; 
    std::cout << " Time Taken for assembly : " << duration <<" ms" << std::endl;
    std::cout << " ----------------- Assembly [Completed] -------------------- " << std::endl;
}


// Applies the dirichlet Boubdary Condition while preserving Symmetricity
// 0 - Dirichlet
// 1 - Neumann
void FE_1D::apply_boundary_condition_symmetric(int a , double Boundval1,int b , double Boundval2)
{
      if(a + b > 1)
    {
        std::cout<< "Cannot provide boundary condition as Neumann on both the edges for a stationary problem " << std::endl;
        exit(0);
    }
   
    if(a == 0) //  DIRICHLET
    {
        int row = 0;
        int end = matrix->rowPtr[1];
        for( int i = 0 ; i < end; i++){
            int j =  matrix->colPtr[i];
            if(row == j)
            {
                matrix->getValues(row,j) = 1;
                FGlobal[j] = Boundval1;
            }
            else
            {
                matrix->getValues(row,j) = 0.;
                FGlobal[j] -= matrix->getValues(j,row)*Boundval1;
                matrix->getValues(j,row) = 0.;
            }
        }
    }     




    if(b == 0){
        int end = matrix->rowPtr[N_DOF ];
        int start = matrix->rowPtr[N_DOF -1];
        int row = N_DOF-1;
        for( int i = 0 ; i < end; i++){
            int j =  matrix->colPtr[i];
            if(row == j)
            {
                matrix->getValues(row,j) = 1;
                FGlobal[j] = Boundval2;
            }
            else
            {
                matrix->getValues(row,j) = 0.;
                FGlobal[j] -= matrix->getValues(j,row)*Boundval2;
                matrix->getValues(j,row) = 0.;
            }
        }
    }
    
    if(a == 1)
        FGlobal[0] += Boundval1;
    
    if(b == 1)
        FGlobal[N_DOF-1] +=Boundval2;
}



void FE_1D::Plot_GNU()
{
    std::ofstream file;
    file.open("GNUPLOT_FILE.dat");
    
    //file.flags( std::ios::dec | std::ios::scientific);
    file.precision(6);

    file<< "#"<<"   "<<"X"<<"   "<<"Y"<<std::endl;

    for ( int i = 0 ; i < N_DOF ; i++)
    {
        file << i*h << "   " << Solution[i]<<std::endl;
    }

    file.close();

    system("gnuplot -e \"plot 'GNUPLOT_FILE.dat' using 1:2 with linespoints; pause -1 ; exit\" ");
}



void FE_1D::Write_VTK(std::string filename)
{
    filename.append(".vtk") ;
    std::ofstream myfile(filename);

    myfile.flags( std::ios::dec | std::ios::scientific);
    myfile.precision(6);
 
    // LIne 1 - VTK Version
    myfile << "# vtk DataFile Version 1.0"<<std::endl;

    //Line 2 - Title
    myfile <<" Solution 1D " << std::endl<<std::endl;

    //Line 3 - Data type ( ASCII / BINARY )
    myfile <<"ASCII" << std::endl;

    //Line 3 - Structured or unstructured grid ( STRUCTURED_GRID / UNSTRUCTURED_GRID )
    myfile <<"DATASET UNSTRUCTURED_GRID" << std::endl;

    //Line 4 - no of data points :: syntax -> POINTS <no of points> <datatype>
    myfile <<"POINTS " <<N_DOF<<" float" <<std::endl;

    for( int i = 0; i < N_DOF; i++)
        myfile <<i*(h)<<" "<<"0.0000000 "<<"0.000000"<<std::endl;
    
    myfile<<std::endl;

    // CELLS -> syntax : CELLS <no of cells> <no of parameters  totally needed to define the total cell points>
    // for eg: in 1D total points in cells is 2 , So the last parameter will be (2+1)* No of cells

    myfile<<"CELLS "  <<N_Cells<< " " << (N_Nodes_Element+1)*N_Cells<<std::endl;

    int node = 0;
    for( int i = 0; i < N_Cells; i++)
        myfile <<"2 "<<node<<" " <<++node<<std::endl;

    myfile<<std::endl;

    // CELL TYPES : syntax: CELL_TYPES <No of cells>
    // cell type for all the cells, 
    // for 1d Element - cell type is 3

    myfile<<"CELL_TYPES "  <<N_Cells<<std::endl;

    for(int i = 0 ; i<N_Cells;i++)
        myfile<<"3 ";
    
    myfile<<std::endl<<std::endl<<std::endl;


    // POINT DATA : syntax - PONT_DATA <no of points>
    // < scalar or vector> < Datatype>
    //"LOOKUPTABLE" < lookuptable type >

    myfile<<"POINT_DATA "<<N_DOF<<std::endl;
    myfile<<"SCALARS "<<"1D_Solution "<<"float"<<std::endl;
    myfile<<"LOOKUP_TABLE "<<"default" <<std::endl;

    for( int i =0 ; i < N_DOF ; i++)
        myfile<<Solution[i]<<std::endl;

    myfile.close();

    std::cout<<"VTK File "<< filename << " has been generated"<<std::endl;
}