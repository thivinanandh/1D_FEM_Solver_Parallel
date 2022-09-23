#include<cstring>
#include<cmath>
#include<cstdlib>
#include <fstream>
#include<iostream>
#include<vector>
#include "MATRIX_1D.h"
// Set the value of pi
double const PI = 4.0*atan(1.0);

using namespace std;

/*---------- End of user input parameters ---------------------------------*/

void memset_2d( double**& arr, int row , int col)
{
    arr = new double*[row] ;
    for (int i = 0 ; i < row ; i ++ )
    {
        arr[i] = new double[col];
        for (int j = 0 ; j< col ; j++ )
            arr[i][j] = 0.0;
    }
}



void memset_1d( double*& arr, int row)
{
    arr = new double[row] ;
    for (int i = 0 ; i < row ; i ++ )
    {
        arr[i] = 0.0;
    }
}

double vector_inner_dot(double* vec1 , double* vec2,int Nx)
{
    double ans = 0.0;
    for(int i = 0 ; i< Nx ; i++)
        ans += vec1[i]*vec2[i];
    return ans;
}

void sp_matrix_vector_multiply(double*& ans_arr,Matrix_1D* arr, double* vec, int size)
{
    for ( int i =0 ; i < size ; i++){
        ans_arr[i] = 0;
        for( int j =0 ; j< size ; j++){
            ans_arr[i] += arr->getValues(i,j) * vec[j] ; }
    }
}

void matrix_vector_multiply(double*& ans_arr,double** arr, double* vec, int size)
{
    for ( int i =0 ; i < size ; i++){
        ans_arr[i] = 0;
        for( int j =0 ; j< size ; j++)
            ans_arr[i] += arr[i][j] * vec[j] ; 
    }
}


void print_matrix(double** arr1,int row , int col , string var_name)
{
    
    std::cout<< var_name << " : " ;
    
    for(int i=0;i<row;i++)
    {
        cout<<endl;
        for (int j=0;j<col;j++){
            cout<<arr1[i][j]<<"\t";
        }
    }
    cout<<endl<<endl;;
}

void print_array(double* arr1,int row, string var_name)
{
    std::cout<< var_name << " : " ;
    for(int i=0;i<row;i++)
            cout<<arr1[i]<<"\t";
    cout<<endl<<endl;;
}

void matrix_multiply( double**& ans_arr, double** arr_a , double** arr_b, int row1, int col1,int col2, int i_dev, int j_dev)
{
    int i,j,k;
    for( i=0 + i_dev ;i < row1-i_dev;i++)
        for ( k=0 + j_dev ; k<col1 - j_dev ;k++)
            for ( j=0 + j_dev ;j < col2 - j_dev ;j++)    
                ans_arr[i-i_dev][j-j_dev]  += (arr_a[i][k] * arr_b[k][j]);
}

double norm_2_matrix(double** b , int n)
{
    double norm = 0.0;
    for (int i =0 ; i< n ; i ++ )
        for(int j = 0 ; j < n ; j++)
            norm += (b[i][j]*b[i][j])/(n*n);
        
    return sqrt(norm);
    
}

double norm_2_vector(double* b , int n)
{
    double norm = 0.0;
    for (int i =0 ; i< n ; i ++ )
            norm += (b[i]*b[i])/(n);
        
    return sqrt(norm);
    
}

double** BiCGSTAB_Matrix(double **A, double **b,double **x,int row, int col , double tolerance , int max_iter)
{
    /* Algorithm taken from SIAM Paper :
    BI-CGSTAB: A FAST AND SMOOTHLY CONVERGING VARIANT OF BI-CG FOR THE SOLUTION OF NONSYMMETRIC LINEAR SYSTEMS* H. A. VAN DER VORSTt 
    2d matrrix details
    A = Input matrix - Square and Real
    b = Result matrix - suare and Real
    x =  intial guess matrix
    tolerance =  error tolerance for convergence check
    max_iter = maxximum Iteration;
    
    Algorithm used:
        r0 = b − Ax0
        Choose an arbitrary vector r̂0 such that (r̂0, r0) ≠ 0, e.g., r̂0 = r0 . Note that notation (x,y) applies for scalar product of vectors (x,y) = <x,y> = x·y = x ' y
        ρ0 = α = ω0 = 1
        v0 = p0 = 0
        For i = 1, 2, 3, …
        ρi = (r̂0, ri−1)
        β = (ρi/ρi−1)(α/ωi−1)
        pi = ri−1 + β(pi−1 − ωi−1vi−1)
        vi = Api
        α = ρi/(r̂0, vi)
        h = xi−1 + αpi
        If h is accurate enough, then set xi = h and quit
        s = ri−1 − αvi
        t = As
        ωi = (t, s)/(t, t)
        xi = h + ωis
        If xi is accurate enough, then quit
        ri = s − ωit
    */
    
    int i,j;
    int i_deviation,j_deviation;
    clock_t BiCGSTAB_start, BiCGSTAB_end;
    // Declaration 
    double **r0,**r0_hat;
    double **p0,**p0_new,**v,**v_new, **H,**x_new, **S, **T, **Ax;
    double rho,alpha,omega, rho_new,omega_new,error,beta;
    
    memset_2d(r0,row,col);memset_2d(r0_hat,row,col);
    memset_2d(p0,row,col);memset_2d(p0_new,row,col);memset_2d(v,row,col);memset_2d(v_new,row,col);
    memset_2d(H,row,col);memset_2d(x_new,row,col);memset_2d(S,row,col);memset_2d(T,row,col);memset_2d(Ax,row,col);
    
    double B_norm = norm_2_matrix(b,row);
    
    
    // Stp -1 ; Calulate R0 , Assign it to R0_new
    // Calculate Ax
     matrix_multiply(Ax,A,x,row,col,col,0,0);


    for (i = 0; i< row ; i++){
        for(j=0;j< col ; j++){
            r0[i][j] = b[i][j] - Ax[i][j];
            r0_hat[i][j] = r0[i][j];
        }
    }
    
    // decalatre values of Constants - Alpha , omega 
    alpha = rho = omega = rho= 1.0;
    
    // Set values of v0 and p0 as 0 Matrix, this would have been  done during declaration itself. 
    BiCGSTAB_start = clock();
    int iteration_count = 1;
    double L2_error = 1000;
    while ( (fabs(L2_error) > tolerance) && (iteration_count  <= max_iter))
    {
        //ρi = (r̂0, ri−1)
        rho_new = 0.0;
        for (i = 0;i<row;i++)
            for (j=0;j<col;j++)
                rho_new += r0_hat[i][j] * r0[i][j];
        
         if( iteration_count  > 1 )
         {
            // β = (ρi/ρi−1)(α/ωi−1)
            beta = (rho_new/rho) * (alpha/omega);
            //Pi = ri−1 + β(Pi−1 − ωi−1Vi−1)
            for (i = 0;i<row;i++)
                for (j=0;j<col;j++)
                    p0_new[i][j] = r0[i][j] + beta*(p0[i][j] -  (omega*v[i][j])) ;
         }
         else
             p0_new  = r0;
        
        // vi = APi
        matrix_multiply(v_new,A,p0_new,row,col,row,0,0);
               
        //α = ρi/(r̂0, vi)clc
        double temp  = 0.0;
        for (i = 0;i<row;i++)
            for (j=0;j<col;j++)
                temp += r0_hat[i][j] * v_new[i][j];
        alpha = rho_new/temp;
        
        for (i = 0;i<row;i++)
            for (j=0;j<col;j++)
                H[i][j] = x[i][j] + alpha*p0_new[i][j];
        
        // Early Exit Criteria , if the Norm of (alpha * p0 ) is less than tolerance , then exit the loop 
        L2_error = 0.0;
        L2_error = (alpha*norm_2_matrix(p0_new,row)/B_norm);
        if( (iteration_count > 1) && ( fabs(L2_error) < tolerance))
        {
            BiCGSTAB_end = clock();
            cout << "********* BICGSTAB  has converged ********** " << endl;
            cout<<"Time Elapsed : "<< double (BiCGSTAB_end - BiCGSTAB_start)/ CLOCKS_PER_SEC<<"  ms "<< endl;
            cout << "Iteration :  " << iteration_count << ".5 " <<endl<<"Error :" <<  L2_error << endl;
            cout<<"Tolerance Provided : " << tolerance <<endl;
            //print_matrix(H,row,col,getName(h));
            return H;
        }

        // s = ri−1 − αvi
        for (i = 0;i<row;i++)
            for (j=0;j<col;j++)
                S[i][j] = r0[i][j] -  alpha*v_new[i][j];
        //print_matrix(S,row,col,getName(S));
        // t = As
        matrix_multiply(T,A,S,row,col,col,0,0);
        //print_matrix(T,row,col,getName(T));
        
        //ωi = (t, s)/(t, t)
        double temp2 = 0.0;
        double temp1 = 0.0;
        for (i = 0;i<row;i++)
            for (j=0;j<col;j++){
                temp1+=T[i][j] * S[i][j];
                temp2+=T[i][j] * T[i][j];
        }
        
        omega = temp1/temp2;
        
        // xi = h + ωis
        for (i = 0;i<row;i++)
            for (j=0;j<col;j++)
                x_new[i][j] = H[i][j] + omega*S[i][j];
        
        //print_matrix(x_new,row,col,getName(x_new));
        //ri = s − ωit
         for (i = 0;i<row;i++)
            for (j=0;j<col;j++)
                r0[i][j] = S[i][j] - omega*T[i][j];
        
        //print_matrix(r0,row,col,getName(r0));
        L2_error = 0.0;
        L2_error = (norm_2_matrix(r0,row)/B_norm);
        if( (iteration_count > 1) && fabs(L2_error) < tolerance)
        {
            BiCGSTAB_end = clock();
            cout << "********* BICGSTAB  has converged **********" << endl;
            cout<<"Time Elapsed : "<< double (BiCGSTAB_end - BiCGSTAB_start)/ CLOCKS_PER_SEC<<"  ms "<< endl;
            cout << "Iteration :  " << iteration_count   <<endl<<" Error :  "<< L2_error << endl;
            cout<<"Tolerance Provided : " << tolerance <<endl;
            return x_new;
        }
        
        // Error Exception handling
        if( omega == 0.0 || rho == 0.0 )
            cout <<" BICGSTAB has not converged , it may be due to singularity in the given 'A' matrix " << endl;
        
        // to print the outout for each 100 Iterations 
        if(iteration_count % 100 == 0 )
            cout<<" BICGSTAB - at Iteration " << iteration_count << " Error is : "<<L2_error<< endl;
        
        // assign x to x-new and rho_new to rho
        rho = rho_new;
        x = x_new;
        iteration_count++;
        //cout<<"   ----- Iteration ended ---- " << endl;
    
    }
    BiCGSTAB_end = clock();
    cout<< " BICGSTAB Iteration has ****NOT**** converged to the Solution "<< endl;
    cout<<" Time Elapsed : "<< double (BiCGSTAB_end - BiCGSTAB_start)/ CLOCKS_PER_SEC<<"  ms "<< endl;
    cout << " Iteration :  " << iteration_count   <<endl<< " Error :  "<< L2_error << endl;
    
    return x_new;
    
    // delete the variables
    delete r0,r0_hat, p0,p0_new,v,v_new, H,x_new, S, T, Ax;
//     double rho,alpha,omega, rho_new,omega_new,error,beta;
    
    
}


void BiCGSTAB_Vector(double*& x_new,double **A, double *b,double *x,int row, int col , double tolerance , int max_iter)
{
    /* Algorithm taken from SIAM Paper :
    BI-CGSTAB: A FAST AND SMOOTHLY CONVERGING VARIANT OF BI-CG FOR THE SOLUTION OF NONSYMMETRIC LINEAR SYSTEMS* H. A. VAN DER VORSTt 
    2d matrrix details
    A = Input matrix - Square and Real
    b = Result matrix - suare and Real
    x =  intial guess matrix
    tolerance =  error tolerance for convergence check
    max_iter = maxximum Iteration;
    
    Algorithm used:
        r0 = b − Ax0
        Choose an arbitrary vector r̂0 such that (r̂0, r0) ≠ 0, e.g., r̂0 = r0 . Note that notation (x,y) applies for scalar product of vectors (x,y) = <x,y> = x·y = x ' y
        ρ0 = α = ω0 = 1
        v0 = p0 = 0
        For i = 1, 2, 3, …
        ρi = (r̂0, ri−1)
        β = (ρi/ρi−1)(α/ωi−1)
        pi = ri−1 + β(pi−1 − ωi−1vi−1)
        vi = Api
        α = ρi/(r̂0, vi)
        h = xi−1 + αpi
        If h is accurate enough, then set xi = h and quit
        s = ri−1 − αvi
        t = As
        ωi = (t, s)/(t, t)
        xi = h + ωis
        If xi is accurate enough, then quit
        ri = s − ωit
    */
    
    int i,j;
    int i_deviation,j_deviation;
    clock_t BiCGSTAB_start, BiCGSTAB_end;
    // Declaration 
    double *r0,*r0_hat;
    double *p0,*p0_new,*v,*v_new, *H, *S, *T, *Ax;
    double rho,alpha,omega, rho_new,omega_new,error,beta;
    std::cout<<" Done 1 " << std::endl;
    memset_1d(r0,col);memset_1d(r0_hat,col);
    memset_1d(p0,col);memset_1d(p0_new,col);memset_1d(v,col);memset_1d(v_new,col);
    memset_1d(H,col);memset_1d(S,col);memset_1d(T,col);memset_1d(Ax,col);
    
    double B_norm = norm_2_vector(b,row);
    
    
    // Stp -1 ; Calulate R0 , Assign it to R0_new
    // Calculate Ax
    matrix_vector_multiply(Ax,A,x,col);
    std::cout<<" Done 1 " << std::endl;
         
    for (i = 0; i< col ; i++){
            r0[i] = b[i] - Ax[i];
            r0_hat[i] = r0[i];
        }

    
    // decalatre values of Constants - Alpha , omega 
    alpha = rho = omega = rho= 1.0;
    
    // Set values of v0 and p0 as 0 Matrix, this would have been  done during declaration itself. 
    BiCGSTAB_start = clock();
    int iteration_count = 1;
    double L2_error = 1000;
    while ( (fabs(L2_error) > tolerance) && (iteration_count  <= max_iter))
    {
        //ρi = (r̂0, ri−1)
        rho_new = 0.0;
        rho_new = vector_inner_dot(r0_hat,r0,col);
        
         if( iteration_count  > 1 )
         {
            // β = (ρi/ρi−1)(α/ωi−1)
            beta = (rho_new/rho) * (alpha/omega);
            
            //Pi = ri−1 + β(Pi−1 − ωi−1Vi−1)
            for (i = 0;i<row;i++)
                    p0_new[i] = r0[i] + beta*(p0[i] -  (omega*v[i])) ;
         }
         else
             p0_new  = r0;
        
        // vi = APi
        matrix_vector_multiply(v_new,A,p0_new,col);
        
       
        //α = ρi/(r̂0, vi)clc
        double temp  = 0.0;
        temp =  vector_inner_dot(r0_hat,v_new,col);
        
        alpha = rho_new/temp;
        
        for (i = 0;i<col;i++)
                H[i] = x[i] + alpha*p0_new[i];
        
                // s = ri−1 − αvi
        for (i = 0;i<col;i++)
                S[i] = r0[i]-  alpha*v_new[i];
        
        // Early Exit Criteria , if the Norm of (S ) is less than tolerance , then exit the loop 
        L2_error = 0.0;
        L2_error = (norm_2_vector(S,col));
        if( ( fabs(L2_error) < tolerance))
        {
            BiCGSTAB_end = clock();
            cout << "********* BICGSTAB  has converged in "<< iteration_count << " in " << double (BiCGSTAB_end - BiCGSTAB_start)/ CLOCKS_PER_SEC << "time ********** " << endl;
            //print_matrix(H,row,col,getName(h));
            return ;
        }
        if( ( fabs(L2_error) > 1e20))
        {
            BiCGSTAB_end = clock();
            cout << "********* BICGSTAB  has NOT converged , Error blown up , Matrix might be Singular "<< endl;
            //print_matrix(H,row,col,getName(h));
            return ;
        }

        //print_matrix(S,row,col,getName(S));
        // t = As
         matrix_vector_multiply(T,A,S,col);
        //print_matrix(T,row,col,getName(T));

        //ωi = (t, s)/(t, t)
        double temp2 = 0.0;
        double temp1 = 0.0;
        temp1 = vector_inner_dot(T,S,col);
        temp2 = vector_inner_dot(T,T,col);
       
        omega = temp1/temp2;
        
        // xi = h + ωis
        for (i = 0;i<row;i++)
                x_new[i] = H[i] + omega*S[i];


        //print_matrix(x_new,row,col,getName(x_new));
        //ri = s − ωit
         for (i = 0;i<row;i++)
                r0[i] = S[i] - omega*T[i];
        
        //print_matrix(r0,row,col,getName(r0));
        L2_error = 0.0;
        L2_error = (norm_2_vector(r0,col)/B_norm);
        
        if( (iteration_count > 1) && fabs(L2_error) < tolerance)
        {
            BiCGSTAB_end = clock();
            cout << "********* BICGSTAB  has converged in "<< iteration_count << " in " << double (BiCGSTAB_end - BiCGSTAB_start)/ CLOCKS_PER_SEC << "time ********** " << endl;
            return ;
        }
        
        // Error Exception handling
        if( omega == 0.0 || rho == 0.0 )
            cout <<" BICGSTAB has not converged , it may be due to singularity in the given 'A' matrix " << endl;
        
        // to print the outout for each 100 Iterations 
        if(iteration_count % 10000 == 0 )
            cout<<" BICGSTAB - at Iteration " << iteration_count << " Error is : "<<L2_error<< endl;
        
        // assign x to x-new and rho_new to rho
        rho = rho_new;
        x = x_new;
        iteration_count++;
        //cout<<"   ----- Iteration ended ---- " << endl;
    }
    BiCGSTAB_end = clock();
    cout<< " BICGSTAB Iteration has ****NOT**** converged to the Solution "<< endl;
    cout<<" Time Elapsed : "<< double (BiCGSTAB_end - BiCGSTAB_start)/ CLOCKS_PER_SEC<<"  ms "<< endl;
    cout << " Iteration :  " << iteration_count   <<endl<< " Error :  "<< L2_error << endl;
    
    return ;
    
    // delete the variables
    delete r0,r0_hat, p0,p0_new,v,v_new, H, S, T, Ax;
}

double** FD_2D_Second_derivative_matrix(int Nx , int Ny)
{
    double **arr;
    double** diagonal;double** identity;double** zero;
    memset_2d(arr,Nx*Ny,Nx*Ny);
    memset_2d(diagonal,Nx,Ny);memset_2d(identity,Nx,Ny);memset_2d(zero,Nx,Ny);
    
    int block_x_id;int block_y_id;int block_x;int block_y;
       
    // fills the block matrices Individually ( A - Diagonal , B - Identity , C - Zero matrix )
    for (int i=0;i<Nx;i++)
        for(int j=0;j<Nx;j++){
            if(i==j){
                diagonal[i][j] = -4.0;
                identity[i][j] = 1.0;
                zero[i][j] = 0.0;
            }
            else if(abs(i-j) == 1){
                diagonal[i][j] = 1.0;
                identity[i][j] = 0.0;
                zero[i][j] = 0.0;
            }
            else{
                diagonal[i][j] = 0.0;
                identity[i][j] = 0.0;
                zero[i][j] = 0.0;
            }
        }
        
    // Now fils the Blocks of the overall  matrix  based on the block id's of the matrix 
    for ( int block_x_id = 0; block_x_id < Nx ; block_x_id++)
    {
        for(int block_y_id = 0 ; block_y_id<Ny; block_y_id++)
        {
            if(block_x_id == block_y_id)
            {
                // Diagonal Block - diagonal matrix
                for(int i = 0 ; i < Nx ; i++)
                    for(int j=0;j<Ny;j++)
                        arr[(block_y_id*Ny) + i][(block_x_id*Nx) + j] = diagonal[i][j];
                    
            }
            else if(abs(block_x_id - block_y_id) == 1)
            {
                 // Sub Diagonal - Identity
                for(int i = 0 ; i < Nx ; i++)
                    for(int j=0;j<Ny;j++){ 
                        
                        arr[(block_y_id*Ny) + i][(block_x_id*Nx) + j] = identity[i][j];
                    }
            }
            
            else
            {
                for(int i = 0 ; i < Nx ; i++)
                    for(int j=0;j<Ny;j++)
                        arr[(block_y_id*Ny) + i][(block_x_id*Nx) + j] = zero[i][j];
            }
        }
    }
    
    return arr;
    delete arr,diagonal,identity,zero;
}


void BiCGSTAB_Vector_CSR(std::vector<double>& x_new,Matrix_1D *A, std::vector<double>& b,std::vector<double>& x,int row, int col , double tolerance , int max_iter)
{
    /* Algorithm taken from SIAM Paper :
    BI-CGSTAB: A FAST AND SMOOTHLY CONVERGING VARIANT OF BI-CG FOR THE SOLUTION OF NONSYMMETRIC LINEAR SYSTEMS* H. A. VAN DER VORSTt 
    2d matrrix details
    A = Input matrix - Square and Real
    b = Result matrix - suare and Real
    x =  intial guess matrix
    tolerance =  error tolerance for convergence check
    max_iter = maxximum Iteration; */

    
    int i,j;
    int i_deviation,j_deviation;
    clock_t BiCGSTAB_start, BiCGSTAB_end;
    // Declaration 
    double *r0,*r0_hat;
    double *p0,*p0_new,*v,*v_new, *H, *S, *T, *Ax;
    double rho,alpha,omega, rho_new,omega_new,error,beta;
    
    memset_1d(r0,col);memset_1d(r0_hat,col);
    memset_1d(p0,col);memset_1d(p0_new,col);memset_1d(v,col);memset_1d(v_new,col);
    memset_1d(H,col);memset_1d(S,col);memset_1d(T,col);memset_1d(Ax,col);
    
    
    double B_norm = norm_2_vector(&b[0],row);
    std::cout<<" DONE 1   col : " << col<<std::endl;   
    
    // Stp -1 ; Calulate R0 , Assign it to R0_new
    // Calculate Ax
    sp_matrix_vector_multiply(Ax,A,&x[0],col);
      
    for (i = 0; i< col ; i++){
            r0[i] = b[i] - Ax[i];
            r0_hat[i] = r0[i];
        }

    
    // decalatre values of Constants - Alpha , omega 
    alpha = rho = omega = rho= 1.0;
    
    // Set values of v0 and p0 as 0 Matrix, this would have been  done during declaration itself. 
    BiCGSTAB_start = clock();
    int iteration_count = 1;
    double L2_error = 1000;
    while ( (fabs(L2_error) > tolerance) && (iteration_count  <= max_iter))
    {
        //ρi = (r̂0, ri−1)
        rho_new = 0.0;
        rho_new = vector_inner_dot(r0_hat,r0,col);
        
         if( iteration_count  > 1 )
         {
            // β = (ρi/ρi−1)(α/ωi−1)
            beta = (rho_new/rho) * (alpha/omega);
            
            //Pi = ri−1 + β(Pi−1 − ωi−1Vi−1)
            for (i = 0;i<row;i++)
                    p0_new[i] = r0[i] + beta*(p0[i] -  (omega*v[i])) ;
         }
         else
             p0_new  = r0;
        
        // vi = APi
        sp_matrix_vector_multiply(v_new,A,p0_new,col);
        
       
        //α = ρi/(r̂0, vi)clc
        double temp  = 0.0;
        temp =  vector_inner_dot(r0_hat,v_new,col);
        
        alpha = rho_new/temp;
        
        for (i = 0;i<col;i++)
                H[i] = x[i] + alpha*p0_new[i];
        
                // s = ri−1 − αvi
        for (i = 0;i<col;i++)
                S[i] = r0[i]-  alpha*v_new[i];
        
        // Early Exit Criteria , if the Norm of (S ) is less than tolerance , then exit the loop 
        L2_error = 0.0;
        L2_error = (norm_2_vector(S,col));
        if( ( fabs(L2_error) < tolerance))
        {
            BiCGSTAB_end = clock();
            cout << "********* BICGSTAB  has converged in "<< iteration_count << " in " << double (BiCGSTAB_end - BiCGSTAB_start)/ CLOCKS_PER_SEC << "time ********** " << endl;
            //print_matrix(H,row,col,getName(h));
            return ;
        }
        if( ( fabs(L2_error) > 1e20))
        {
            BiCGSTAB_end = clock();
            cout << "********* BICGSTAB  has NOT converged , Error blown up , Matrix might be Singular "<< endl;
            //print_matrix(H,row,col,getName(h));
            return ;
        }

        //print_matrix(S,row,col,getName(S));
        // t = As
        sp_matrix_vector_multiply(T,A,S,col);
        //print_matrix(T,row,col,getName(T));

        //ωi = (t, s)/(t, t)
        double temp2 = 0.0;
        double temp1 = 0.0;
        temp1 = vector_inner_dot(T,S,col);
        temp2 = vector_inner_dot(T,T,col);
       
        omega = temp1/temp2;
        
        // xi = h + ωis
        for (i = 0;i<row;i++)
                x_new[i] = H[i] + omega*S[i];


        //print_matrix(x_new,row,col,getName(x_new));
        //ri = s − ωit
         for (i = 0;i<row;i++)
                r0[i] = S[i] - omega*T[i];
        
        //print_matrix(r0,row,col,getName(r0));
        L2_error = 0.0;
        L2_error = (norm_2_vector(r0,col)/B_norm);
        
        if( (iteration_count > 1) && fabs(L2_error) < tolerance)
        {
            BiCGSTAB_end = clock();
            cout << "********* BICGSTAB  has converged in "<< iteration_count << " in " << double (BiCGSTAB_end - BiCGSTAB_start)/ CLOCKS_PER_SEC << "time ********** " << endl;
            return ;
        }
        
        // Error Exception handling
        if( omega == 0.0 || rho == 0.0 )
            cout <<" BICGSTAB has not converged , it may be due to singularity in the given 'A' matrix " << endl;
        
        // to print the outout for each 100 Iterations 
        if(iteration_count % 10000 == 0 )
            cout<<" BICGSTAB - at Iteration " << iteration_count << " Error is : "<<L2_error<< endl;
        
        // assign x to x-new and rho_new to rho
        rho = rho_new;
        x = x_new;
        iteration_count++;
        //cout<<"   ----- Iteration ended ---- " << endl;
    }
    BiCGSTAB_end = clock();
    cout<< " BICGSTAB Iteration has ****NOT**** converged to the Solution "<< endl;
    cout<<" Time Elapsed : "<< double (BiCGSTAB_end - BiCGSTAB_start)/ CLOCKS_PER_SEC<<"  ms "<< endl;
    cout << " Iteration :  " << iteration_count   <<endl<< " Error :  "<< L2_error << endl;
    
    return ;
    
    // delete the variables
    delete r0,r0_hat, p0,p0_new,v,v_new, H, S, T, Ax;
}