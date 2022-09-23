#include "Sparse_Matrix.h"
#include <fstream>
#include<algorithm>


double& Sparse_Matrix::getValues(int i , int j)
{
    int start = rowPtr[i];
    int end = rowPtr[i+1];
    for (int index = start ; index < end ; index++)
        if(colPtr[index]== j)
            return values[index];
    return zero;
}   

double& Sparse_Matrix::operator()(int i, int j)
{
    std::cout<<" ----- Entering operator"<<std::endl;
    getValues(i,j);
}


void Sparse_Matrix::setSparsityPattern(int N_Nodes_Element, int N_Cells, std::vector<int> local2Global )
{
    sizeRowPtr = (N_Nodes_Element-1)*N_Cells + 2 ;
    int N_DOF = sizeRowPtr - 1;
    int order = N_Nodes_Element - 1;

    // all edge nodes =  5, INternal nodes = NNE 
    // Edge Nodes =  (N_cells - 1) -2 ( For start and End Correction)    - No of Shared Nodes = 2 * NNE - 1
    // INternal Nodes =  N_cells * (NNE-1) + 2 ( For start and the end )  - No of Shared Nodes =   NNE
    int EdgeNodes =  N_Cells + 1 - 2;
    int IntNodes = (N_Cells * (order-1)) + 2 ;
    NNZSize = EdgeNodes*((2*N_Nodes_Element) - 1 );
    NNZSize += IntNodes * N_Nodes_Element;

    // resize the ColPtr and the values array 
    colPtr.resize(NNZSize,0);
    values.resize(NNZSize,0);
    rowPtr.resize(sizeRowPtr,0);
   // rowPtr[0] = 0;
    
    
    // Loop through all the DOF's
    int index = 0;
    int cell = 0;
    int start =0;
    int nodesSharedbyEdgenode =  (2*N_Nodes_Element) - 1 ;  // NUmber of DOF shared by edge nodes
    for ( int node = 0 ; node < N_DOF; node ++)
    {
       cell = int(node / order);
       if(node%order == 0 && node !=0) cell -= 1;
       index = N_Nodes_Element*cell;
       start = rowPtr[node];
        if(node == 0 || node == N_DOF -1 || node % order !=0  ) // Internal Nodes
        {
            rowPtr[node+1] = start + N_Nodes_Element;

            for(int i=0;i<N_Nodes_Element;i++){
                colPtr[start + i] = local2Global[index + i];
            }
        }

        else // Edge Nodes
        {
            rowPtr[node+1] = start + nodesSharedbyEdgenode;
            int j=0;
            for(int i=0;i<=nodesSharedbyEdgenode;i++){
                colPtr[start + j] = local2Global[index + i];
                j++;
                if(i==order)
                    i++;
            }

            
        }
    }   
    setMatrixDimensions(N_DOF,N_DOF,NNZSize);


    // for ( int i=0;i<sizeRowPtr;i++)
    //     std::cout<<rowPtr[i]<<"\t";
    
    // std::cout<<std::endl;

    // for ( int i=0;i<NNZSize;i++)
    //     std::cout<<colPtr[i]<<"\t";
    
    // std::cout<<std::endl; 

    // Output :
    std::cout<<" Sparsity Pattern Set for Matrix " <<std::endl;

    // SORT THE COLUMN INDEX of the system 
    for (int i = 0 ; i < N_DOF;i++){
        std::sort(colPtr.data() + rowPtr[i],colPtr.data() + rowPtr[i+1] );
    }
    
    //     for ( int i=0;i<sizeRowPtr;i++)
    //     std::cout<<rowPtr[i]<<"\t";
    
    // std::cout<<std::endl;

    // for ( int i=0;i<NNZSize;i++)
    //     std::cout<<colPtr[i]<<"\t";
    
    // std::cout<<std::endl; 

}

void Sparse_Matrix::deleteMatrix()
{
    rowPtr.clear();
    colPtr.clear();
    values.clear();
}

void Sparse_Matrix::setMatrixDimensions(int _row, int _col)
{
    row =_row;
    col =_col;
}

void Sparse_Matrix::setMatrixDimensions(int _row , int _col,int _nnz)
{
    row =_row;
    col =_col;
    NNZSize = _nnz;

}

void Sparse_Matrix::Display_matrix()
{
    // Rowptr
    std::cout << std::endl;

    std::cout << " ROW : " << row<<std::endl;

    std::cout << " COL : " << col<<std::endl;

    for ( int i = 0 ; i < row+1 ; i ++)
        std::cout << rowPtr[i] <<"     ";
    std::cout << std::endl;

    for ( int i = 0 ; i < NNZSize ; i ++)
        std::cout << colPtr[i] <<"     ";
    std::cout << std::endl;

    for ( int i = 0 ; i < NNZSize ; i ++)
        std::cout << values[i] <<"     ";
    std::cout << std::endl;
    std::cout << std::endl;

    // Display as Matrix 
    for ( int i = 0 ; i < row ; i ++)
    {
        for ( int j = 0 ; j < col ; j++)
        {
            std::cout << getValues(i,j)<< "\t";
        }
        std::cout<< std::endl;
    }

}

void Sparse_Matrix::Display_matrix_arr()
{
    // Rowptr
    std::cout << std::endl;

    std::cout << " ROW : " << row<<std::endl;

    std::cout << " COL : " << col<<std::endl;

    for ( int i = 0 ; i < row+1 ; i ++)
        std::cout << rowPtr_c[i] <<"     ";
    std::cout << std::endl;

    for ( int i = 0 ; i < NNZSize ; i ++)
        std::cout << colPtr_c[i] <<"     ";
    std::cout << std::endl;

    for ( int i = 0 ; i < NNZSize ; i ++)
        std::cout << values_c[i] <<"     ";
    std::cout << std::endl;
    std::cout << std::endl;

    // Display as Matrix 
    for ( int i = 0 ; i < row ; i ++)
    {
        for ( int j = 0 ; j < col ; j++)
        {
            std::cout << getValues(i,j)<< "\t";
        }
        std::cout<< std::endl;
    }

}



// Insert Values into Row Pointer and Column Array
// CAUTION , THIS HAS TO BE FOLLOWED BY AN ROW REDUCTION OPERATION
void Sparse_Matrix::insertValues(int _i , int _j,double _val, int _tmp_nnz )
{
   
    rowPtr[_i+1] += 1;
    colPtr[_tmp_nnz] = _j;
    values[_tmp_nnz] = _val;
    
}

void Sparse_Matrix::UpdateRowPointer(int _row, int _col,int _nnz)
{  

    NNZSize = _nnz;
    values.resize(_nnz);
    colPtr.resize(_nnz);
    rowPtr.resize(row+1);
    rowPtr[0] = 0;
    for ( int i = 1 ; i <= row ; i++){
        rowPtr[i] += rowPtr[i-1]; 
        //if(rowPtr[i] == rowPtr[i-1] ) row = row -1;
    }
   
}

void Sparse_Matrix::importValuesfromFile(std::string filename)
{
    std::ifstream file;
    file.open(filename);
    int _row;
    int _col;
    int _nnz;
    file>>_row;  file>>_col;  file>>_nnz;

    int k;

    row = _row;
    col = _col;
    NNZSize = _nnz;
    colPtr.resize(NNZSize);
    values.resize(NNZSize);
    rowPtr.resize(row+1);

    for(int i = 0;i < _nnz;i++)
    {
        file>>k;
        rowPtr[k+1]++;
        file>>colPtr[i];
        file>>values[i];
    }

    for(int i = 0;i < _row;i++)
    {
       rowPtr[i+1] += rowPtr[i];
    }

 

    file.close();

}