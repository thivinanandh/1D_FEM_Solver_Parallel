#include<vector>
#include<cmath>
#include<iostream>


#ifndef __Sparse_Matrix__
#define __Sparse_Matrix__
class Sparse_Matrix
{
private:
    // int sizeRowPtr;

    double zero = 0;
    int temp_NNZ = 0;
    int temp_row = 0;
    int temp_col = 0;

public:
    int sizeRowPtr ;
    std::vector<int> rowPtr;
    std::vector<int> colPtr;
    std::vector<double> values;
    long int NNZSize;
    int row;
    int col;
    int* rowPtr_c;
    int* colPtr_c;
    double* values_c;



    // CONSTRUCTOR
    Sparse_Matrix()
    {
        
    }

    Sparse_Matrix(const Sparse_Matrix &p2 )
    {
        std::cout << " Called Copy Sonctryctor " <<std::endl;
        sizeRowPtr = p2.sizeRowPtr;
        rowPtr = p2.rowPtr;
        colPtr = p2.colPtr;
        values = p2.values;
        NNZSize = p2.NNZSize;
        row = p2.row;
        col = p2.col;
    }

    Sparse_Matrix(int _row, int _col, int _NNZ )
    {
        row = _row;
        col = _col;
        NNZSize = _NNZ;
    }

    double& getValues(int i , int j);

    double& operator()(int row, int col);   

    void setSparsityPattern(int N_Nodes_Element, int N_Cells, std::vector<int> local2Global);

    void deleteMatrix();

    void setMatrixDimensions(int _row , int _col);

    void setMatrixDimensions(int _row , int _col,int _nnz);

    void insertValues(int i , int j,double val,int tmp_nnz );

    void UpdateRowPointer(int row,int col,int nnz);
    
    void importValuesfromFile(std::string filename);

    //Check Matrix
    void Display_matrix();

    void Display_matrix_arr();

};


#endif
