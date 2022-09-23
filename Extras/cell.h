#ifndef __MATRIX_1D__
#define __MATRIX_1D__
class Matrix_1D
{
private:
    /* data */
public:
    int sizeRowPtr;
    std::vector<int> rowPtr;
    std::vector<int> colPtr;
    std::vector<double> values;
    

    double& getValues(int i , int j);

   // double& operator()(int row, int col);   

    void setSparsityPattern(int N_Nodes_Element, int N_Cells, std::vector<int> local2Global);
};

double& Matrix_1D::getValues(int i , int j)
{
    double zero = 0;
    int start = rowPtr[i];
    int end = rowPtr[i+1];
    for (int index = start ; index < end ; index++)
        if(colPtr[index]== j)
            return values[index];
    return zero;
}   

void Matrix_1D::setSparsityPattern(int N_Nodes_Element, int N_Cells, std::vector<int> local2Global )
{
    sizeRowPtr = (N_Nodes_Element-1)*N_Cells ;
    int N_DOF = sizeRowPtr -1;
    int order = N_Nodes_Element - 1;

    // all edge nodes =  5, INternal nodes = NNE 
    // Edge Nodes =  (N_cells - 1) -2 ( For start and End Correction)    - No of Shared Nodes = 2 * NNE - 1
    // INternal Nodes =  N_cells * (NNE-1) + 2 ( For start and the end )  - No of Shared Nodes =   NNE
    int EdgeNodes =  N_Cells - 1 - 2;
    int IntNodes = (N_Cells * order) + 2 ;
    int NNZSize = EdgeNodes*((2*N_Nodes_Element) - 1 );
    NNZSize += IntNodes * N_Nodes_Element;

    // resize the ColPtr and the values array 
    colPtr.resize(NNZSize,0);
    values.resize(NNZSize,0);
    rowPtr.resize(sizeRowPtr,0);
    rowPtr[0] = 0;
    
    std::cout<< " N_DOF " << N_DOF;
    std::cout<<"NNZ Values : " << NNZSize << std::endl;
    
    // Loop through all the DOF's
    int index = 0;
    int cell = 0;
    int start =0;
    int nodesSharedbyEdgenode =  (2*N_Nodes_Element) - 1 ;  // NUmber of DOF shared by edge nodes
    for ( int node = 0 ; node < N_DOF; node ++)
    {
       cell = ceil(node / order) - 1;
       index = N_Nodes_Element*cell;
       start = rowPtr[node-1];
        if(node == 0 || node == N_DOF -1 || node % order !=0  ) // Internal Nodes
        {
            rowPtr[node] = start + N_Nodes_Element;
            for(int i=0;i<N_Nodes_Element;i++)
                colPtr[start + i] = local2Global[index + i];
        }

        else // Edge Nodes
        {
            rowPtr[node] = rowPtr[node-1] + nodesSharedbyEdgenode;
            for(int i=0;i<=nodesSharedbyEdgenode;i++){
                colPtr[start + i] = local2Global[index + i];
                if(i== N_Nodes_Element)
                    i++;
        }

    }
}

#endif
