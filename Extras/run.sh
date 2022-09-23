rm Prob_1.o
rm Matrix_1D.o
icc -g -c Prob_1.cpp
icc -g -c Matrix_1D.cpp
icc -g -c FE_1D.cpp
icc -g -c Iterative_Solvers.cpp
icc -g -o prob1 Prob_1.o Matrix_1D.o Iterative_Solvers.o -lmkl_rt
