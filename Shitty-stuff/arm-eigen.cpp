#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "armadillo"

#define DIM (1000)

using namespace arma;
using namespace std;

int main(int argc, char *argv[]) {
    
    // for matrices with real elements
    
    mat A = randu<mat>(DIM,DIM);
    mat B = A.t()*A;  // generate a symmetric matrix
    
    vec eigval;
    mat eigvec;
    
    eig_sym(eigval, eigvec, B);


    
//    cout << "Eigenvalues\n";
//    eigval.print();
//    vec eigv1 = eigvec.col(0);
    
//    cout << "Eigenvector\n";
//    eigv1.print();
//    //takes the eigenvalues vector, the eigenvectors matrix and the "to-be-diag" matrix
    
    return 0;
}
