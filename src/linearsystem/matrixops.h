/*
 * matrixops.h
 *
 *  Created on: Jul 10, 2013
 *      Author: lordvon
 */

#ifndef MATRIXOPS_H_
#define MATRIXOPS_H_

double * identity(int dim){
	double * I=(double *)malloc(sizeof(double)*dim*dim);
	int i,j,jj;
	for(j=0;j<dim;j++){
		jj=j*dim;
		for(i=0;i<dim;i++){
			I[i+jj]=0;
			if(i==j){
				I[i+jj]=1;
			}
		}
	}
	return I;
}
int checkSymmetry(double * mat, int dim){
	int i,j,symmetric=1;
	for(j=0;j<dim;j++){
		for(i=0;i<dim;i++){
			if(mat[i+j*dim]!=mat[j+i*dim]){
				symmetric=0;
				printf("Asymmetry at %d, %d\n",i,j);
			}
		}
	}
	return symmetric;
}
double * makeAT(double * A, int rows, int cols){
	int i,j;
	int a;
	double * AT=(double *)malloc(sizeof(double)*cols*rows);
	for(j=0;j<rows;j++){
		a=j*cols;
		for(i=0;i<cols;i++){
			AT[j+rows*i]=A[i+a];
		}
	}
	return AT;
}
double * matmult(double * mat1, double * mat2, int outer1, int inner, int outer2){
	//multiplies two matrices.
	double * result=(double *)malloc(sizeof(double)*outer1*outer2);
	int i,j,k;
	double total;
	for(j=0;j<outer1;j++){
		for(i=0;i<outer2;i++){
			total=0;
			for(k=0;k<inner;k++){
				total+=mat1[j*inner+k]*mat2[i+k*outer2];
			}
			result[i+j*outer2]=total;
		}
	}
	return result;
}
double * makeATA(double * A, int rows, int cols){
	double * AT=makeAT(A,rows,cols);
	double * ATA=matmult(AT,A,cols,rows,cols);
	return ATA;
}
double * makeAAT(double * A, int rows, int cols){
	double * AT=makeAT(A,rows,cols);
	double * AAT=matmult(A,AT,rows,cols,rows);
	return AAT;
}

#endif /* MATRIXOPS_H_ */
