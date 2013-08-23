void csr(CSR*csr,double*x,double*b,int length){
	//Adds on to b.
	int i;
	initializeVector(b,length,0);
	for(i=0;i<csr->n;i++){
		b[csr->ri[i]]+=csr->v[i]*x[csr->ci[i]];
	}
}
void csrTranspose(CSR*csr,double*x,double*b,int length){
	//Adds on to b.
	int i;
	initializeVector(b,length,0);
	for(i=0;i<csr->n;i++){
		b[csr->ci[i]]+=csr->v[i]*x[csr->ri[i]];
	}
}
void csrSymmetric(CSR*csr,double*x,double*b){
	//Adds on to b.
	int i;
	int row,col;
	for(i=0;i<csr->n;i++){
		row=csr->ri[i];col=csr->ci[i];
		b[row]+=csr->v[i]*x[col];
		b[col]+=csr->v[i]*x[row];
	}
}
