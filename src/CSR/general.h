void csrmult1(CSR*csr,double*x,double*b){
	//Assumes b already contains valuable information.
	int i;
	for(i=0;i<csr->n;i++){
		b[csr->ri[i]]+=csr->v[i]*x[csr->ci[i]];
	}
}
void csrmult1symmetric(CSR*csr,double*x,double*b){
	//Assumes b already contains valuable information.
	int i;
	int row,col;
	for(i=0;i<csr->n;i++){
		row=csr->ri[i];col=csr->ci[i];
		b[row]+=csr->v[i]*x[col];
		b[col]+=csr->v[i]*x[row];
	}
}
