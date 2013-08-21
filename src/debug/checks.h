/*
 * checks.h
 *
 *  Created on: Jul 16, 2013
 *      Author: lordvon
 */

#ifndef CHECKS_H_
#define CHECKS_H_

void compareCsrSlA(SymLap*slcsr,SymLap*slA){
	int i,r,c,dim;
	double tol=1e-3,diff;
	dim=slcsr->dim;
	if(slcsr->dim!=slA->dim){
		printf("sl dims do not match!!!\n");
		return;
	}
	for(i=0;i<slcsr->csr.nn;i++){
		r=slcsr->csr.ri[i];
		c=slcsr->csr.ci[i];
		diff=slcsr->csr.v[i]-slA->A[r*dim+c];
		if(fabs(diff)>tol){
			printf("CSR does not match sl.A!!! (%.2e vs. %.2e, element %d)\n",slcsr->csr.v[i],slA->A[r*dim+c],i);
		}
	}
	printf("sl.A / CSR comparison completed.\n");
}


void printvec(double * vec, int len, char * name){
	int j;
	printf("Printing \"%s\":\n", name);
	for(j=0;j<len;j++){
		printf("Row %d:\t",j);
		printf("%e\n",vec[j]);
	}
}
void printInterfaces(Interfaces* interfaces){
	int i;
	for(i=0;i<interfaces->totalinterfaces;i++){
		printf("Interface %d: %d %d %d %d %d\n",i,interfaces->block[i],interfaces->side[i],interfaces->adjacentBlock[i],interfaces->xsize[i],interfaces->rowstart[i]);
	}
}
void printDim(Grid* grid){
	int b;
	for(b=0;b<grid->totalblocks;b++){
		printf("Block %d:\n",b);
		printf("xdim: %d %d %d %d\n",grid->xdim[b][0],grid->xdim[b][1],grid->xdim[b][2],grid->xdim[b][3]);
		printf("vdim: %d %d %d %d\n",grid->vdim[b][0],grid->vdim[b][1],grid->vdim[b][2],grid->vdim[b][3]);
		printf("hdim: %d %d %d %d\n",grid->hdim[b][0],grid->hdim[b][1],grid->hdim[b][2],grid->hdim[b][3]);
		printf("cdim: %d %d %d %d\n",grid->cdim[b][0],grid->cdim[b][1],grid->cdim[b][2],grid->cdim[b][3]);
	}
}
void printBlockNames(Grid* grid){
	int b;
	printf("Block File Names:\n");
	for(b=0;b<grid->totalblocks;b++){
		printf("%s\n",grid->blockNames[b]);
	}
}
void printbc(Grid* grid,BoundaryConditions* bc){
	int b;
	printf("Block BCs:\n");
	for(b=0;b<grid->totalblocks;b++){
		printf("Block %d: %d %d %d %d\n",b,bc->id[b][0],bc->id[b][1],bc->id[b][2],bc->id[b][3]);
	}
}
void printinterfacepoints(Interfaces* interfaces){
	int i,ip;
	for(i=0;i<interfaces->totalinterfaces;i++){
		printf("Interface %d:\n",i);
		for(ip=0;ip<interfaces->xsize[i];ip++){
			printf("Point %d: %d %d\n",ip,interfaces->blockx[i][ip],interfaces->adjacentBlockx[i][ip]);
		}
	}
}
void printinterfaceedges(Interfaces* interfaces){
	int i,ie;
	for(i=0;i<interfaces->totalinterfaces;i++){
		printf("Interface %d:\n",i);
		for(ie=0;ie<interfaces->xsize[i]-1;ie++){
			printf("Edge %d: %d %d\n",ie,interfaces->e[i][ie],interfaces->ae[i][ie]);
		}
	}
}
int comparemaindiag(double *CTC,LinearSystem* ls){
	int i;
	int match=1;
	for(i=0;i<ls->ccols;i++){
		if(CTC[i+i*ls->ccols]!=ls->ctcmain[i]){
			match=0;
			printf("Main diag did not match at %d\n",i);
			break;
		}
	}
	return match;
}
int comparectcmult(double *CTC,Grid* grid,Interfaces* interfaces,LinearSystem* ls){
	int xi,match=1;
	double *s=malloc(grid->totalnodes*sizeof(double));
	for(xi=0;xi<grid->totalnodes;xi++){
		s[xi]=xi;
	}
	double *r1=matmult(CTC,s,grid->totalnodes,grid->totalnodes,1);
	double *r2=malloc(grid->totalnodes*sizeof(double));
	initv(r2,grid->totalnodes,0);
	ctcmult(s,r2,grid,interfaces,ls);
	for(xi=0;xi<grid->totalnodes;xi++){
		if(r1[xi]!=r2[xi]){
			printf("Mismatch at %d (%g vs %g)\n",xi,r1[xi],r2[xi]);
			match=0;
			break;
		}
	}
	free(s);free(r1);free(r2);
	return match;
}
int comparectmult(double *CT,Grid* grid,Interfaces* interfaces,LinearSystem* ls){
	int xi,match=1;
	double *s=malloc(ls->crows*sizeof(double));
	for(xi=0;xi<ls->crows;xi++){
		s[xi]=xi;
	}
	double *r1=matmult(CT,s,ls->ccols,ls->crows,1);
	double *r2=malloc(ls->ccols*sizeof(double));
	initv(r2,ls->ccols,0);
	ctmult(s,r2,grid,interfaces,ls);
	for(xi=0;xi<ls->ccols;xi++){
		if(r1[xi]!=r2[xi]){
			printf("Mismatch at %d (%g vs %g)\n",xi,r1[xi],r2[xi]);
			match=0;
			break;
		}
	}
	free(s);free(r1);free(r2);
	return match;
}
int comparecmult(double *C,Grid* g,Interfaces* is,LinearSystem* ls){
	int xi,match=1;
	double *s=malloc(ls->ccols*sizeof(double));
	for(xi=0;xi<ls->ccols;xi++){
		s[xi]=xi;
	}
	double *r1=matmult(C,s,ls->crows,ls->ccols,1);
	double *r2=malloc(ls->crows*sizeof(double));
	initv(r2,ls->crows,0);
	cmult(s,r2,g,is,ls);
	for(xi=0;xi<ls->crows;xi++){
		if(r1[xi]!=r2[xi]){
			printf("Mismatch at %d (%g vs %g)\n",xi,r1[xi],r2[xi]);
			match=0;
			break;
		}
	}
	free(s);free(r1);free(r2);
	return match;
}
int compareslmult(Grid* g,Interfaces* is,SymLap*sl){
	int ci,match=1;
	int d=sl->dim;
	for(ci=0;ci<d;ci++){ sl->x[ci]=ci; }
	double *r1=matmult(sl->A,sl->x,d,d,1);
	initv(sl->b,d,0);
	slmultCSR(sl->x,sl->b,sl,g,is);
	for(ci=0;ci<d;ci++){
		if(r1[ci]!=sl->b[ci]){
			printf("Mismatch at %d (%.5e vs %.5e)\n",ci,r1[ci],sl->b[ci]);
			//printf("Mismatch at %d (%g vs %g)\n",ci,r1[ci],sl->b[ci]);
			match=0;
		}
	}
	free(r1);
	return match;
}
void gridinfo(Grid*g){

}

#endif /* CHECKS_H_ */
