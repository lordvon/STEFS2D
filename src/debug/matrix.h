/*
 * matrix.h
 *
 *  Created on: Jul 23, 2013
 *      Author: lordvon
 */

#ifndef MATRIX_H_
#define MATRIX_H_


void printmi(double * m, int rows, int cols, char * name){
	int i,j;
	printf("Printing \"%s\":\n", name);
	for(j=0;j<cols;j++){
		printf("\tCol%d",j);
	}
	printf("\n");
	for(j=0;j<rows;j++){
		printf("Row%d:\t",j);
		for(i=0;i<cols;i++){
			printf("%d\t",(int)m[j*cols+i]);
		}
		printf("\n");
	}
}
void printJxx(Grid* g){
	int b,i,j,xi;
	int xd0,xd1;
	printf("PRINTING Jxx\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		xd0=g->xdim[b][0];
		xd1=g->xdim[b][1];
		for(j=0;j<xd1;j++){
			for(i=0;i<xd0;i++){
				xi=i+j*xd0;
				printf("\t%.3g",g->Jxx[b][xi]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printhareas(Grid* g){
	int b,i,j,hi;
	int hd0,hd1;
	printf("PRINTING hareas\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		hd0=g->hdim[b][0];
		hd1=g->hdim[b][1];
		for(j=0;j<hd1;j++){
			for(i=0;i<hd0;i++){
				hi=i+j*hd0;
				printf("\t%.3g",g->hareas[b][hi]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printvareas(Grid* g){
	int b,i,j,vi;
	int vd0,vd1;
	printf("PRINTING vareas\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				printf("\t%.3g",g->vareas[b][vi]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printu(State* s,Grid* g){
	int b,i,j,vi;
	int vd0,vd1;
	printf("PRINTING u\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				printf("\t%.3g",s->u[b][vi]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printv(State* s,Grid* g){
	int b,i,j,hi;
	int hd0,hd1;
	printf("PRINTING v\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		hd0=g->hdim[b][0];
		hd1=g->hdim[b][1];
		for(j=0;j<hd1;j++){
			for(i=0;i<hd0;i++){
				hi=i+j*hd0;
				printf("\t%.3g",s->v[b][hi]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printuxx(Momentum* mm,Grid* g){
	int b,i,j,xi;
	int xd0,xd1;
	printf("PRINTING uxx\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		xd0=g->xdim[b][0];
		xd1=g->xdim[b][1];
		for(j=0;j<xd1;j++){
			for(i=0;i<xd0;i++){
				xi=i+j*xd0;
				printf("\t%.3g",mm->uxx[b][xi]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printvxx(Momentum* mm,Grid* g){
	int b,i,j,xi;
	int xd0,xd1;
	printf("PRINTING vxx\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		xd0=g->xdim[b][0];
		xd1=g->xdim[b][1];
		for(j=0;j<xd1;j++){
			for(i=0;i<xd0;i++){
				xi=i+j*xd0;
				printf("\t%.3g",mm->vxx[b][xi]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printucartxx(Momentum* mm,Grid* g){
	int b,i,j,xi;
	int xd0,xd1;
	printf("PRINTING ucartxx\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		xd0=g->xdim[b][0];
		xd1=g->xdim[b][1];
		for(j=0;j<xd1;j++){
			for(i=0;i<xd0;i++){
				xi=i+j*xd0;
				printf("\t%.3g",mm->ucartxx[b][xi]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printvcartxx(Momentum* mm,Grid* g){
	int b,i,j,xi;
	int xd0,xd1;
	printf("PRINTING vcartxx\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		xd0=g->xdim[b][0];
		xd1=g->xdim[b][1];
		for(j=0;j<xd1;j++){
			for(i=0;i<xd0;i++){
				xi=i+j*xd0;
				printf("\t%.3g",mm->vcartxx[b][xi]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printuxxx(Momentum* mm,Grid* g){
	int b,i,j,xi;
	int xd0,xd1;
	printf("PRINTING uxxx\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		xd0=g->xdim[b][0];
		xd1=g->xdim[b][1];
		for(j=0;j<xd1;j++){
			for(i=0;i<xd0;i++){
				xi=i+j*xd0;
				printf("\t%.3g",mm->uxxx[b][xi]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printcontrac(Grid* g){
	int b,i,j,xi;
	int xd0,xd1;
	printf("PRINTING contrac\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		xd0=g->xdim[b][0];
		xd1=g->xdim[b][1];
		for(j=0;j<xd1;j++){
			for(i=0;i<xd0;i++){
				xi=i+j*xd0;
				printf("\t%.3g",g->contrac11[b][xi]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printcartdiffvv(Momentum* mm,Grid* g){
	int b,i,j,vi;
	int vd0,vd1;
	printf("PRINTING cartvv\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				printf("\t%.3g",mm->vcartvv[b][vi]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printvvv(Momentum* mm,Grid* g){
	int b,i,j,vi;
	int vd0,vd1;
	printf("PRINTING vvv\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				printf("\t%.3g",mm->vvv[b][vi]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printcarthh(Momentum* mm,Grid* g){
	int b,i,j,hi;
	int hd0,hd1;
	printf("PRINTING carthh\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		hd0=g->hdim[b][0];
		hd1=g->hdim[b][1];
		for(j=0;j<hd1;j++){
			for(i=0;i<hd0;i++){
				hi=i+j*hd0;
				printf("\t%.3g",mm->vcarthh[b][hi]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printcartdiffxx(Momentum* mm,Grid* g){
	int b,i,j,xi;
	int xd0,xd1;
	printf("PRINTING cartdiff\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		xd0=g->xdim[b][0];
		xd1=g->xdim[b][1];
		for(j=0;j<xd1;j++){
			for(i=0;i<xd0;i++){
				xi=i+j*xd0;
				printf("\t%.3g",mm->ucarteta[b][xi]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printuhh(Momentum* mm,Grid* g){
	int b,i,j,hi;
	int hd0,hd1;
	printf("PRINTING uhh\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		hd0=g->hdim[b][0];
		hd1=g->hdim[b][1];
		for(j=0;j<hd1;j++){
			for(i=0;i<hd0;i++){
				hi=i+j*hd0;
				printf("\t%.3g",mm->uhh[b][hi]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printcartvv(Momentum* mm,Grid* g){
	int b,i,j,vi;
	int vd0,vd1;
	printf("PRINTING cartvv\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				printf("\t%.3g",mm->vcartvv[b][vi]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printI1(Momentum* mm,Grid* g){
	int b,i,j,vi;
	int vd0,vd1;
	printf("PRINTING I11_1\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				printf("\t%.3g",mm->I11_1[b][vi]);
			}
			printf("\n");
		}
		printf("\n");
	}
	printf("PRINTING I11_2\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				printf("\t%.3g",mm->I11_2[b][vi]);
			}
			printf("\n");
		}
		printf("\n");
	}
	printf("PRINTING I12_1\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				printf("\t%.3g",mm->I12_1[b][vi]);
			}
			printf("\n");
		}
		printf("\n");
	}
	printf("PRINTING I12_2\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				printf("\t%.3g",mm->I12_2[b][vi]);
			}
			printf("\n");
		}
		printf("\n");
	}
	printf("PRINTING I14_1\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				printf("\t%.3g",mm->I14_1[b][vi]);
			}
			printf("\n");
		}
		printf("\n");
	}
	printf("PRINTING I14_2\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				printf("\t%.3g",mm->I14_2[b][vi]);
			}
			printf("\n");
		}
		printf("\n");
	}
	printf("PRINTING I15_1\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				printf("\t%.3g",mm->I15_1[b][vi]);
			}
			printf("\n");
		}
		printf("\n");
	}
	printf("PRINTING I15_2\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				printf("\t%.3g",mm->I15_2[b][vi]);
			}
			printf("\n");
		}
		printf("\n");
	}
}

void printI2(Momentum* mm,Grid* g){
	int b,i,j,vi;
	int vd0,vd1;
	printf("PRINTING I21_1\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				printf("\t%.3g",mm->I21_1[b][vi]);
			}
			printf("\n");
		}
		printf("\n");
	}
	printf("PRINTING I21_2\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				printf("\t%.3g",mm->I21_2[b][vi]);
			}
			printf("\n");
		}
		printf("\n");
	}
	printf("PRINTING I22_1\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				printf("\t%.3g",mm->I22_1[b][vi]);
			}
			printf("\n");
		}
		printf("\n");
	}
	printf("PRINTING I22_2\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				printf("\t%.3g",mm->I22_2[b][vi]);
			}
			printf("\n");
		}
		printf("\n");
	}
	printf("PRINTING I24_1\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				printf("\t%.3g",mm->I24_1[b][vi]);
			}
			printf("\n");
		}
		printf("\n");
	}
	printf("PRINTING I24_2\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				printf("\t%.3g",mm->I24_2[b][vi]);
			}
			printf("\n");
		}
		printf("\n");
	}
	printf("PRINTING I25_1\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				printf("\t%.3g",mm->I25_1[b][vi]);
			}
			printf("\n");
		}
		printf("\n");
	}
	printf("PRINTING I25_2\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				printf("\t%.3g",mm->I25_2[b][vi]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printucartcc(Momentum* mm,Grid* g){
	int b,i,j,ci;
	int cd0,cd1;
	printf("PRINTING ucartcc\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		cd0=g->cdim[b][0];
		cd1=g->cdim[b][1];
		for(j=0;j<cd1;j++){
			for(i=0;i<cd0;i++){
				ci=i+j*cd0;
				printf("\t%.3g",mm->ucartcc[b][ci]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printvcartcc(Momentum* mm,Grid* g){
	int b,i,j,ci;
	int cd0,cd1;
	printf("PRINTING vcartcc\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		cd0=g->cdim[b][0];
		cd1=g->cdim[b][1];
		for(j=0;j<cd1;j++){
			for(i=0;i<cd0;i++){
				ci=i+j*cd0;
				printf("\t%.3g",mm->vcartcc[b][ci]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printcc(double**cc,Grid*g){
	int b,i,j,ci;
	int cd0,cd1;
	printf("PRINTING cc field:\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		cd0=g->cdim[b][0];
		cd1=g->cdim[b][1];
		for(j=0;j<cd1;j++){
			for(i=0;i<cd0;i++){
				ci=i+j*cd0;
				printf("\t%.3g",cc[b][ci]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printucc(Momentum* mm,Grid* g){
	int b,i,j,ci;
	int cd0,cd1;
	printf("PRINTING ucc\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		cd0=g->cdim[b][0];
		cd1=g->cdim[b][1];
		for(j=0;j<cd1;j++){
			for(i=0;i<cd0;i++){
				ci=i+j*cd0;
				printf("\t%.3g",mm->ucc[b][ci]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printvcc(Momentum* mm,Grid* g){
	int b,i,j,ci;
	int cd0,cd1;
	printf("PRINTING vcc\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		cd0=g->cdim[b][0];
		cd1=g->cdim[b][1];
		for(j=0;j<cd1;j++){
			for(i=0;i<cd0;i++){
				ci=i+j*cd0;
				printf("\t%.3g",mm->vcc[b][ci]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printrhs(LinearSystem* ls,Grid* g){
	int ii,e;
	printf("PRINTING rhs\n");
	e=ls->crows;
	for(ii=0;ii<e;ii++){
		printf("\t%.3g\n",ls->rhs[ii]);
	}
}
void printctrhs(LinearSystem* ls,Grid* g){
	int ii,e;
	printf("PRINTING ctrhs\n");
	e=ls->ccols;
	for(ii=0;ii<e;ii++){
		printf("\t%.3g\n",ls->ctrhs[ii]);
	}
}
void printut(State*s,Momentum* mm,Grid* g){
	int b,i,j,vi;
	int vd0,vd1;
	printf("PRINTING u\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				printf("\t%.3g",s->ut[b][vi]);
			}
			printf("\n");
		}
		printf("\n");
	}
}
void printm(double * m, int rows, int cols, char * name){
	int i,j;
	printf("Printing \"%s\":\n", name);
	for(j=0;j<cols;j++){
		printf("\tCol%d\t",j);
	}
	printf("\n");
	for(j=0;j<rows;j++){
		printf("Row%d:\t",j);
		for(i=0;i<cols;i++){
			printf("%.2e\t",m[j*cols+i]);
		}
		printf("\n");
	}
}
void printmfile(double * m, int rows, int cols, char * name){
	int i,j;
	FILE*file=fopen("out/matrix.txt","w");
	fprintf(file,"Printing \"%s\":\n", name);
	for(j=0;j<cols;j++){
		fprintf(file,"\tCol%d\t",j);
	}
	fprintf(file,"\n");
	for(j=0;j<rows;j++){
		fprintf(file,"Row%d:\t",j);
		for(i=0;i<cols;i++){
			fprintf(file,"%.2e\t",m[j*cols+i]);
		}
		fprintf(file,"\n");
	}
	fclose(file);
}
void printcsrfile(int nn,double*v,int*ri,int*ci){//CSR*csr){
	int i;
	FILE*file=fopen("out/csr.txt","w");
	fprintf(file,"Printing CSR:\n");
	fprintf(file,"Val\t\t\tRow\tCol\n");
	for(i=0;i<nn;i++){
		fprintf(file,"%.2e\t%d\t%d\n",v[i],ri[i],ci[i]);
	}
	//for(ci=0;ci<csr->nn;ci++){
	//	fprintf(file,"%g\t%d\t%d\n",csr->v[ci],csr->ri[ci],csr->ci[ci]);
	//}
	fclose(file);
}
void printsl(SymLap* s,Grid* g){
	int b,ci;
	int cd0,cd3;
	printf("PRINTING u\n");
	printf("sl0\t\tsl1\t\tsl2\t\tsl3\t\tsl4\n");
	for(b=0;b<g->totalblocks;b++){
		printf("Block %d\n",b);
		cd0=g->cdim[b][0];
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			printf("%.2e\t",s->d0[b][ci]);
			if(ci<cd3-1){ printf("%.2e\t",s->d1[b][ci]); }
			if(ci<cd3-cd0+1){ printf("%.2e\t",s->d2[b][ci]); }
			if(ci<cd3-cd0){ printf("%.2e\t",s->d3[b][ci]); }
			if(ci<cd3-cd0-1){ printf("%.2e\t",s->d4[b][ci]); }
			printf("\n");
		}
		printf("\n");
	}
}
void checkcsr(SymLap*sl){
	CSR csr=sl->csr;
	int i;
	printf("v\t\tci\t\tri\n");
	for(i=0;i<csr.nn;i++){
		printf("%.2e\t%d\t%d\n",csr.v[i],csr.ci[i],csr.ri[i]);
	}
}

#endif /* MATRIX_H_ */
