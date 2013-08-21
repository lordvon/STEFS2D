/*
 * solution.h
 *
 *  Created on: Jul 25, 2013
 *      Author: lordvon
 */

#ifndef SOLUTION_H_
#define SOLUTION_H_

void subWriteSoln(FILE *file,double *soln,int b,Grid* g){
	int i,j,k,xi;
	for (k=0;k<g->xdim[b][2];k++){
		for (j=0;j<g->xdim[b][1];j++){
			for (i=0;i<g->xdim[b][0];i++){
				xi=i+j*g->xdim[b][0];
				fprintf(file,"\t%f",soln[xi]);
			}
			fprintf(file,"\n");
		}
		fprintf(file,"\n");
	}
}
void writeMultiBlockStateSolution(char *name,Grid* g,State*s,Momentum* mm,
		Interfaces* is,BoundaryConditions* bc,LinearSystem* ls){
	//Writes out formatted (ASCII) solution in whole, multi-block PLOT3D format.
	FILE * file = fopen(name, "w");
	int b;
	//first convert latest state to cartesian form.
	//update(s,g,ls,bc);
	fillxx(s,mm,g,bc,is);
	fillCartxx(g,mm);
	//Information headers.
	fprintf(file,"%d\n",g->totalblocks);
	for (b=0;b<g->totalblocks;b++){
		fprintf(file,"%d\t%d\t%d\n",g->xdim[b][0],g->xdim[b][1],g->xdim[b][2]);
	}
	//dummy field definitions.
	double **wv;
	double **e;
	wv=(double **)malloc(sizeof(double)*g->totalblocks);
	e=(double **)malloc(sizeof(double)*g->totalblocks);
	for (b=0;b<g->totalblocks;b++){
		wv[b]=(double *)malloc(sizeof(double)*g->xdim[b][3]);
		initv(wv[b],g->xdim[b][3],0);
		e[b]=(double *)malloc(sizeof(double)*g->xdim[b][3]);
		initv(e[b],g->xdim[b][3],0);
	}
	//Write fields.
	double mach=0, alpha=0, reyn=0, time=0;
	for (b=0;b<g->totalblocks;b++){
		fprintf(file,"%f\t%f\t%f\t%f\n",mach,alpha,reyn,time);
		subWriteSoln(file,e[b],b,g);
		subWriteSoln(file,mm->ucartxx[b],b,g);
		subWriteSoln(file,mm->vcartxx[b],b,g);
		subWriteSoln(file,wv[b],b,g);
		subWriteSoln(file,e[b],b,g);
	}
	fclose(file);
	//free
	for (b=0;b<g->totalblocks;b++){
		free(wv[b]);
		free(e[b]);
	}
	free(wv);
	free(e);
}
void writeMultiBlockCustomSolution(char *name,Grid* g,double **xxsolution){
	//Writes out formatted (ASCII) solution in whole, multi-block PLOT3D format.
	FILE * file = fopen(name, "w");
	int b;
	//Information headers.
	fprintf(file,"%d\n",g->totalblocks);
	for (b=0;b<g->totalblocks;b++){
		fprintf(file,"%d\t%d\t%d\n",g->xdim[b][0],g->xdim[b][1],g->xdim[b][2]);
	}
	//dummy field definitions.
	double **zero=(double **)malloc(sizeof(double)*g->totalblocks);
	for (b=0;b<g->totalblocks;b++){
		zero[b]=(double *)malloc(sizeof(double)*g->xdim[b][3]);
		initv(zero[b],g->xdim[b][3],0);
	}
	//Write fields.
	double mach=0, alpha=0, reyn=0, time=0;
	for (b=0;b<g->totalblocks;b++){
		fprintf(file,"%f\t%f\t%f\t%f\n",mach,alpha,reyn,time);
		subWriteSoln(file,xxsolution[b],b,g);
		subWriteSoln(file,zero[b],b,g);
		subWriteSoln(file,zero[b],b,g);
		subWriteSoln(file,zero[b],b,g);
		subWriteSoln(file,zero[b],b,g);
	}
	fclose(file);
	//free
	for (b=0;b<g->totalblocks;b++){
		free(zero[b]);
	}
	free(zero);
}

#endif /* SOLUTION_H_ */
