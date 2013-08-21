/*
 * efs.h
 *
 *  Created on: Jul 10, 2013
 *      Author: lordvon
 */

#ifndef EFS_H_
#define EFS_H_

void uncoupled(int edge0,int node0,int *xdim,double *Cmb,int Cmbnodes){
	int hi,vi,xi,ei;
	int vdim0=xdim[0];
	int vdim3=vdim0*(xdim[1]-1);
	int hdim0=xdim[0]-1;
	int hdim3=hdim0*xdim[1];
	int edge,pos,neg;
	ei=edge0;
	for(vi=0;vi<vdim3;vi++){
		xi=vi;
		edge=ei+vi;
		pos=node0+xi+xdim[0];
		neg=node0+xi;
		Cmb[edge*Cmbnodes+neg]=-1;
		Cmb[edge*Cmbnodes+pos]=1;
	}
	ei+=vdim3;
	for(hi=0;hi<hdim3;hi++){
		xi=hi+hi/hdim0;
		edge=ei+hi;
		pos=node0+xi;
		neg=node0+xi+1;
		Cmb[edge*Cmbnodes+neg]=-1;
		Cmb[edge*Cmbnodes+pos]=1;
	}
}
void stitch(int interface,Grid* grid,Interfaces* interfaces,double *Cmb,int Cmbnodes){
	int i,xi1,xi2,rowbase,base1,base2;
	int iv1,iv2;
	int sign;
	int block1=interfaces->block[interface];
	int side=interfaces->side[interface];
	int block2=interfaces->adjacentBlock[interface];

	if(side==0){
		iv1=iv2=1;sign=1;
		base1=0;
		base2=grid->xdim[block2][3]-grid->xdim[block2][0];
	} else if(side==1){
		iv1=grid->xdim[block1][0];iv2=grid->xdim[block2][0];sign=1;
		base1=grid->xdim[block1][0]-1;
		base2=0;
	} else if(side==2){
		iv1=iv2=1;sign=-1;
		base1=grid->xdim[block1][3]-grid->xdim[block1][0];
		base2=0;
	} else if(side==3){
		iv1=grid->xdim[block1][0];iv2=grid->xdim[block2][0];sign=-1;
		base1=0;
		base2=grid->xdim[block2][0]-1;
	}

	for(i=0;i<interfaces->xsize[interface];i++){
		rowbase=(interfaces->rowstart[interface]+i)*Cmbnodes;
		xi1=rowbase+grid->xstart[block1]+base1+i*iv1;
		xi2=rowbase+grid->xstart[block2]+base2+i*iv2;
		Cmb[xi1]=sign;
		Cmb[xi2]=-sign;
	}
}
double *constructCmb(Grid* grid,Interfaces* interfaces){
	int b,i;
	int totalrows=grid->totaledges+interfaces->totalinterfacepoints;
	double *Cmb=(double *)malloc(sizeof(double)*grid->totalnodes*totalrows);
	initv(Cmb,grid->totalnodes*totalrows,0);
	for(b=0;b<grid->totalblocks;b++){
		printf("cmb construction: %d %d\n",grid->estart[b],grid->xstart[b]);
		uncoupled(grid->estart[b],grid->xstart[b],grid->xdim[b],Cmb,grid->totalnodes);
	}
	for(i=0;i<interfaces->totalinterfaces;i++){
		stitch(i,grid,interfaces,Cmb,grid->totalnodes);
	}
	return Cmb;
}


#endif /* EFS_H_ */
