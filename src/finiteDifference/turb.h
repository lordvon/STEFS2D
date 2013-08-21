/*
 * turb.h
 *
 *  Created on: Jul 21, 2013
 *      Author: lordvon
 */

#ifndef TURB_H_
#define TURB_H_

void CentralDiffCcFieldInterfaces(Interfaces*is,Grid*g,
		double**cc,double**ccxi,double**cceta){
	int i,ip,b1,b2,xi01,xi02,ci01,ci02,ci1,ci2,civ1,civ2,adj1,adj2,s1;
	int cd01,cd02;
	for(i=0;i<is->totalinterfaces;i++){
		b1=is->block[i];
		s1=is->side[i];
		b2=is->adjacentBlock[i];
		cd01=g->cdim[b1][0];
		cd02=g->cdim[b2][0];
		xi01=is->blockx[i][0]-g->xstart[b1];
		xi02=is->adjacentBlockx[i][0]-g->xstart[b2];
		switch(s1){
		case 0: adj1=0;adj2=-cd02; break;
		case 1: adj1=-1;adj2=0; break;
		case 2: adj1=-cd01;adj2=0; break;
		case 3: adj1=0;adj2=-1; break;
		}
		ci01=xi01-xi01/g->xdim[b1][0]+adj1;
		ci02=xi02-xi02/g->xdim[b2][0]+adj2;
		if(s1%2==0){
			civ1=civ2=1;
		} else {
			civ1=cd01;
			civ2=cd02;
		}
		for(ip=0;ip<is->xsize[i]-1;ip++){
			ci1=ci01+ip*civ1;
			ci2=ci02+ip*civ2;
			switch(s1){
			case 0:
				cceta[b1][ci1]=0.5*(cc[b1][ci1+cd01]-	cc[b2][ci2]);
				cceta[b2][ci2]=0.5*(cc[b1][ci1]-		cc[b2][ci2-cd02]);
				break;
			case 2:
				cceta[b1][ci1]=0.5*(cc[b2][ci2]-		cc[b1][ci1-cd01]);
				cceta[b2][ci2]=0.5*(cc[b2][ci2+cd02]-	cc[b1][ci1]);
				break;
			case 1:
				ccxi[b1][ci1]=0.5*(cc[b2][ci2]-cc[b1][ci1-1]);
				ccxi[b2][ci2]=0.5*(cc[b2][ci2+1]-cc[b1][ci1]);
				break;
			case 3:
				ccxi[b1][ci1]=0.5*(cc[b1][ci1+1]-cc[b2][ci2]);
				ccxi[b2][ci2]=0.5*(cc[b1][ci1]-cc[b2][ci2-1]);
				break;
			}
		}
	}
}
void CentralDiffCcField(Interfaces*is,Grid*g,
		char**bcid,double*fixed,
		double**cc,double**ccxi,double**cceta){
	//first order derivative of second order accuracy.
	//zero gradient assumes boundary orthogonality!
	int b,cd0,cd1,ci,i,j,s;
	double interp;
	for(b=0;b<g->totalblocks;b++){
		cd0=g->cdim[b][0];
		cd1=g->cdim[b][1];
		//interior
		for(j=0;j<cd1;j++){
			for(i=1;i<cd0-1;i++){
				ci=i+j*cd0;
				ccxi[b][ci]=0.5*(cc[b][ci+1]-cc[b][ci-1]);
			}
		}
		for(j=1;j<cd1-1;j++){
			for(i=0;i<cd0;i++){
				ci=i+j*cd0;
				cceta[b][ci]=0.5*(cc[b][ci+cd0]-cc[b][ci-cd0]);
			}
		}
		//boundaries.
		for(i=0;i<cd0;i++){
			j=0;s=0;ci=i+j*cd0;
			interp=0.5*(cc[b][ci]+cc[b][ci+cd0]);
			switch(bcid[b][s]){
			case 'w': cceta[b][ci]=interp; break;
			case 'f': cceta[b][ci]=interp-fixed[b]; break;
			case 'z': cceta[b][ci]=0; break;
			}
			j=cd1-1;s=2;ci=i+j*cd0;
			interp=0.5*(cc[b][ci]+cc[b][ci-cd0]);
			switch(bcid[b][s]){
			case 'w': cceta[b][ci]=-interp; break;
			case 'f': cceta[b][ci]=fixed[b]-interp; break;
			case 'z': cceta[b][ci]=0; break;
			}
		}
		for(j=0;j<cd1;j++){
			i=0;s=3;ci=i+j*cd0;
			interp=0.5*(cc[b][ci]+cc[b][ci+1]);
			switch(bcid[b][s]){
			case 'w': ccxi[b][ci]=interp; break;
			case 'f': ccxi[b][ci]=interp-fixed[b]; break;
			case 'z': ccxi[b][ci]=0; break;
			}
			i=cd0-1;s=1;ci=i+j*cd0;
			interp=0.5*(cc[b][ci]+cc[b][ci-1]);
			switch(bcid[b][s]){
			case 'w': ccxi[b][ci]=-interp; break;
			case 'f': ccxi[b][ci]=fixed[b]-interp; break;
			case 'z': ccxi[b][ci]=0; break;
			}
		}
	}
	CentralDiffCcFieldInterfaces(is,g,cc,ccxi,cceta);
}
void CentralDiffCcFieldFOB(Interfaces*is,Grid*g,
		double**cc,double**ccxi,double**cceta){
	//FOB: first-order accurate boundary.
	//first order derivative of second order accuracy.
	//currently only for wall and fixed.
	int b,cd0,cd1,ci,i,j;
	for(b=0;b<g->totalblocks;b++){
		cd0=g->cdim[b][0];
		cd1=g->cdim[b][1];
		//interior
		for(j=0;j<cd1;j++){
			for(i=1;i<cd0-1;i++){
				ci=i+j*cd0;
				ccxi[b][ci]=0.5*(cc[b][ci+1]-cc[b][ci-1]);
			}
		}
		for(j=1;j<cd1-1;j++){
			for(i=0;i<cd0;i++){
				ci=i+j*cd0;
				cceta[b][ci]=0.5*(cc[b][ci+cd0]-cc[b][ci-cd0]);
			}
		}
		//boundaries.
		for(i=0;i<cd0;i++){
			j=0;ci=i+j*cd0;
			cceta[b][ci]=cc[b][ci+cd0]-cc[b][ci];
			j=cd1-1;ci=i+j*cd0;
			cceta[b][ci]=cc[b][ci]-cc[b][ci-cd0];
		}
		for(j=0;j<cd1;j++){
			i=0;ci=i+j*cd0;
			ccxi[b][ci]=cc[b][ci+1]-cc[b][ci];
			i=cd0-1;ci=i+j*cd0;
			ccxi[b][ci]=cc[b][ci]-cc[b][ci-1];
		}
	}
	CentralDiffCcFieldInterfaces(is,g,cc,ccxi,cceta);
}
void fillCartDerivativeCc(Grid*g,double **dx,double**dy,double**dxi,double**deta){
	int b,ci,cd3;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			dx[b][ci]=	g->xixcc[b][ci]*dxi[b][ci]+
					g->etaxcc[b][ci]*deta[b][ci];
			dy[b][ci]=	g->xiycc[b][ci]*dxi[b][ci]+
					g->etaycc[b][ci]*deta[b][ci];
		}
	}
}
void fillCartDerivativeXCc(Grid*g,double **dx,double**dxi,double**deta){
	int b,ci,cd3;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			dx[b][ci]=	g->xixcc[b][ci]*dxi[b][ci]+
					g->etaxcc[b][ci]*deta[b][ci];
		}
	}
}
void fillCartDerivativeYCc(Grid*g,double**dy,double**dxi,double**deta){
	int b,ci,cd3;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			dy[b][ci]=	g->xiycc[b][ci]*dxi[b][ci]+
					g->etaycc[b][ci]*deta[b][ci];
		}
	}
}

#endif /* TURB_H_ */
