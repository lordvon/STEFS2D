/*
 * cmbops.h
 *
 *  Created on: Jul 17, 2013
 *      Author: lordvon
 */

#ifndef CMBOPS_H_
#define CMBOPS_H_

void posdiag(int rowstart,int colstart,int length,int stride,
		double *s,double *ctrhs){
	//stride must be at least one longer than length for the case of no stride.
	int i;
	ctrhs[rowstart]+=s[colstart];
	for(i=1;i<length;i++){
		if(((i+1)%stride)!=0){
			ctrhs[rowstart+i]+=s[colstart+i];
		}
	}
}
void negdiag(int rowstart,int colstart,int length,int stride,
		double *s,double *ctrhs){
	//stride must be at least one longer than length for the case of no stride.
	int i;
	ctrhs[rowstart]-=s[colstart];
	for(i=1;i<length;i++){
		if(((i+1)%stride)!=0){
			ctrhs[rowstart+i]-=s[colstart+i];
		}
	}
}
void symmnegdiag(int upperrowstart,int uppercolstart,int length,int stride,
		double *s,double *ctrhs){
	negdiag(upperrowstart,uppercolstart,length,stride,s,ctrhs);
	negdiag(uppercolstart,upperrowstart,length,stride,s,ctrhs);
}
void fillmaindiag(Grid* grid,Interfaces* interfaces,LinearSystem* ls){
	int i,ip,j,xi,b,base,xd0,xd1,hside,vside;
	for(b=0;b<grid->totalblocks;b++){
		base=grid->xstart[b];
		xd0=grid->xdim[b][0];
		xd1=grid->xdim[b][1];
		for(j=0;j<xd1;j++){
			for(i=0;i<xd0;i++){
				hside=(j%(xd1-1))==0;
				vside=(i%(xd0-1))==0;
				xi=base+i+j*xd0;
				if( hside & vside ){
					ls->ctcmain[xi]=2;
				} else if( hside | vside){
					ls->ctcmain[xi]=3;
				} else {
					ls->ctcmain[xi]=4;
				}
			}
		}
	}
	for(i=0;i<interfaces->totalinterfaces;i++){
		for(ip=0;ip<interfaces->xsize[i];ip++){
			ls->ctcmain[interfaces->blockx[i][ip]]+=1;
			ls->ctcmain[interfaces->adjacentBlockx[i][ip]]+=1;
		}
	}
}
void ctcmult(double *s, double *newvec,Grid *grid,Interfaces* interfaces,LinearSystem* ls){
	int xi,i,ip,b,xd0,xd3;
	initv(newvec,ls->ccols,0);
	for(xi=0;xi<ls->ccols;xi++){
		newvec[xi]+=s[xi]*ls->ctcmain[xi];
	}
	for(b=0;b<grid->totalblocks;b++){
		xd0=grid->xdim[b][0];
		xd3=grid->xdim[b][3];
		symmnegdiag(grid->xstart[b],grid->xstart[b]+xd0,xd3-xd0,xd3-xd0+1,s,newvec);
		symmnegdiag(grid->xstart[b],grid->xstart[b]+1,xd3-1,xd0,s,newvec);
	}
	for(i=0;i<interfaces->totalinterfaces;i++){
		for(ip=0;ip<interfaces->xsize[i];ip++){
			symmnegdiag(interfaces->blockx[i][ip],interfaces->adjacentBlockx[i][ip],
					1,2,s,newvec);
		}
	}
}
void ctmult(double *s, double *ctrhs,Grid *g,Interfaces* is,LinearSystem* ls){
	int i,ip,b,xd0,vd3,hd0,xs,es,hrow,posorneg;
	initv(ctrhs,ls->ccols,0);
	for(b=0;b<g->totalblocks;b++){
		xs=g->xstart[b];
		es=g->estart[b];
		vd3=g->vdim[b][3];
		hd0=g->hdim[b][0];
		xd0=g->xdim[b][0];
		negdiag(xs,		es,		vd3,vd3+1,s,ctrhs);
		posdiag(xs+xd0,	es,		vd3,vd3+1,s,ctrhs);
		for(hrow=0;hrow<g->hdim[b][1];hrow++){
			posdiag(xs+hrow*hd0+hrow,	es+vd3+hrow*hd0,	hd0,hd0+1,s,ctrhs);
			negdiag(xs+hrow*hd0+hrow+1,	es+vd3+hrow*hd0,	hd0,hd0+1,s,ctrhs);
		}
	}
	for(i=0;i<is->totalinterfaces;i++){
		posorneg=(is->side[i]==0) | (is->side[i]==1);
		for(ip=0;ip<is->xsize[i];ip++){
			es=is->rowstart[i]+ip;
			if(posorneg){
				posdiag(is->blockx[i][ip],es,
						1,2,s,ctrhs);
				negdiag(is->adjacentBlockx[i][ip],es,
						1,2,s,ctrhs);
			} else {
				negdiag(is->blockx[i][ip],es,
						1,2,s,ctrhs);
				posdiag(is->adjacentBlockx[i][ip],es,
						1,2,s,ctrhs);
			}
		}
	}
}
void cmult(double*s,double*cs,Grid*g,Interfaces*is,LinearSystem*ls){
	int i,ip,b,xd0,vd3,hd0,xs,es,hrow,posorneg;
	initv(cs,ls->crows,0);
	for(b=0;b<g->totalblocks;b++){
		xs=g->xstart[b];
		es=g->estart[b];
		vd3=g->vdim[b][3];
		hd0=g->hdim[b][0];
		xd0=g->xdim[b][0];
		negdiag(es,		xs,		vd3,vd3+1,s,cs);
		posdiag(es,		xs+xd0,	vd3,vd3+1,s,cs);
		for(hrow=0;hrow<g->hdim[b][1];hrow++){
			posdiag(es+vd3+hrow*hd0,	xs+hrow*hd0+hrow,		hd0,hd0+1,s,cs);
			negdiag(es+vd3+hrow*hd0,	xs+hrow*hd0+hrow+1,		hd0,hd0+1,s,cs);
		}
	}
	for(i=0;i<is->totalinterfaces;i++){
		posorneg=(is->side[i]==0) | (is->side[i]==1);
		for(ip=0;ip<is->xsize[i];ip++){
			es=is->rowstart[i]+ip;
			if(posorneg){
				posdiag(es,is->blockx[i][ip],
						1,2,s,cs);
				negdiag(es,is->adjacentBlockx[i][ip],
						1,2,s,cs);
			} else {
				negdiag(es,is->blockx[i][ip],
						1,2,s,cs);
				posdiag(es,is->adjacentBlockx[i][ip],
						1,2,s,cs);
			}
		}
	}
}

#endif /* CMBOPS_H_ */
