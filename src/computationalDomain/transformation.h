/*
 * transformation.h
 *
 *  Created on: Jul 18, 2013
 *      Author: lordvon
 */

#ifndef TRANSFORMATION_H_
#define TRANSFORMATION_H_

void computeCovariant(Grid* g){
	int b,xi,vi,hi;
	int hd3,hd0,xd0,vd3;
	for(b=0;b<g->totalblocks;b++){
		hd0=g->hdim[b][0];
		hd3=g->hdim[b][3];
		vd3=g->vdim[b][3];
		xd0=g->xdim[b][0];
		for(hi=0;hi<hd3;hi++){
			xi=hi+hi/hd0;
			g->a11[b][hi]=g->x[b][xi+1]-g->x[b][xi];
			g->a12[b][hi]=g->y[b][xi+1]-g->y[b][xi];
		}
		for(vi=0;vi<vd3;vi++){
			xi=vi;
			g->a21[b][vi]=g->x[b][xi+xd0]-g->x[b][xi];
			g->a22[b][vi]=g->y[b][xi+xd0]-g->y[b][xi];
		}
	}
}
void fillCovariantcc(Grid* g){
	int b,ci,hi,vi;
	int cd3,cd0,hd0;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		cd0=g->cdim[b][0];
		hd0=g->hdim[b][0];
		for(ci=0;ci<cd3;ci++){
			hi=ci;
			vi=ci+ci/cd0;
			g->a11cc[b][ci]=0.5*(g->a11[b][hi]+g->a11[b][hi+hd0]);
			g->a12cc[b][ci]=0.5*(g->a12[b][hi]+g->a12[b][hi+hd0]);
			g->a21cc[b][ci]=0.5*(g->a21[b][vi]+g->a21[b][vi+1]);
			g->a22cc[b][ci]=0.5*(g->a22[b][vi]+g->a22[b][vi+1]);
		}
	}
}
void fillCovariantxx(Grid* g,BoundaryConditions* bc){
	int b,i,j,xi,hi,vi;
	int xd1,xd0,vd0;
	for(b=0;b<g->totalblocks;b++){
		xd1=g->xdim[b][1];
		xd0=g->xdim[b][0];
		vd0=g->vdim[b][0];
		for(j=0;j<xd1;j++){
			for(i=0;i<xd0;i++){
				xi=i+j*xd0;
				hi=xi-j;
				vi=xi;
				if(i==0){
					g->a11xx[b][xi]=-0.5*(3*g->x[b][xi]-4*g->x[b][xi+1]+g->x[b][xi+2]);
					g->a12xx[b][xi]=-0.5*(3*g->y[b][xi]-4*g->y[b][xi+1]+g->y[b][xi+2]);
				} else if(i==xd0-1){
					g->a11xx[b][xi]=0.5*(3*g->x[b][xi]-4*g->x[b][xi-1]+g->x[b][xi-2]);
					g->a12xx[b][xi]=0.5*(3*g->y[b][xi]-4*g->y[b][xi-1]+g->y[b][xi-2]);
				} else {
					g->a11xx[b][xi]=0.5*(g->a11[b][hi]+g->a11[b][hi-1]);
					g->a12xx[b][xi]=0.5*(g->a12[b][hi]+g->a12[b][hi-1]);
				}
				if((j==0)){
					g->a21xx[b][xi]=-0.5*(3*g->x[b][xi]-4*g->x[b][xi+xd0]+g->x[b][xi+2*xd0]);
					g->a22xx[b][xi]=-0.5*(3*g->y[b][xi]-4*g->y[b][xi+xd0]+g->y[b][xi+2*xd0]);
				} else if(j==xd1-1){
					g->a21xx[b][xi]=0.5*(3*g->x[b][xi]-4*g->x[b][xi-xd0]+g->x[b][xi-2*xd0]);
					g->a22xx[b][xi]=0.5*(3*g->y[b][xi]-4*g->y[b][xi-xd0]+g->y[b][xi-2*xd0]);
				} else {
					g->a21xx[b][xi]=0.5*(g->a21[b][vi]+g->a21[b][vi-vd0]);
					g->a22xx[b][xi]=0.5*(g->a22[b][vi]+g->a22[b][vi-vd0]);
				}
			}
		}
	}
}
void fillCovariantxxinterfaces(Grid* g,Interfaces* is){
	int i,ip,side;
	int b1,b2,xi1,xi2,vi1,vi2,hi1,hi2;
	for(i=0;i<is->totalinterfaces;i++){
		side=is->side[i];
		b1=is->block[i];
		b2=is->adjacentBlock[i];
		for(ip=0;ip<is->xsize[i];ip++){
			xi1=is->blockx[i][ip]-g->xstart[b1];
			xi2=is->adjacentBlockx[i][ip]-g->xstart[b2];
			vi1=xi1;
			vi2=xi2;
			hi1=xi1-xi1/g->xdim[b1][0];
			hi2=xi2-xi2/g->xdim[b2][0];
			switch(side){
			case 0:
				g->a21xx[b1][xi1]=g->a21xx[b2][xi2]=0.5*(
						g->a21[b1][vi1]+
						g->a21[b2][vi2-g->vdim[b2][0]]);
				g->a22xx[b1][xi1]=g->a22xx[b2][xi2]=0.5*(
						g->a22[b1][vi1]+
						g->a22[b2][vi2-g->vdim[b2][0]]);
				break;
			case 1:
				g->a11xx[b1][xi1]=g->a11xx[b2][xi2]=0.5*(
						g->a11[b1][hi1-1]+
						g->a11[b2][hi2]);
				g->a12xx[b1][xi1]=g->a12xx[b2][xi2]=0.5*(
						g->a12[b1][hi1-1]+
						g->a12[b2][hi2]);
				break;
			case 2:
				g->a21xx[b1][xi1]=g->a21xx[b2][xi2]=0.5*(
						g->a21[b1][vi1-g->vdim[b1][0]]+
						g->a21[b2][vi2]);
				g->a22xx[b1][xi1]=g->a22xx[b2][xi2]=0.5*(
						g->a22[b1][vi1-g->vdim[b1][0]]+
						g->a22[b2][vi2]);
				break;
			case 3:
				g->a11xx[b1][xi1]=g->a11xx[b2][xi2]=0.5*(
						g->a11[b1][hi1]+
						g->a11[b2][hi2-1]);
				g->a12xx[b1][xi1]=g->a12xx[b2][xi2]=0.5*(
						g->a12[b1][hi1]+
						g->a12[b2][hi2-1]);
				break;
			}
		}
	}
}
void fillCovarianthh(Grid* g,BoundaryConditions* bc){
	int b,i,j,hi,vi;
	int hd1,hd0,vd0;
	for(b=0;b<g->totalblocks;b++){
		hd1=g->hdim[b][1];
		hd0=g->hdim[b][0];
		vd0=g->vdim[b][0];
		for(j=0;j<hd1;j++){
			for(i=0;i<hd0;i++){
				hi=i+j*hd0;
				vi=hi+j;
				if(j==0){
					g->a21hh[b][hi]=0.5*(g->a21[b][vi]+g->a21[b][vi+1]);
					g->a22hh[b][hi]=0.5*(g->a22[b][vi]+g->a22[b][vi+1]);
				} else if(j==hd1-1){
					g->a21hh[b][hi]=0.5*(g->a21[b][vi-vd0]+g->a21[b][vi+1-vd0]);
					g->a22hh[b][hi]=0.5*(g->a22[b][vi-vd0]+g->a22[b][vi+1-vd0]);
				} else {
					g->a21hh[b][hi]=0.25*(g->a21[b][vi]+g->a21[b][vi+1]+g->a21[b][vi-vd0]+g->a21[b][vi+1-vd0]);
					g->a22hh[b][hi]=0.25*(g->a22[b][vi]+g->a22[b][vi+1]+g->a22[b][vi-vd0]+g->a22[b][vi+1-vd0]);
				}
			}
		}
	}
}
void fillCovarianthhinterfaces(Grid* g,Interfaces* is){
	int i,ip,side;
	int b1,b2,xi1,xi2,vi1,vi2,hi1,hi2;
	for(i=0;i<is->totalinterfaces;i++){
		side=is->side[i];
		b1=is->block[i];
		b2=is->adjacentBlock[i];
		for(ip=0;ip<is->xsize[i]-1;ip++){
			xi1=is->blockx[i][ip]-g->xstart[b1];
			xi2=is->adjacentBlockx[i][ip]-g->xstart[b2];
			vi1=xi1;
			vi2=xi2;
			hi1=xi1-xi1/g->xdim[b1][0];
			hi2=xi2-xi2/g->xdim[b2][0];
			switch(side){
			case 0:
				g->a21hh[b1][hi1]=g->a21hh[b2][hi2]=0.25*(
						g->a21[b1][vi1]+
						g->a21[b1][vi1+1]+
						g->a21[b2][vi2-g->vdim[b2][0]]+
						g->a21[b2][vi2-g->vdim[b2][0]+1]);
				g->a22hh[b1][hi1]=g->a22hh[b2][hi2]=0.25*(
						g->a22[b1][vi1]+
						g->a22[b1][vi1+1]+
						g->a22[b2][vi2-g->vdim[b2][0]]+
						g->a22[b2][vi2-g->vdim[b2][0]+1]);
				break;
			case 2:
				g->a21hh[b1][hi1]=g->a21hh[b2][hi2]=0.25*(
						g->a21[b1][vi1-g->vdim[b1][0]]+
						g->a21[b1][vi1-g->vdim[b1][0]+1]+
						g->a21[b2][vi2]+
						g->a21[b2][vi2+1]);
				g->a22hh[b1][hi1]=g->a22hh[b2][hi2]=0.25*(
						g->a22[b1][vi1-g->vdim[b1][0]]+
						g->a22[b1][vi1-g->vdim[b1][0]+1]+
						g->a22[b2][vi2]+
						g->a22[b2][vi2+1]);
				break;
			}
		}
	}
}
void fillCovariantvv(Grid* g,BoundaryConditions* bc){
	int b,i,j,hi,vi;
	int vd1,vd0,hd0;
	for(b=0;b<g->totalblocks;b++){
		vd1=g->vdim[b][1];
		vd0=g->vdim[b][0];
		hd0=g->hdim[b][0];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				hi=vi-j;
				if(i==0){
					g->a11vv[b][vi]=0.5*(g->a11[b][hi]+g->a11[b][hi+hd0]);
					g->a12vv[b][vi]=0.5*(g->a12[b][hi]+g->a12[b][hi+hd0]);
				} else if(i==vd0-1){
					g->a11vv[b][vi]=0.5*(g->a11[b][hi-1]+g->a11[b][hi-1+hd0]);
					g->a12vv[b][vi]=0.5*(g->a12[b][hi-1]+g->a12[b][hi-1+hd0]);
				} else {
					g->a11vv[b][vi]=0.25*(g->a11[b][hi]+g->a11[b][hi-1]+g->a11[b][hi+hd0]+g->a11[b][hi-1+hd0]);
					g->a12vv[b][vi]=0.25*(g->a12[b][hi]+g->a12[b][hi-1]+g->a12[b][hi+hd0]+g->a12[b][hi-1+hd0]);
				}
			}
		}
	}
}
void fillCovariantvvinterfaces(Grid* g,Interfaces* is){
	int i,ip,side;
	int b1,b2,xi1,xi2,vi1,vi2,hi1,hi2;
	for(i=0;i<is->totalinterfaces;i++){
		side=is->side[i];
		b1=is->block[i];
		b2=is->adjacentBlock[i];
		for(ip=0;ip<is->xsize[i]-1;ip++){
			xi1=is->blockx[i][ip]-g->xstart[b1];
			xi2=is->adjacentBlockx[i][ip]-g->xstart[b2];
			vi1=xi1;
			vi2=xi2;
			hi1=xi1-xi1/g->xdim[b1][0];
			hi2=xi2-xi2/g->xdim[b2][0];
			switch(side){
			case 1:
				g->a11vv[b1][vi1]=g->a11vv[b2][vi2]=0.25*(
						g->a11[b1][hi1-1]+
						g->a11[b1][hi1+g->hdim[b1][0]-1]+
						g->a11[b2][hi2]+
						g->a11[b2][hi2+g->hdim[b2][0]]);
				g->a12vv[b1][vi1]=g->a12vv[b2][vi2]=0.25*(
						g->a12[b1][hi1-1]+
						g->a12[b1][hi1+g->hdim[b1][0]-1]+
						g->a12[b2][hi2]+
						g->a12[b2][hi2+g->hdim[b2][0]]);
				break;
			case 3:
				g->a11vv[b1][vi1]=g->a11vv[b2][vi2]=0.25*(
						g->a11[b1][hi1]+
						g->a11[b1][hi1+g->hdim[b1][0]]+
						g->a11[b2][hi2-1]+
						g->a11[b2][hi2+g->hdim[b2][0]-1]);
				g->a12vv[b1][vi1]=g->a12vv[b2][vi2]=0.25*(
						g->a12[b1][hi1]+
						g->a12[b1][hi1+g->hdim[b1][0]]+
						g->a12[b2][hi2-1]+
						g->a12[b2][hi2+g->hdim[b2][0]-1]);
				break;
			}
		}
	}
}
void fillJcc(Grid* g){
	int b,ci;
	int cd3;
	double a1[2],a2[2];
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			a1[0]=g->a11cc[b][ci];
			a1[1]=g->a12cc[b][ci];
			a2[0]=g->a21cc[b][ci];
			a2[1]=g->a22cc[b][ci];
			g->Jcc[b][ci]=cross(a1,a2);
		}
	}
}
void fillJxx(Grid* g){
	int b,xi;
	int xd3;
	double a1[2],a2[2];
	for(b=0;b<g->totalblocks;b++){
		xd3=g->xdim[b][3];
		for(xi=0;xi<xd3;xi++){
			a1[0]=g->a11xx[b][xi];
			a1[1]=g->a12xx[b][xi];
			a2[0]=g->a21xx[b][xi];
			a2[1]=g->a22xx[b][xi];
			g->Jxx[b][xi]=cross(a1,a2);
		}
	}
}
void fillJvv(Grid* g){
	int b,vi,vd3;
	double a1[2],a2[2];
	for(b=0;b<g->totalblocks;b++){
		vd3=g->vdim[b][3];
		for(vi=0;vi<vd3;vi++){
			a1[0]=g->a11vv[b][vi];
			a1[1]=g->a12vv[b][vi];
			a2[0]=g->a21[b][vi];
			a2[1]=g->a22[b][vi];
			g->Jvv[b][vi]=cross(a1,a2);
		}
	}
}
void fillJhh(Grid* g){
	int b,hi,hd3;
	double a1[2],a2[2];
	for(b=0;b<g->totalblocks;b++){
		hd3=g->hdim[b][3];
		for(hi=0;hi<hd3;hi++){
			a1[0]=g->a11[b][hi];
			a1[1]=g->a12[b][hi];
			a2[0]=g->a21hh[b][hi];
			a2[1]=g->a22hh[b][hi];
			g->Jhh[b][hi]=cross(a1,a2);
		}
	}
}
void fillContravariantcc(Grid* g){
	//need covariants (interpolated here), and Jacobians at cell center.
	int b,ci,cd3;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			g->xixcc[b][ci]=g->a22cc[b][ci]/g->Jcc[b][ci];
			g->xiycc[b][ci]=-g->a21cc[b][ci]/g->Jcc[b][ci];
			g->etaxcc[b][ci]=-g->a12cc[b][ci]/g->Jcc[b][ci];
			g->etaycc[b][ci]=g->a11cc[b][ci]/g->Jcc[b][ci];
		}
	}
}
void fillContravariantxx(Grid* g){
	int b,xi,xd3;
	for(b=0;b<g->totalblocks;b++){
		xd3=g->xdim[b][3];
		for(xi=0;xi<xd3;xi++){
			g->xixxx[b][xi]=g->a22xx[b][xi]/g->Jxx[b][xi];
			g->xiyxx[b][xi]=-g->a21xx[b][xi]/g->Jxx[b][xi];
			g->etaxxx[b][xi]=-g->a12xx[b][xi]/g->Jxx[b][xi];
			g->etayxx[b][xi]=g->a11xx[b][xi]/g->Jxx[b][xi];
		}
	}
}
void computevareas(Grid* g){
	int i,j;
	int vi,hi,ci;
	int b,vd0,vd1;
	double a1[2],b1[2],a2[2],b2[2];
	double area1,area2;
	for(b=0;b<g->totalblocks;b++){
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for (j=0;j<vd1;j++){
			for (i=1;i<vd0-1;i++){
				vi=i+j*vd0;
				hi=ci=vi-j;
				a1[0]=0.5*g->a11[b][hi-1]+g->a21[b][vi];
				a1[1]=0.5*g->a12[b][hi-1]+g->a22[b][vi];
				b1[0]=-0.5*g->a11[b][hi-1]+g->a21cc[b][ci-1];
				b1[1]=-0.5*g->a12[b][hi-1]+g->a22cc[b][ci-1];
				a2[0]=0.5*g->a11[b][hi]+g->a21cc[b][ci];
				a2[1]=0.5*g->a12[b][hi]+g->a22cc[b][ci];
				b2[0]=-0.5*g->a11[b][hi]+g->a21[b][vi];
				b2[1]=-0.5*g->a12[b][hi]+g->a22[b][vi];
				area1=0.5*cross(a1,b1);
				area2=0.5*cross(a2,b2);
				if((area1<0) | (area2<0)){ printf("Negative varea detected!"); }
				g->vareas[b][vi]=area1+area2;
			}
		}
		for (j=0;j<vd1;j++){
			i=0;vi=i+j*vd0;
			g->vareas[b][vi]=g->Jvv[b][vi];
			i=vd0-1;vi=i+j*vd0;
			g->vareas[b][vi]=g->Jvv[b][vi];
		}
	}
}
void computevareasinterfaces(Grid* g,Interfaces* is){
	int i,ip,side;
	int b1,b2,xi1,xi2,vi1,vi2,hi1,hi2,ci1,ci2;
	double a[2],b[2],area1,area2;
	for(i=0;i<is->totalinterfaces;i++){
		side=is->side[i];
		b1=is->block[i];
		b2=is->adjacentBlock[i];
		for(ip=0;ip<is->xsize[i]-1;ip++){
			vi1=xi1=is->blockx[i][ip]-g->xstart[b1];
			vi2=xi2=is->adjacentBlockx[i][ip]-g->xstart[b2];
			hi1=ci1=xi1-xi1/g->xdim[b1][0];
			hi2=ci2=xi2-xi2/g->xdim[b2][0];
			switch(side){
			case 1:
				a[0]=0.5*g->a11[b1][hi1-1]+g->a21[b1][vi1];
				a[1]=0.5*g->a12[b1][hi1-1]+g->a22[b1][vi1];
				b[0]=-0.5*g->a11[b1][hi1-1]+g->a21cc[b1][ci1-1];
				b[1]=-0.5*g->a12[b1][hi1-1]+g->a22cc[b1][ci1-1];
				area1=0.5*cross(a,b);
				a[0]=0.5*g->a11[b2][hi2]+g->a21cc[b2][ci2];
				a[1]=0.5*g->a12[b2][hi2]+g->a22cc[b2][ci2];
				b[0]=-0.5*g->a11[b2][hi2]+g->a21[b2][vi2];
				b[1]=-0.5*g->a12[b2][hi2]+g->a22[b2][vi2];
				area2=0.5*cross(a,b);
				if((area1<0) | (area2<0)){ printf("Negative varea detected!"); }
				g->vareas[b1][vi1]=g->vareas[b2][vi2]=area1+area2;
				break;
			case 3:
				a[0]=0.5*g->a11[b2][hi2-1]+g->a21[b2][vi2];
				a[1]=0.5*g->a12[b2][hi2-1]+g->a22[b2][vi2];
				b[0]=-0.5*g->a11[b2][hi2-1]+g->a21cc[b2][ci2-1];
				b[1]=-0.5*g->a12[b2][hi2-1]+g->a22cc[b2][ci2-1];
				area2=0.5*cross(a,b);
				a[0]=0.5*g->a11[b1][hi1]+g->a21cc[b1][ci1];
				a[1]=0.5*g->a12[b1][hi1]+g->a22cc[b1][ci1];
				b[0]=-0.5*g->a11[b1][hi1]+g->a21[b1][vi1];
				b[1]=-0.5*g->a12[b1][hi1]+g->a22[b1][vi1];
				area1=0.5*cross(a,b);
				if((area1<0) | (area2<0)){ printf("Negative varea detected!"); }
				g->vareas[b1][vi1]=g->vareas[b2][vi2]=area1+area2;
				break;
			}
		}
	}
}
void computehareas(Grid* g){
	int b,i,j;
	int hi,vi,ci;
	int hd1,hd0,cd0,vd0;
	double a1[2],b1[2],a2[2],b2[2];
	double area1,area2;
	for(b=0;b<g->totalblocks;b++){
		hd0=g->hdim[b][0];
		hd1=g->hdim[b][1];
		cd0=g->cdim[b][0];
		vd0=g->vdim[b][0];
		for (j=1;j<hd1-1;j++){
			for (i=0;i<hd0;i++){
				hi=ci=i+j*hd0;
				vi=hi+j;
				a1[0]=g->a11[b][hi]+0.5*g->a21[b][vi+1];
				a1[1]=g->a12[b][hi]+0.5*g->a22[b][vi+1];
				b1[0]=-g->a11[b][hi]+0.5*g->a21[b][vi];
				b1[1]=-g->a12[b][hi]+0.5*g->a22[b][vi];
				a2[0]=g->a11cc[b][ci-cd0]+0.5*g->a21[b][vi+1-vd0];
				a2[1]=g->a12cc[b][ci-cd0]+0.5*g->a22[b][vi+1-vd0];
				b2[0]=-g->a11cc[b][ci-cd0]+0.5*g->a21[b][vi-vd0];
				b2[1]=-g->a12cc[b][ci-cd0]+0.5*g->a22[b][vi-vd0];
				area1=0.5*cross(a1,b1);
				area2=0.5*cross(a2,b2);
				if((area1<0) | (area2<0)){ printf("Negative harea detected!"); }
				g->hareas[b][hi]=area1+area2;
			}
		}
		for (i=0;i<hd0;i++){
			j=0;hi=i+j*hd0;
			g->hareas[b][hi]=g->Jhh[b][hi];
			j=hd1-1;hi=i+j*hd0;
			g->hareas[b][hi]=g->Jhh[b][hi];
		}
	}
}
void computehareasinterfaces(Grid* g,Interfaces* is){
	int i,ip,side;
	int b1,b2,xi1,xi2,vi1,vi2,hi1,hi2,ci1,ci2;
	int cd0,vd0;
	double a[2],b[2],area1,area2;
	for(i=0;i<is->totalinterfaces;i++){
		side=is->side[i];
		b1=is->block[i];
		b2=is->adjacentBlock[i];
		for(ip=0;ip<is->xsize[i]-1;ip++){
			vi1=xi1=is->blockx[i][ip]-g->xstart[b1];
			vi2=xi2=is->adjacentBlockx[i][ip]-g->xstart[b2];
			hi1=ci1=xi1-xi1/g->xdim[b1][0];
			hi2=ci2=xi2-xi2/g->xdim[b2][0];
			switch(side){
			case 0:
				a[0]=g->a11[b1][hi1]+0.5*g->a21[b1][vi1+1];
				a[1]=g->a12[b1][hi1]+0.5*g->a22[b1][vi1+1];
				b[0]=-g->a11[b1][hi1]+0.5*g->a21[b1][vi1];
				b[1]=-g->a12[b1][hi1]+0.5*g->a22[b1][vi1];
				area1=0.5*cross(a,b);
				vd0=g->vdim[b2][0];
				cd0=g->cdim[b2][0];
				a[0]=g->a11cc[b2][ci2-cd0]+0.5*g->a21[b2][vi2+1-vd0];
				a[1]=g->a12cc[b2][ci2-cd0]+0.5*g->a22[b2][vi2+1-vd0];
				b[0]=-g->a11cc[b2][ci2-cd0]+0.5*g->a21[b2][vi2-vd0];
				b[1]=-g->a12cc[b2][ci2-cd0]+0.5*g->a22[b2][vi2-vd0];
				area2=0.5*cross(a,b);
				if((area1<0) | (area2<0)){ printf("Negative harea detected!"); }
				g->hareas[b1][hi1]=g->hareas[b2][hi2]=area1+area2;
				break;
			case 2:
				a[0]=g->a11[b2][hi2]+0.5*g->a21[b2][vi2+1];
				a[1]=g->a12[b2][hi2]+0.5*g->a22[b2][vi2+1];
				b[0]=-g->a11[b2][hi2]+0.5*g->a21[b2][vi2];
				b[1]=-g->a12[b2][hi2]+0.5*g->a22[b2][vi2];
				area2=0.5*cross(a,b);
				vd0=g->vdim[b1][0];
				cd0=g->cdim[b1][0];
				a[0]=g->a11cc[b1][ci1-cd0]+0.5*g->a21[b1][vi1+1-vd0];
				a[1]=g->a12cc[b1][ci1-cd0]+0.5*g->a22[b1][vi1+1-vd0];
				b[0]=-g->a11cc[b1][ci1-cd0]+0.5*g->a21[b1][vi1-vd0];
				b[1]=-g->a12cc[b1][ci1-cd0]+0.5*g->a22[b1][vi1-vd0];
				area1=0.5*cross(a,b);
				if((area1<0) | (area2<0)){ printf("Negative harea detected!"); }
				g->hareas[b1][hi1]=g->hareas[b2][hi2]=area1+area2;
				break;
			}
		}
	}
}
void fillcovc(Grid* g,BoundaryConditions* bc){
	int b,i,j,xi,i1,i2;
	int xd1,xd0;
	double factor;
	for(b=0;b<g->totalblocks;b++){
		xd0=g->xdim[b][0];
		xd1=g->xdim[b][1];
		for(j=0;j<xd1;j++){
			for(i=0;i<xd0;i++){
				xi=i+j*xd0;
				if(i==0){
					i1=xi;
					i2=xi+1;
					factor=1;
				} else if(i==xd0-1){
					i1=xi-1;
					i2=xi;
					factor=1;
				} else {
					i1=xi-1;
					i2=xi+1;
					factor=0.5;
				}
				g->covc11[b][xi]=
						factor*(g->x[b][i2]-g->x[b][i1]);
				g->covc12[b][xi]=
						factor*(g->y[b][i2]-g->y[b][i1]);
				if(j==0){
					i1=xi;
					i2=xi+xd0;
					factor=1;
				} else if(j==xd1-1){
					i1=xi-xd0;
					i2=xi;
					factor=1;
				} else {
					i1=xi-xd0;
					i2=xi+xd0;
					factor=0.5;
				}
				g->covc21[b][xi]=factor*(g->x[b][i2]-g->x[b][i1]);
				g->covc22[b][xi]=factor*(g->y[b][i2]-g->y[b][i1]);
			}
		}
	}
}
void fillcovcinterfaces(Grid* g,Interfaces* is){
	int i,ip,side;
	int b1,b2,xi1,xi2;
	int xd01,xd02;
	for(i=0;i<is->totalinterfaces;i++){
		side=is->side[i];
		b1=is->block[i];
		b2=is->adjacentBlock[i];
		for(ip=0;ip<is->xsize[i];ip++){
			xi1=is->blockx[i][ip]-g->xstart[b1];
			xi2=is->adjacentBlockx[i][ip]-g->xstart[b2];
			xd01=g->xdim[b1][0];
			xd02=g->xdim[b2][0];
			switch(side){
			case 0:
				g->covc21[b1][xi1]=g->covc21[b2][xi2]=
						0.5*(g->x[b1][xi1+xd01]-g->x[b2][xi2-xd02]);
				g->covc22[b1][xi1]=g->covc22[b2][xi2]=
						0.5*(g->y[b1][xi1+xd01]-g->y[b2][xi2-xd02]);
				break;
			case 1:
				g->covc11[b1][xi1]=g->covc11[b2][xi2]=
						0.5*(g->x[b2][xi2+1]-g->x[b1][xi1-1]);
				g->covc12[b1][xi1]=g->covc12[b2][xi2]=
						0.5*(g->y[b2][xi2+1]-g->y[b1][xi1-1]);
				break;
			case 2:
				g->covc21[b1][xi1]=g->covc21[b2][xi2]=
						0.5*(g->x[b2][xi2+xd02]-g->x[b1][xi1-xd01]);
				g->covc22[b1][xi1]=g->covc22[b2][xi2]=
						0.5*(g->y[b2][xi2+xd02]-g->y[b1][xi1-xd01]);
				break;
			case 3:
				g->covc11[b1][xi1]=g->covc11[b2][xi2]=
						0.5*(g->x[b1][xi1+1]-g->x[b2][xi2-1]);
				g->covc12[b1][xi1]=g->covc12[b2][xi2]=
						0.5*(g->y[b1][xi1+1]-g->y[b2][xi2-1]);
				break;
			}
		}
	}
}
void fillcontrac(Grid* g){
	int b,xi;
	double c1[2],c2[2],Cinv;
	int xd3;
	for(b=0;b<g->totalblocks;b++){
		xd3=g->xdim[b][3];
		for(xi=0;xi<xd3;xi++){
			c1[0]=g->covc11[b][xi];
			c1[1]=g->covc12[b][xi];
			c2[0]=g->covc21[b][xi];
			c2[1]=g->covc22[b][xi];
			Cinv=1/cross(c1,c2);
			g->contrac11[b][xi]=c2[1]*Cinv;
			g->contrac12[b][xi]=-c2[0]*Cinv;
			g->contrac21[b][xi]=-c1[1]*Cinv;
			g->contrac22[b][xi]=c1[0]*Cinv;
		}
	}
}
void fillBTerms(Grid*g,BoundaryConditions*bc){
	//bcs: bottom,right,top,left
	//bc: 0:fixed,1:zerogradient
	int b,i,j,vi,hi,xi;
	int st,iv,en,side;
	int vd0,vd1,hd0,hd1,xd0,xd1,xd3;
	char bc0,bc1,bc2,bc3;
	double xi_x, xi_y, eta_x, eta_y;//for non-nodal locations.
	for(b=0;b<g->totalblocks;b++){
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		hd0=g->hdim[b][0];
		hd1=g->hdim[b][1];
		xd0=g->xdim[b][0];
		xd1=g->xdim[b][1];
		xd3=g->xdim[b][3];
		//B11
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				xi_x=g->a22[b][vi]/g->Jvv[b][vi];
				xi_y=-g->a21[b][vi]/g->Jvv[b][vi];
				g->B11[b][vi]=g->Jvv[b][vi]*(pow(xi_x,2)+pow(xi_y,2));
			}
		}
		//B22
		for(j=0;j<hd1;j++){
			for(i=0;i<hd0;i++){
				hi=i+j*hd0;
				eta_x=-g->a12[b][hi]/g->Jhh[b][hi];
				eta_y=g->a11[b][hi]/g->Jhh[b][hi];
				g->B22[b][hi]=g->Jhh[b][hi]*(pow(eta_x,2)+pow(eta_y,2));
			}
		}
		//B12
		for(j=0;j<xd1;j++){
			for(i=0;i<xd0;i++){
				xi=i+j*xd0;
				g->B12[b][xi]=g->Jxx[b][xi]*(g->xixxx[b][xi]*g->etaxxx[b][xi]+
						g->xiyxx[b][xi]*g->etayxx[b][xi]);
			}
		}
		//boundaries zerogradient
		for(side=0;side<4;side++){
			st=bc->side[b][side][0];
			iv=bc->side[b][side][1];
			en=bc->side[b][side][2];
			switch(bc->id[b][side]){
			case 'z':
				for(xi=st+iv;xi<en;xi+=iv){ g->B12[b][xi]=0; }
				if(side%2==0){
					for(xi=st;xi<en;xi+=iv){
						hi=xi-xi/xd0;
						g->B22[b][hi]=0;
					}
				} else {
					for(xi=st;xi<en;xi+=iv){
						vi=xi;
						g->B11[b][vi]=0;
					}
				}
				break;
			}
		}
		//corners
		bc0=bc->id[b][0];
		bc1=bc->id[b][1];
		bc2=bc->id[b][2];
		bc3=bc->id[b][3];
		if((bc3=='z') & (bc0=='z')){
			g->B12[b][0]=0;
		}
		if((bc1=='z') & (bc0=='z')){
			g->B12[b][xd0-1]=0;
		}
		if((bc1=='z') & (bc2=='z')){
			g->B12[b][xd3-1]=0;
		}
		if((bc3=='z') & (bc2=='z')){
			g->B12[b][xd3-xd0]=0;
		}
	}
}
void transformation(Grid* g,BoundaryConditions* bc,Interfaces* is){
	mallocTransformation(g);
	computeCovariant(g);
	fillCovariantcc(g);
	fillCovariantxx(g,bc);
	fillCovariantxxinterfaces(g,is);
	fillCovariantvv(g,bc);
	fillCovariantvvinterfaces(g,is);
	fillCovarianthh(g,bc);
	fillCovarianthhinterfaces(g,is);
	fillJcc(g);
	fillJxx(g);
	fillJvv(g);
	fillJhh(g);
	fillContravariantcc(g);
	fillContravariantxx(g);
	computevareas(g);
	computevareasinterfaces(g,is);
	computehareas(g);
	computehareasinterfaces(g,is);
	fillcovc(g,bc);
	fillcovcinterfaces(g,is);
	fillcontrac(g);

	fillBTerms(g,bc);
}


#endif /* TRANSFORMATION_H_ */
