/*
 * old.h
 *
 *  Created on: Jul 15, 2013
 *      Author: lordvon
 */

#ifndef OLD_H_
#define OLD_H_

/*
 *
 *
		//vertical boundaries
		os=vd0-1;
		for(j=0;j<vd1;j++){
			vi1=j*vd0;
			vi2=j*vd0+os;
			a2[0]=g->a21[b][vi1];
			a2[1]=g->a22[b][vi1];
			if(bc3=='f'){ s->u[b][vi1]=cross(fixedcart,a2); }
			else if(bc3=='w'){ s->u[b][vi1]=0; }
			a2[0]=g->a21[b][vi2];
			a2[1]=g->a22[b][vi2];
			if(bc1=='f'){ s->u[b][vi2]=cross(fixedcart,a2); }
			else if(bc1=='w'){ s->u[b][vi2]=0; }
		}
		//horizontal boundaries
		os=hd3-hd0;
		for(i=0;i<hd0;i++){
			hi1=i;
			hi2=i+os;
			a1[0]=g->a11[b][hi1];
			a1[1]=g->a12[b][hi1];
			if(bc0=='f'){ s->v[b][hi1]=-cross(fixedcart,a1); }
			else if(bc0=='w'){ s->v[b][hi1]=0; }
			a1[0]=g->a11[b][hi2];
			a1[1]=g->a12[b][hi2];
			if(bc2=='f'){ s->v[b][hi2]=-cross(fixedcart,a1); }
			else if(bc2=='w'){ s->v[b][hi2]=0; }
		}
void update(State* s,Grid* g,LinearSystem* ls,BoundaryConditions* bc){
	//Updates u and v from C*s.
	int b,i,j,vi,hi,es,side;
	int vd0,vd1,vd3,hd0,hd1;
	for(b=0;b<g->totalblocks;b++){
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		vd3=g->vdim[b][3];
		hd0=g->hdim[b][0];
		hd1=g->hdim[b][1];
		es=g->estart[b];
		for(j=0;j<vd1;j++){
			for(i=0;i<vd0;i++){
				vi=i+j*vd0;
				side=-1;
				if(i==0){ side=3; }
				else if(i==vd0-1){ side=1; }
				if(side>=0){
					if((bc->id[b][side]=='f') | (bc->id[b][side]=='w')){
						continue;//dont update if fixed or wall bc.
					}
				}
				s->u[b][vi]=ls->cs[es+vi];
			}
		}
		for(j=0;j<hd1;j++){
			for(i=0;i<hd0;i++){
				hi=i+j*hd0;
				side=-1;
				if(j==0){ side=0; }
				else if(j==hd1-1){ side=2; }
				if(side>=0){
					if((bc->id[b][side]=='f') | (bc->id[b][side]=='w')){
						continue;//dont update if fixed or wall bc.
					}
				}
				s->v[b][hi]=ls->cs[es+vd3+hi];
			}
		}
	}
}
void explicitRhs(Grid* g,Momentum* mm,LinearSystem* ls,
		Numerics* n, State* s,Interfaces* is){
	int b,vi,hi,es,ip;
	int hd3,vd3;
	for(b=0;b<g->totalblocks;b++){
		vd3=g->vdim[b][3];
		hd3=g->hdim[b][3];
		es=g->estart[b];
		for(vi=0;vi<vd3;vi++){
			ls->rhs[es+vi]=s->u[b][vi]+n->dt*mm->ut[b][vi];
		}
		for(hi=0;hi<hd3;hi++){
			ls->rhs[es+vd3+hi]=s->v[b][hi]+n->dt*mm->vt[b][hi];
		}
	}
	for(ip=0;ip<is->totalinterfacepoints;ip++){
		ls->rhs[g->totaledges+ip]=0;
	}
}
void explicitEulerInterior(Grid* g,Momentum* mm,BoundaryConditions* bc){
	//Explicit Euler (first-order) time derivatives.
	//Finite-volume formulation.

	//Zero gradient and fixed time derivatives are not filled.
	int b,i,j,vi,hi;
	int vd0,vd1,hd0,hd1,hd3,os;
	int imin,jmin,imax,jmax;
	int vi1,vi2,hi1,hi2;
	char bc0,bc1,bc2,bc3;
	double umom[2],vmom[2],a2[2],a1[2];
	//Cycle through blocks.
	for(b=0;b<g->totalblocks;b++){
		//Get dimensions.
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		hd0=g->hdim[b][0];
		hd1=g->hdim[b][1];
		hd3=g->hdim[b][3];
		bc0=bc->id[b][0];
		bc1=bc->id[b][1];
		bc2=bc->id[b][2];
		bc3=bc->id[b][3];
		//vertical edges
		imin=1;
		imax=vd0-2;
		if(bc->id[b][3]=='i'){ imin=0; }
		if(bc->id[b][1]=='i'){ imax=vd0-1; }
		//
		for(j=0;j<vd1;j++){
			for(i=imin;i<=imax;i++){
				vi=i+j*vd0;
				a2[0]=g->a21[b][vi];
				a2[1]=g->a22[b][vi];
				umom[0]=-mm->I11_1[b][vi]-mm->I12_1[b][vi]+
						mm->I14_1[b][vi]+mm->I15_1[b][vi];
				umom[1]=-mm->I11_2[b][vi]-mm->I12_2[b][vi]+
						mm->I14_2[b][vi]+mm->I15_2[b][vi];
				mm->ut[b][vi]=-cross(a2,umom)/g->vareas[b][vi];
			}
		}
		//vertical boundaries zerogradient
		os=vd0-1;
		for(j=0;j<vd1;j++){
			vi1=j*vd0;
			vi2=j*vd0+os;
			if(bc3=='z'){ mm->ut[b][vi1]=mm->ut[b][vi1+1]; }
			if(bc1=='z'){ mm->ut[b][vi2]=mm->ut[b][vi2-1]; }
		}
		//horizontal edges
		jmin=1;
		jmax=hd1-2;
		if(bc->id[b][0]=='i'){ jmin=0; }
		if(bc->id[b][2]=='i'){ jmax=hd1-1; }
		for(j=jmin;j<=jmax;j++){
			for(i=0;i<hd0;i++){
				hi=i+j*hd0;
				a1[0]=g->a11[b][hi];
				a1[1]=g->a12[b][hi];
				vmom[0]=-mm->I21_1[b][hi]-mm->I22_1[b][hi]+
						mm->I24_1[b][hi]+mm->I25_1[b][hi];
				vmom[1]=-mm->I21_2[b][hi]-mm->I22_2[b][hi]+
						mm->I24_2[b][hi]+mm->I25_2[b][hi];
				mm->vt[b][hi]=cross(a1,vmom)/g->hareas[b][hi];
			}
		}
		//horizontal boundaries zerogradient
		os=hd3-hd0;
		for(i=0;i<hd0;i++){
			hi1=i;
			hi2=i+os;
			if(bc0=='z'){ mm->vt[b][hi1]=mm->vt[b][hi1+hd0]; }
			if(bc2=='z'){ mm->vt[b][hi2]=mm->vt[b][hi2-hd0]; }
		}
	}
}
void fillCartDiffxxInterface(Grid* g,Momentum* mm,Interfaces* is){
	int i,ip;
	int xi1,xi2,vi1,vi2,hi1,hi2;
	int b1,b2;
	for(i=0;i<is->totalinterfaces;i++){
		b1=is->block[i];
		b2=is->adjacentBlock[i];
		switch(is->side[i]){
		case 0:
			for(ip=0;ip<is->xsize[i];ip++){
				xi1=is->blockx[i][ip]-g->xstart[b1];
				xi2=is->adjacentBlockx[i][ip]-g->xstart[b2];
				vi1=xi1;
				vi2=xi2;
				mm->ucarteta[b1][xi1]=mm->ucarteta[b2][xi2]=
						mm->ucartvv[b1][vi1]-mm->ucartvv[b2][vi2-g->vdim[b2][0]];
				mm->vcarteta[b1][xi1]=mm->vcarteta[b2][xi2]=
						mm->vcartvv[b1][vi1]-mm->vcartvv[b2][vi2-g->vdim[b2][0]];
			}
			break;
		case 1:
			for(ip=0;ip<is->xsize[i];ip++){
				xi1=is->blockx[i][ip]-g->xstart[b1];
				xi2=is->adjacentBlockx[i][ip]-g->xstart[b2];
				hi1=xi1-xi1/g->xdim[b1][0];
				hi2=xi2-xi2/g->xdim[b2][0];
				mm->ucartxi[b1][xi1]=mm->ucartxi[b2][xi2]=
						mm->ucarthh[b2][hi2]-mm->ucarthh[b1][hi1-1];
				mm->vcartxi[b1][xi1]=mm->vcartxi[b2][xi2]=
						mm->vcarthh[b2][hi2]-mm->vcarthh[b1][hi1-1];
			}
			break;
		case 2:
			for(ip=0;ip<is->xsize[i];ip++){
				xi1=is->blockx[i][ip]-g->xstart[b1];
				xi2=is->adjacentBlockx[i][ip]-g->xstart[b2];
				vi1=xi1;
				vi2=xi2;
				mm->ucarteta[b1][xi1]=mm->ucarteta[b2][xi2]=
						mm->ucartvv[b2][vi2]-mm->ucartvv[b1][vi1-g->vdim[b1][0]];
				mm->vcarteta[b1][xi1]=mm->vcarteta[b2][xi2]=
						mm->vcartvv[b2][vi2]-mm->vcartvv[b1][vi1-g->vdim[b1][0]];
			}
			break;
		case 3:
			for(ip=0;ip<is->xsize[i];ip++){
				xi1=is->blockx[i][ip]-g->xstart[b1];
				xi2=is->adjacentBlockx[i][ip]-g->xstart[b2];
				hi1=xi1-xi1/g->xdim[b1][0];
				hi2=xi2-xi2/g->xdim[b2][0];
				//printf("uhh: %g %g\n",mm->ucarthh[b1][hi1],mm->ucarthh[b2][hi2-1]);
				mm->ucartxi[b1][xi1]=mm->ucartxi[b2][xi2]=
						mm->ucarthh[b1][hi1]-mm->ucarthh[b2][hi2-1];
				mm->vcartxi[b1][xi1]=mm->vcartxi[b2][xi2]=
						mm->vcarthh[b1][hi1]-mm->vcarthh[b2][hi2-1];
			}
			break;
		default:
			printf("Error! Side not valid in 'fillcartdiffinterface()'");
			break;
		}
	}
}
void enforcecornercartdiffxx(Grid* g,Momentum* mm,Interfaces* is){
	//to be called just before fillcartderivativexx is called.
	int i,xs;
	int xc1_1,xc2_1,xc1_2,xc2_2;
	int b1,b2;
	for(i=0;i<is->totalinterfaces;i++){
		b1=is->block[i];
		b2=is->adjacentBlock[i];
		xs=is->xsize[i];
		xc1_1=is->blockx[i][0]-g->xstart[b1];
		xc2_1=is->adjacentBlockx[i][0]-g->xstart[b2];
		xc1_2=is->blockx[i][xs-1]-g->xstart[b1];
		xc2_2=is->adjacentBlockx[i][xs-1]-g->xstart[b2];
		mm->ucarteta[b2][xc2_1]=mm->ucarteta[b1][xc1_1];
		mm->vcarteta[b2][xc2_1]=mm->vcarteta[b1][xc1_1];
		mm->ucartxi[b2][xc2_1]=mm->ucartxi[b1][xc1_1];
		mm->vcartxi[b2][xc2_1]=mm->vcartxi[b1][xc1_1];
		mm->ucarteta[b2][xc2_2]=mm->ucarteta[b1][xc1_2];
		mm->vcarteta[b2][xc2_2]=mm->vcarteta[b1][xc1_2];
		mm->ucartxi[b2][xc2_2]=mm->ucartxi[b1][xc1_2];
		mm->vcartxi[b2][xc2_2]=mm->vcartxi[b1][xc1_2];
	}
}
void enforcecartxxcornerssub(Grid* g,Momentum* mm,BoundaryConditions* bc,Interfaces* is,
		int b,int s1,int s2,int bc1,int bc2,int xc,
		double *cart){
	int wall=0;
	if((bc1=='w') | (bc2=='w')){ wall=1; }
	if(bc1=='i'){
		if(bc->id[is->conn[b][s1]][s2]=='w'){ wall=1; }
	}
	if(bc2=='i'){
		if(bc->id[is->conn[b][s2]][s1]=='w'){ wall=1; }
	}
	if(wall) {
		mm->ucartxx[b][xc]=0;
		mm->vcartxx[b][xc]=0;
	} else if((bc1=='f') | (bc2=='f')) {
		mm->ucartxx[b][xc]=cart[0];
		mm->vcartxx[b][xc]=cart[1];
	}
}
void enforcecartxxcorners(State* s,Momentum* mm,Grid* g,BoundaryConditions* bc,
		Interfaces* is){
	//to be called after the xxinterfaces function
	int xd0,xd3;
	int b;
	int xc0,xc1,xc2,xc3;
	char bc0,bc1,bc2,bc3;
	double cart[2];
	for(b=0;b<g->totalblocks;b++){
		xd0=g->xdim[b][0];
		xd3=g->xdim[b][3];
		xc0=0;
		xc1=xd0-1;
		xc2=xd3-1;
		xc3=xd3-xd0;
		bc0=bc->id[b][0];
		bc1=bc->id[b][1];
		bc2=bc->id[b][2];
		bc3=bc->id[b][3];
		if((bc0=='f') | (bc1=='f') | (bc2=='f') | (bc3=='f')){
			cart[0]=bc->f[b][0];
			cart[1]=bc->f[b][1];
		}
		enforcecartxxcornerssub(g,mm,bc,is,b,3,0,bc3,bc0,xc0,cart);
		enforcecartxxcornerssub(g,mm,bc,is,b,1,0,bc1,bc0,xc1,cart);
		enforcecartxxcornerssub(g,mm,bc,is,b,1,2,bc1,bc2,xc2,cart);
		enforcecartxxcornerssub(g,mm,bc,is,b,3,2,bc3,bc2,xc3,cart);
	}
}
void boundarycartxx(Grid* g,Momentum* mm,BoundaryConditions* bc){
	int b,xi;
	int in,iv,st,en,si;
	for(b=0;b<g->totalblocks;b++){
		for(si=0;si<4;si++){
			st=bc->side[b][si][0];
			en=bc->side[b][si][2];
			iv=bc->side[b][si][1];
			in=bc->side[b][si][3];
			//printf("%d,%d,%d,%d\n",st,en,iv,in);
			switch(bc->id[b][si]){
			case 'i':
				for(xi=st;xi<=en;xi+=iv){
					mm->ucartxx[b][xi]=(g->a11xx[b][xi]*mm->uxx[b][xi]+g->a21xx[b][xi]*mm->vxx[b][xi])/g->Jxx[b][xi];
					mm->vcartxx[b][xi]=(g->a12xx[b][xi]*mm->uxx[b][xi]+g->a22xx[b][xi]*mm->vxx[b][xi])/g->Jxx[b][xi];
				}
				break;
			case 'w':
				for(xi=st;xi<=en;xi+=iv){
					mm->ucartxx[b][xi]=0;
					mm->vcartxx[b][xi]=0;
				}
				break;
			case 'f':
				for(xi=st;xi<=en;xi+=iv){
					mm->ucartxx[b][xi]=bc->f[b][0];
					mm->vcartxx[b][xi]=bc->f[b][1];
				}
				break;
			case 'z'://assume boundary-grid orthogonality.
				for(xi=st;xi<=en;xi+=iv){
					mm->ucartxx[b][xi]=mm->ucartxx[b][xi+in];
					mm->vcartxx[b][xi]=mm->vcartxx[b][xi+in];
				}
				break;
			}
		}
	}
}
void enforcexxcornerssub(Grid* g,Momentum* mm,BoundaryConditions* bc,Interfaces* is,
		int b,int s1,int s2,int bc1,int bc2,int xc,
		double *cart){
	int wall=0;
	double a[2];
	if((bc1=='w') | (bc2=='w')){ wall=1; }
	if(bc1=='i'){
		if(bc->id[is->conn[b][s1]][s2]=='w'){ wall=1; }
	}
	if(bc2=='i'){
		if(bc->id[is->conn[b][s2]][s1]=='w'){ wall=1; }
	}
	if(wall) {
		mm->uxx[b][xc]=0;
		mm->vxx[b][xc]=0;
	} else if((bc1=='f') | (bc2=='f')) {
		a[0]=g->a21xx[b][xc];
		a[1]=g->a22xx[b][xc];
		mm->uxx[b][xc]=cross(cart,a);
		a[0]=g->a11xx[b][xc];
		a[1]=g->a12xx[b][xc];
		mm->vxx[b][xc]=-cross(cart,a);
	}
}
void enforcexxcorners(State* s,Momentum* mm,Grid* g,BoundaryConditions* bc,
		Interfaces* is){
	//to be called after the xxinterfaces function
	int xd0,xd3;
	int b;
	int xc0,xc1,xc2,xc3;
	char bc0,bc1,bc2,bc3;
	double cart[2];
	for(b=0;b<g->totalblocks;b++){
		xd0=g->xdim[b][0];
		xd3=g->xdim[b][3];
		xc0=0;
		xc1=xd0-1;
		xc2=xd3-1;
		xc3=xd3-xd0;
		bc0=bc->id[b][0];
		bc1=bc->id[b][1];
		bc2=bc->id[b][2];
		bc3=bc->id[b][3];
		if((bc0=='f') | (bc1=='f') | (bc2=='f') | (bc3=='f')){
			cart[0]=bc->f[b][0];
			cart[1]=bc->f[b][1];
		}
		enforcexxcornerssub(g,mm,bc,is,b,3,0,bc3,bc0,xc0,cart);
		enforcexxcornerssub(g,mm,bc,is,b,1,0,bc1,bc0,xc1,cart);
		enforcexxcornerssub(g,mm,bc,is,b,1,2,bc1,bc2,xc2,cart);
		enforcexxcornerssub(g,mm,bc,is,b,3,2,bc3,bc2,xc3,cart);
	}
}
void symlapBoundaryCorners2(SymLap*s,Grid*g,Interfaces*is){
	//fixed/wall boundary corner enforcement in blocks that have walls at corners (but no wall sides, rather interfaces)
	int ic,c,b,couch,oppositeb,oppositec;
	int mc,mr,entryi;
	int cornerci,cornerxi,oppositecornerci;
	double sign;
	printf("symlap diagonal adjacency corners: %d\n",countSymlapBoundaryCorners(is));

	s->csr2.nn=countSymlapBoundaryCorners(is);
	s->csr2.ri=malloc(sizeof(int)*s->csr2.nn);
	s->csr2.v=malloc(sizeof(double)*s->csr2.nn);
	s->csr2.ci=malloc(sizeof(int)*s->csr2.nn);
	int *bchecklist=malloc(sizeof(int)*g->totalblocks);
	for(b=0;b<g->totalblocks;b++){ bchecklist[b]=0; }
	int **cornerblocks=malloc(sizeof(int *)*s->csr2.nn);
	for(c=0;c<s->csr2.nn;c++){ cornerblocks[c]=malloc(sizeof(int)*2); }
	int tic=is->totalicorners;
	//cycle through interface intersections.
	entryi=0;
	for(ic=0;ic<tic;ic++){
		//determine couched cells.
		for(c=0;c<4;c++){
			b=is->icorners[ic][c];
			if(b>=0){
				oppositec=(c+2)%4;
				oppositeb=is->icorners[ic][oppositec];
				couch=checkCouch(ic,c,is);
				cornerci=getcornercellci(g,b,c);
				cornerxi=getcornerxi(g,b,c);
				if((couch>=0) & (oppositeb<0)){//case that wall is at corner.
					//fprintf(is->out,"Corner procedure reached.\n");
					if(c%2==0){ sign=-1; } else { sign=1; }
					s->d0[b][cornerci]-=0.5*sign*
							g->B12[b][cornerxi];
				} else if (oppositeb>=0){//case that another non-adjacent block is at corner.
					if((bchecklist[oppositeb]==0) & (bchecklist[b]==0)){
						oppositecornerci=getcornercellci(g,oppositeb,oppositec);
						mc=g->cstart[b]+cornerci;
						mr=g->cstart[oppositeb]+oppositecornerci;
						s->csr2.ri[entryi]=mr;
						s->csr2.ci[entryi]=mc;
						if(c%2==0){ sign=1; } else { sign=-1; }
						s->csr2.v[entryi]=0.5*sign*g->B12[b][cornerxi];
						bchecklist[oppositeb]=1;
						bchecklist[b]=1;
						entryi++;
						printf("Entry %d\n",entryi);
					}
				}
			}
		}
	}
	free(bchecklist);
	free(cornerblocks);
}
int counticornerssub2(Interfaces* is,int i,int cornercount,int end,
		int***ifromb,int**ichecklist){
	//end: 0 or 1. Corresponds to each end of the interface line.
	int b1,b2,adjacentSide,s1,i1,i2,icheck1,icheck2,is1,is2;
	b1=is->block[i];
	b2=is->adjacentBlock[i];
	s1=is->side[i];
	adjacentSide=(s1+1+2*end)%4;//adjacent side to interface (same number for both blocks)
	i1=ifromb[b1][adjacentSide][0];
	i2=ifromb[b2][adjacentSide][0];
	//printf("i: %d, i1: %d, i2: %d\n",i,i1,i2);
	//printf("ichecklist: %d, %d, %d, %d\n",ichecklist[0],ichecklist[1],ichecklist[2],ichecklist[3]);
	//check if the adjacent sides are interfaces.
	is1=i1 >= 0;
	is2=i2 >= 0;
	//check to see if they have been checked off.
	icheck1=icheck2=1;
	if(is1){ icheck1=ichecklist[i1][end]<1; }
	if(is2){ icheck2=ichecklist[i2][end]<1; }
	if((is1 | is2) & icheck1 & icheck2){//at least one is an interface, and both are not checked off.
		if(is1){ ichecklist[i1][end]=1; }
		if(is2){ ichecklist[i2][end]=1; }
		cornercount++;
		printf("Counted here %d %d %d\n",b1,b2,end);
	}
	return cornercount;
}
int counticorners2(Interfaces*is,int***ifromb){
	int i;
	int cornercount=0;
	int ti=is->totalinterfaces;
	int**ichecklist=malloc(sizeof(int)*ti);//boolean
	for(i=0;i<ti;i++){ ichecklist[i]=malloc(sizeof(int)*2);ichecklist[i][0]=0;ichecklist[i][1]=0; }//initialize check list.
	for(i=0;i<ti;i++){
		//end 1
		cornercount=counticornerssub2(is,i,cornercount,1,ifromb,ichecklist);
		//end 0
		cornercount=counticornerssub2(is,i,cornercount,0,ifromb,ichecklist);
	}
	return cornercount;
}
void fillnutxx(State* s,Grid* g,Momentum* mm){
	int b,i,j,xi,ci,cd0;
	for(b=0;b<g->totalblocks;b++){
		for(j=0;j<g->xdim[b][1]-1;j++){
			for(i=0;i<g->xdim[b][0]-1;i++){
				xi=i+j*g->xdim[b][0];
				ci=xi-j;
				cd0=g->cdim[b][0];
				mm->nutxx[b][xi]=0.25*(
						s->nut[b][ci]+
						s->nut[b][ci-1]+
						s->nut[b][ci-cd0]+
						s->nut[b][ci-cd0-1]);
			}
		}
	}
}
void fillnutxxboundary(State* s,Momentum* mm,Grid* g,Interfaces* is){
	int i,ip;
	int b1,b2,xi1,xi2,vi1,vi2,hi1,hi2;
	double avgv,avgu;
	for(i=0;i<is->totalinterfaces;i++){
		b1=is->block[i];
		b2=is->adjacentBlock[i];
		for(ip=1;ip<is->xsize[i]-1;ip++){
			xi1=is->blockx[i][ip]-g->xstart[b1];
			xi2=is->adjacentBlockx[i][ip]-g->xstart[b2];
			switch(is->side[i]){
			case 0:
				vi1=xi1;
				vi2=xi2-g->xdim[b2][0];
				break;
			case 2:
				vi1=xi1-g->xdim[b1][0];
				vi2=xi2;
				break;
			case 1:
				hi1=xi1-xi1/g->xdim[b1][0]-1;
				hi2=xi2-xi2/g->xdim[b2][0];
				break;
			case 3:
				hi1=xi1-xi1/g->xdim[b1][0];
				hi2=xi2-xi2/g->xdim[b2][0]-1;
				break;
			}
			if(is->side[i]%2!=0){
				avgv=0.5*(s->v[b1][hi1]+s->v[b2][hi2]);
				mm->vxx[b1][xi1]=avgv;
				mm->vxx[b2][xi2]=avgv;
			} else {
				avgu=0.5*(s->u[b1][vi1]+s->u[b2][vi2]);
				mm->uxx[b1][xi1]=avgu;
				mm->uxx[b2][xi2]=avgu;
			}
		}
	}
}
double interpnut(State* s,Momentum* mm,int xd0,int xd1,int cd0,int b,int i, int j){
	//Does not account for corner case.
	int xi,ci;
	double returnnu;
	xi=i+j*xd0;
	ci=xi-j;
	if(i==0){
		returnnu=0.5*(s->nut[b][ci]+s->nut[b][ci-cd0]);
	} else if(i==xd0-1){
		returnnu=0.5*(s->nut[b][ci-1]+s->nut[b][ci-1-cd0]);
	} else if(j==0){
		returnnu=0.5*(s->nut[b][ci]+s->nut[b][ci-1]);
	} else if(j==xd1-1){
		returnnu=0.5*(s->nut[b][ci-cd0]+s->nut[b][ci-cd0-1]);
	} else {
		returnnu=0.25*(s->nut[b][ci]+s->nut[b][ci-1]+s->nut[b][ci-cd0]+s->nut[b][ci-cd0-1]);
	}
	return returnnu;
}
	//
	s->iloc=malloc(sizeof(int**)*ti);
	s->valcol=malloc(sizeof(double**)*ti);
	for(i=0;i<ti;i++){
		cs=is->xsize[i]-1;
		s->iloc[i]=malloc(sizeof(int*)*cs);
		s->valcol[i]=malloc(sizeof(double*)*cs);
		for(j=0;j<cs;j++){
			s->iloc[i][j]=malloc(sizeof(int)*3);
			s->valcol[i][j]=malloc(sizeof(double)*3);
		}
	}
void symlapInterface(SymLap*s,Grid*g,Interfaces* is){
	int i,b1,b2,side,ip;
	int xi1,xi2,vi1,vi2,ci1,ci2,hi1,hi2;
	int civ2,xiv,adj1,adj2;
	double bbb;
	int i1,i2;
	//int tc,civ1,rowi,coli,slin;
	//tc=g->totalcells;
	//interfaces, currently ignores corners of multiple interfaces.
	for(i=0;i<is->totalinterfaces;i++){
		b1=is->block[i];
		side=is->side[i];
		b2=is->adjacentBlock[i];
		for(ip=0;ip<is->xsize[i]-1;ip++){
			xi1=vi1=is->blockx[i][ip]-g->xstart[b1];
			xi2=vi2=is->adjacentBlockx[i][ip]-g->xstart[b2];
			ci1=hi1=xi1-xi1/g->xdim[b1][0];
			ci2=hi2=xi2-xi2/g->xdim[b2][0];
			if(side%2==0){
				xiv=1;
				//civ1=civ2=1;
				bbb=g->B22[b1][hi1];
			} else {
				xiv=g->xdim[b1][0];
				//civ1=g->cdim[b1][0];
				civ2=g->cdim[b2][0];
				bbb=g->B11[b1][vi1];
			}
			if(side==0){
				adj1=0;
				adj2=-g->cdim[b2][0];
			} else if(side==2){
				adj1=-g->cdim[b1][0];
				adj2=0;
			} else if(side==1){
				adj1=-1;
				adj2=0;
			} else if(side==3){
				adj1=0;
				adj2=-1;
			}
			//column method.
			i1=g->cstart[b1]+ci1+adj1;
			i2=g->cstart[b2]+ci2+adj2;
			s->iloc[i][ip][0]=i1;
			s->iloc[i][ip][1]=i2;
			if((ip>0) & (ip<is->xsize[i]-2)){
				s->iloc[i][ip][2]=civ2;
				s->valcol[i][ip][0]=g->B12[b1][xi1];
				s->valcol[i][ip][2]=g->B12[b1][xi1+xiv];
				//s->A[i1+(i2-civ2)*tc]=g->B12[b1][xi1];
				//s->A[i1+(i2+civ2)*tc]=g->B12[b1][xi1+xiv];
				//s->A[i1*tc+i2-civ2]=g->B12[b1][xi1];
				//s->A[i1*tc+i2+civ2]=g->B12[b1][xi1+xiv];
			} else {
				s->iloc[i][ip][2]=0;
			}
			s->valcol[i][ip][1]=bbb;
			//s->A[i1+i2*tc]=bbb;
			//s->A[i1*tc+i2]=bbb;
		}
	}
}
void symlapInterfaceMult(SymLap*s,Interfaces*is){
	//adds on to s->b
	int i,c,i1,i2,stride;
	for(i=0;i<is->totalinterfaces;is++){
		for(c=0;c<is->xsize[i]-1;c++){
			i1=s->iloc[i][c][0];
			i2=s->iloc[i][c][1];
			stride=s->iloc[i][c][2];
			if(stride>0){
				s->b[i2-stride]	+=s->x[i1]			*s->valcol[i][c][0];
				s->b[i1]		+=s->x[i2-stride]	*s->valcol[i][c][0];
				s->b[i2+stride]	+=s->x[i1]			*s->valcol[i][c][2];
				s->b[i1]		+=s->x[i2+stride]	*s->valcol[i][c][2];
			}
			s->b[i2]+=s->x[i1]*s->valcol[i][c][1];
			s->b[i1]+=s->x[i2]*s->valcol[i][c][1];
		}
	}
}
void slmult(SymLap*s,Grid*g,Interfaces*is){
	int b,i,cd0,cd3,cs,offset;
	int row,col;
	for(b=0;b<g->totalblocks;b++){
		cs=g->cstart[b];
		cd0=g->cdim[b][0];
		cd3=g->cdim[b][3];
		offset=0;
		for(i=0;i<cd3-offset;i++){
			row=cs+i;col=row+offset;
			s->b[row]+=s->x[col]*s->d0[b][i];
			s->b[col]+=s->x[row]*s->d0[b][i];
		}
		offset=1;
		for(i=0;i<cd3-offset;i++){
			row=cs+i;col=row+offset;
			s->b[row]+=s->x[col]*s->d1[b][i];
			s->b[col]+=s->x[row]*s->d1[b][i];
		}
		offset=cd0-1;
		for(i=0;i<cd3-offset;i++){
			row=cs+i;col=row+offset;
			s->b[row]+=s->x[col]*s->d2[b][i];
			s->b[col]+=s->x[row]*s->d2[b][i];
		}
		offset=cd0;
		for(i=0;i<cd3-offset;i++){
			row=cs+i;col=row+offset;
			s->b[row]+=s->x[col]*s->d3[b][i];
			s->b[col]+=s->x[row]*s->d3[b][i];
		}
		offset=cd0+1;
		for(i=0;i<cd3-offset;i++){
			row=cs+i;col=row+offset;
			s->b[row]+=s->x[col]*s->d4[b][i];
			s->b[col]+=s->x[row]*s->d4[b][i];
		}
	}
	symlapInterfaceMult(s,is);
}
void fillcartderivativecc(Grid* g,Momentum* mm,BoundaryConditions* bc){
	double ucarteta,vcarteta,ucartxi,vcartxi;
	double Jinv;
	int b,ci,vi,hi;
	int cd3,cd0,hd0;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		cd0=g->cdim[b][0];
		hd0=g->hdim[b][0];
		for(ci=0;ci<cd3;ci++){
			vi=ci+ci/cd0;
			hi=ci;
			ucarteta=mm->ucarthh[b][hi+hd0]-mm->ucarthh[b][hi];
			vcarteta=mm->vcarthh[b][hi+hd0]-mm->vcarthh[b][hi];
			ucartxi=mm->ucartvv[b][vi+1]-mm->ucartvv[b][vi];
			vcartxi=mm->vcartvv[b][vi+1]-mm->vcartvv[b][vi];
			Jinv=1/g->Jcc[b][ci];
			mm->uxcc[b][ci]=Jinv*(g->a22cc[b][ci]*ucartxi-g->a12cc[b][ci]*ucarteta);
			mm->vxcc[b][ci]=Jinv*(g->a22cc[b][ci]*vcartxi-g->a12cc[b][ci]*vcarteta);
			mm->uycc[b][ci]=Jinv*(-g->a21cc[b][ci]*ucartxi+g->a11cc[b][ci]*ucarteta);
			mm->vycc[b][ci]=Jinv*(-g->a21cc[b][ci]*vcartxi+g->a11cc[b][ci]*vcarteta);
		}
	}
}
void computevareas2(Grid* g){
	int i,j;
	int vi,xi;
	int b,vd0,vd1,xd0;
	double x1[2],x2[2],x3[2],x4[2];
	double x2p[2],x3p[2];
	double a1[2],b1[2],a2[2],b2[2];
	double area1,area2;
	for(b=0;b<g->totalblocks;b++){
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		xd0=g->xdim[b][0];
		for (j=0;j<vd1;j++){
			for (i=1;i<vd0-1;i++){
				vi=i+j*vd0;
				xi=vi;
				x1[0]=g->x[b][xi-1];
				x1[1]=g->y[b][xi-1];
				x2[0]=g->x[b][xi];
				x2[1]=g->y[b][xi];
				x3[0]=g->x[b][xi+xd0];
				x3[1]=g->y[b][xi+xd0];
				x4[0]=g->x[b][xi-1+xd0];
				x4[1]=g->y[b][xi-1+xd0];
				x2p[0]=g->x[b][xi+1];
				x2p[1]=g->y[b][xi+1];
				x3p[0]=g->x[b][xi+1+xd0];
				x3p[1]=g->y[b][xi+1+xd0];
				a1[0]=x3[0]-(x1[0]+x2[0])/2;
				a1[1]=x3[1]-(x1[1]+x2[1])/2;
				b1[0]=(x4[0]+x3[0])/2-x2[0];
				b1[1]=(x4[1]+x3[1])/2-x2[1];
				a2[0]=(x3p[0]+x3[0])/2-x2[0];
				a2[1]=(x3p[1]+x3[1])/2-x2[1];
				b2[0]=x3[0]-(x2p[0]+x2[0])/2;
				b2[1]=x3[1]-(x2p[1]+x2[1])/2;
				area1=0.5*cross(a1,b1);
				area2=0.5*cross(a2,b2);
				if((area1<0) | (area2<0)){
					printf("Negative varea detected!");
				}
				g->vareas[b][vi]=area1+area2;
			}
		}
	}
}
void computehareas2(Grid* g){
	int b,i,j;
	int hi,xi;
	int hd1,hd0,xd0;
	double x1[2],x2[2],x3[2],x4[2];
	double x3p[2],x4p[2];
	double a1[2],b1[2],a2[2],b2[2];
	double area1,area2;
	for(b=0;b<g->totalblocks;b++){
		hd0=g->hdim[b][0];
		hd1=g->hdim[b][1];
		xd0=g->xdim[b][0];
		for (j=1;j<hd1-1;j++){
			for (i=0;i<hd0;i++){
				hi=i+j*hd0;
				xi=hi+j;
				x1[0]=g->x[b][xi-xd0];
				x1[1]=g->y[b][xi-xd0];
				x2[0]=g->x[b][xi+1-xd0];
				x2[1]=g->y[b][xi+1-xd0];
				x3[0]=g->x[b][xi+1];
				x3[1]=g->y[b][xi+1];
				x4[0]=g->x[b][xi];
				x4[1]=g->y[b][xi];
				x3p[0]=g->x[b][xi+1+xd0];
				x3p[1]=g->y[b][xi+1+xd0];
				x4p[0]=g->x[b][xi+xd0];
				x4p[1]=g->y[b][xi+xd0];
				a1[0]=x3[0]-(x1[0]+x4[0])/2;
				a1[1]=x3[1]-(x1[1]+x4[1])/2;
				b1[0]=x4[0]-(x2[0]+x3[0])/2;
				b1[1]=x4[1]-(x2[1]+x3[1])/2;
				a2[0]=(x3p[0]+x3[0])/2-x4[0];
				a2[1]=(x3p[1]+x3[1])/2-x4[1];
				b2[0]=(x4p[0]+x4[0])/2-x3[0];
				b2[1]=(x4p[1]+x4[1])/2-x3[1];
				area1=0.5*cross(a1,b1);
				area2=0.5*cross(a2,b2);
				if((area1<0) | (area2<0)){
					printf("Negative harea detected!");
				}
				g->hareas[b][hi]=area1+area2;
			}
		}
	}
}
void fillcontrac2(Grid* g,BoundaryConditions* bc){
	int b,i,j,xi,i1,i2;
	double c1[2],c2[2],Cinv,factor;
	int xd1,xd0;
	for(b=0;b<g->totalblocks;b++){
		xd0=g->xdim[b][0];
		xd1=g->xdim[b][1];
		for(j=0;j<xd1;j++){
			for(i=0;i<xd0;i++){
				xi=i+j*xd0;
				if(i==0){
					if(bc->id[b][3]!='i'){
						i1=xi;
						i2=xi+1;
						factor=1;
					}
				} else if(i==xd0-1){
					if(bc->id[b][1]!='i'){
						i1=xi-1;
						i2=xi;
						factor=1;
					}
				} else {
					i1=xi-1;
					i2=xi+1;
					factor=0.5;
				}
				c1[0]=factor*(g->x[b][i2]-g->x[b][i1]);
				c1[1]=factor*(g->y[b][i2]-g->y[b][i1]);
				if(j==0){
					if(bc->id[b][0]!='i'){
						i1=xi;
						i2=xi+xd0;
						factor=1;
					}
				} else if(j==xd1-1){
					if(bc->id[b][2]!='i'){
						i1=xi-xd0;
						i2=xi;
						factor=1;
					}
				} else {
					i1=xi-xd0;
					i2=xi+xd0;
					factor=0.5;
				}
				c2[0]=factor*(g->x[b][i2]-g->x[b][i1]);
				c2[1]=factor*(g->y[b][i2]-g->y[b][i1]);
				Cinv=1/cross(c1,c2);
				g->contrac11[b][xi]=c2[1]*Cinv;
				g->contrac12[b][xi]=-c2[0]*Cinv;
				g->contrac21[b][xi]=-c1[1]*Cinv;
				g->contrac22[b][xi]=c1[0]*Cinv;
			}
		}
	}
}
double getcartxi(int xd0,double *carthh, double *cartxx, int i, int j){
	int xi=i+j*xd0, hi=xi-j;
	double cart1,cart2,factor;
	if(i==0){
		cart1=cartxx[xi+1];factor=1;
		//cart1=carthh[hi];factor=2;
		cart2=cartxx[xi];
	} else if(i==xd0-1){
		cart1=cartxx[xi];
		//cart2=carthh[hi-1];factor=2;
		cart2=cartxx[xi-1];factor=1;
	} else {
		cart1=carthh[hi];
		cart2=carthh[hi-1];
		factor=1;
	}
	return factor*(cart1-cart2);
}
double getcarteta(int xd0,int xd1,int vd0,
		double *cartvv, double *cartxx,
		int i, int j){
	int xi=i+j*xd0, vi=xi;
	double cart1,cart2,factor;
	if(j==0){
		cart1=cartxx[xi+xd0];factor=1;
		//cart1=cartvv[vi];factor=2;
		cart2=cartxx[xi];
	} else if(j==xd1-1){
		cart1=cartxx[xi];
		cart2=cartxx[xi-xd0];factor=1;
		//cart2=cartvv[vi-vdim[0]];factor=2;
	} else {
		cart1=cartvv[vi];
		cart2=cartvv[vi-vd0];
		factor=1;
	}
	return factor*(cart1-cart2);
}
void fillcartderivativexx(Grid* g,Momentum* mm,BoundaryConditions* bc){
	//ucarthh,vcarthh,ucartvv,vcartvv,ucartxx,vcartxx
	//
	//uxxx,vxxx,uyxx,vyxx
	int b,i,j,xi;
	double ucarteta,vcarteta,ucartxi,vcartxi;
	int xd1,xd0,vd0;
	for(b=0;b<g->totalblocks;b++){
		xd1=g->xdim[b][1];
		xd0=g->xdim[b][0];
		vd0=g->vdim[b][0];
		for(j=0;j<xd1;j++){
			for(i=0;i<xd0;i++){
				xi=i+j*xd0;
				ucartxi=getcartxi(xd0,mm->ucarthh[b],mm->ucartxx[b],i,j);
				ucarteta=getcarteta(xd0,xd1,vd0,mm->ucartvv[b],mm->ucartxx[b],i,j);
				vcartxi=getcartxi(xd0,mm->vcarthh[b],mm->vcartxx[b],i,j);
				vcarteta=getcarteta(xd0,xd1,vd0,mm->vcartvv[b],mm->vcartxx[b],i,j);
				mm->uxxx[b][xi]=ucartxi*g->contrac11[b][xi]+ucarteta*g->contrac21[b][xi];
				mm->uyxx[b][xi]=ucartxi*g->contrac12[b][xi]+ucarteta*g->contrac22[b][xi];
				mm->vxxx[b][xi]=vcartxi*g->contrac11[b][xi]+vcarteta*g->contrac21[b][xi];
				mm->vyxx[b][xi]=vcartxi*g->contrac12[b][xi]+vcarteta*g->contrac22[b][xi];
			}
		}
	}
}
void readBlockDict(Grid* grid,char *inputfile){
	int maxchars=200;
	int blocknumber,totalblocks,interfacecount,i;
	char bc0,bc1,bc2,bc3;
	int v0,v1,v2,v3;
	char line[maxchars];
	char filename[maxchars];
	FILE * file = fopen (inputfile, "rt");
	int success;
	//Read total number of blocks.
	while(fgets(line, maxchars, file) != NULL){
		success = sscanf(line, "%d", &totalblocks);
		if(success>0){ break; }
	}
	grid->totalblocks=totalblocks;
	grid->bc=charmalloc(grid->totalblocks,4);
	grid->bcval=mbinfo(grid->totalblocks);
	grid->blockNames=charmalloc(grid->totalblocks,200);
	//Read file names for each block.
	interfacecount=0;
	while(fgets(line, maxchars, file) != NULL){
		success = sscanf (line, "%d %s %c=%d %c=%d %c=%d %c=%d",&blocknumber,filename,
				&bc0,&v0,&bc1,&v1,&bc2,&v2,&bc3,&v3);
		if(success>0){
			//printf("%s\n",filename);
			strcpy(grid->blockNames[blocknumber],filename);
					//strcat(blockdirectory,strcat(filename,blockextension)));
			grid->bc[blocknumber][0]=bc0;
			grid->bc[blocknumber][1]=bc1;
			grid->bc[blocknumber][2]=bc2;
			grid->bc[blocknumber][3]=bc3;
			grid->bcval[blocknumber][0]=v0;
			grid->bcval[blocknumber][1]=v1;
			grid->bcval[blocknumber][2]=v2;
			grid->bcval[blocknumber][3]=v3;
			for(i=0;i<4;i++){
				if(grid->bc[blocknumber][i]=='i'){ interfacecount++; }
			}
		}
	}
	if((interfacecount%2)!=0){
		printf("Error: Individual interface declarations not even!\n");
	}
	grid->totalinterfaces=interfacecount/2;
}
int determineSide(int b1,int b2,double matchtol,
		int **xdim,double **gridX,double **gridY){
	int side=-1;
	int corner0,corner1,corner2,corner3;
	corner0=blockMatch(0,						xdim[b2][3]-xdim[b2][0],b1,b2,gridX,gridY,matchtol);
	if(corner0>0){ side=0; }
	corner1=blockMatch(xdim[b1][0]-1,			0,b1,b2,gridX,gridY,matchtol);
	if(corner1>0){ side=1; }
	corner2=blockMatch(xdim[b1][3]-xdim[b1][0],	0,b1,b2,gridX,gridY,matchtol);
	if(corner2>0){ side=2; }
	corner3=blockMatch(0,						xdim[b2][0]-1,b1,b2,gridX,gridY,matchtol);
	if(corner3>0){ side=3; }
	if(side==-1){
		printf("Interface not matched (between blocks %d and %d)!\n",b1,b2);
	}
	return side;
}
int getTotalblocks(char *inputfile){
	int maxchars=200,totalblocks;
	char line[maxchars];
	FILE * file = fopen (inputfile, "rt");
	while(fgets(line, maxchars, file) != NULL){
		if(sscanf(line, "%d", &totalblocks)>0){ break; }
	}
	fclose(file);
	return totalblocks;
}
void readSingleBlockGrid(char *name,int *xdim,double *gridX,double *gridY){
	int maxchars=200;
	char line[maxchars];
	double xcoord=0,ycoord=0;
	int success;

	FILE * file = fopen (name, "rt");
	//Get grid dimensions.
	int dim[3];
	while(fgets(line, maxchars, file) != NULL)
	{
		success = sscanf (line, "%d %d %d", &dim[0], &dim[1], &dim[2]);
		if(success>0){
			printf("%s dimensions: %d %d %d.\n",name,dim[0],dim[1],dim[2]);
			break;
		}
	}
	xdim[0]=dim[0];
	xdim[1]=dim[1];
	xdim[2]=dim[2];
	xdim[3]=xdim[0]*xdim[1]*xdim[2];

	//Read grid.
	gridX=(double *)malloc(sizeof(double)*xdim[3]);
	gridY=(double *)malloc(sizeof(double)*xdim[3]);
	int xi=0;
	while(fgets(line, maxchars, file) != NULL) {
		success = sscanf (line, "%lf %lf", &xcoord, &ycoord);
		if(success>0){
			if(xi>=xdim[3]){ printf("Warning: More grid points were available than specified.\n"); }
			gridX[xi]=xcoord;
			gridY[xi]=ycoord;
			xi++;
		}
	}
	fclose(file);
}
int **establishInterfaces(int totalblocks,int totalinterfaces,int totaledges,
		char **bcs,char **bcvalues,int **xdim){
	int b,otherb,s,i,n,othern,row;
	int **interfaces=mbinfo(totalinterfaces);
	i=0;row=totaledges;
	for(b=0;b<totalblocks;b++){
		for(s=0;s<4;s++){
			if(bcs[b][s]=='i'){
				otherb=bcvalues[b][s];
				if(otherb>=b){//add to interfaces
					interfaces[i][0]=b;
					interfaces[i][1]=s;
					interfaces[i][2]=bcvalues[b][s];
					n=getSideXdim(s,b,xdim);
					othern=getSideXdim(s,otherb,xdim);
					if(n!=othern){
						printf("Error: Interface dimension does not match (interface %d)!\n",i);
					}
					interfaces[i][3]=n;
					interfaces[i][4]=row;
					i++;
					row+=n;
				}
			}
		}
	}
	return interfaces;
}
int **readConnectivities(int *cmbinfo,int **xdim){
	int maxchars=200;
	int totalinterfaces;
	int refb,side,otherb;
	char line[maxchars];
	FILE * file = fopen ("connectivities", "rt");
	int success;
	//Read total number of interfaces.
	while(fgets(line, maxchars, file) != NULL){
		success = sscanf(line, "%d", &totalinterfaces);
		if(success>0){ break; }
	}
	int **interfaces=mbinfo(totalinterfaces);
	int interfaceid=0,interfacelength;
	int row=cmbinfo[1];
	//Read connectivity.
	while(fgets(line, maxchars, file) != NULL){
		success = sscanf (line, "%d %d %d", &refb, &side, &otherb);
		if(success>0){
			interfaces[interfaceid][0]=refb;
			interfaces[interfaceid][1]=side;
			interfaces[interfaceid][2]=otherb;
			if((side==0) | (side==2)){ interfacelength=xdim[refb][0]; }
			else { interfacelength=xdim[refb][1]; }
			interfaces[interfaceid][3]=interfacelength;
			interfaces[interfaceid][4]=row;
			row+=interfacelength;
			interfaceid++;
		}
	}
	return interfaces;
}
 */

#endif /* OLD_H_ */
