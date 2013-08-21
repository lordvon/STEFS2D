void fillxxBoundaryWallSub(Momentum*mm,Grid*g,int block,int side){
	int xi;
	Boundary*bo=&g->sidexx[block].b[side];
	for(xi=bo->st;xi<=bo->en;xi+=bo->iv){
		mm->uxx[block][xi]=0;
	}
	for(xi=bo->st;xi<=bo->en;xi+=bo->iv){
		mm->vxx[block][xi]=0;
	}
}
void fillxxBoundaryWall(Momentum*mm,Grid*g,BoundaryConditions*bc){
	int b,side;
	for(b=0;b<g->totalblocks;b++){
		for(side=0;side<4;side++){
			if(bc->id[b][side]=='w'){
				fillxxBoundaryWallSub(mm,g,b,side);
			}
		}
	}
}
void fillxxBoundaryFixedSub(Grid*g,Momentum*mm,int block,int side,double cart[2]){
	int xi;
	double a[2];
	Boundary*bo=&g->sidexx[block].b[side];
	for(xi=bo->st;xi<=bo->en;xi+=bo->iv){
		a[0]=g->a21xx[block][xi];
		a[1]=g->a22xx[block][xi];
		mm->uxx[block][xi]=cross(cart,a);
	}
	for(xi=bo->st;xi<=bo->en;xi+=bo->iv){
		a[0]=g->a11xx[block][xi];
		a[1]=g->a12xx[block][xi];
		mm->vxx[block][xi]=-cross(cart,a);
	}
}
void fillxxBoundaryFixed(Momentum*mm,Grid*g,BoundaryConditions*bc){
	int b,side;
	for(b=0;b<g->totalblocks;b++){
		for(side=0;side<4;side++){
			if(bc->id[b][side]=='f'){
				fillxxBoundaryFixedSub(g,mm,b,side,bc->fv[b]);
			}
		}
	}
}
void fillxxBoundaryZeroGradientSub(Grid*g,State*s,Momentum*mm,int block,int side){
	int xi,adj=0;
	Boundary*bo=&g->sidexx[block].b[side];
	if(side%2==0){
		int vi;
		if(side==2){ adj=-g->vdim[block][0]; }
		for(xi=bo->st;xi<=bo->en;xi+=bo->iv){
			vi=xi+adj;
			mm->uxx[block][xi]=s->u[block][vi];
		}
	} else {
		int hi;
		if(side==1){ adj=-1; }
		for(xi=bo->st;xi<=bo->en;xi+=bo->iv){
			hi=xi-xi/g->xdim[block][0]+adj;
			mm->vxx[block][xi]=s->v[block][hi];
		}
	}
}
void fillxxBoundaryZeroGradient(State*s,Momentum*mm,Grid*g,BoundaryConditions*bc){
	int b,side;
	for(b=0;b<g->totalblocks;b++){
		for(side=0;side<4;side++){
			if(bc->id[b][side]=='z'){
				fillxxBoundaryZeroGradientSub(g,s,mm,b,side);
			}
		}
	}
}
void fillxxBoundaryInterface02(State* s,Momentum* mm,Grid* g,Interfaces* is,int interface){
	int ip;
	int adj1,adj2;
	int b1,b2,side,xi1,xi2,vi1,vi2;
	b1=is->block[interface];
	side=is->side[interface];
	b2=is->adjacentBlock[interface];
	for(ip=0;ip<is->xsize[interface];ip++){
		xi1=is->blockx[interface][ip]-g->xstart[b1];
		xi2=is->adjacentBlockx[interface][ip]-g->xstart[b2];
		if(side==0){
			adj1=0;adj2=-g->xdim[b2][0];
		} else {
			adj1=-g->xdim[b1][0];adj2=0;
		}
		vi1=xi1+adj1;
		vi2=xi2+adj2;
		mm->uxx[b1][xi1]=mm->uxx[b2][xi2]=
				0.5*(s->u[b1][vi1]+s->u[b2][vi2]);
	}
}
void fillxxBoundaryInterface13(State* s,Momentum* mm,Grid* g,Interfaces* is,int interface){
	int ip;
	int adj1,adj2;
	int b1,b2,side,xi1,xi2,hi1,hi2;
	b1=is->block[interface];
	side=is->side[interface];
	b2=is->adjacentBlock[interface];
	for(ip=0;ip<is->xsize[interface];ip++){
		xi1=is->blockx[interface][ip]-g->xstart[b1];
		xi2=is->adjacentBlockx[interface][ip]-g->xstart[b2];
		if(side==1){
			adj1=-1;adj2=0;
		} else {
			adj1=0;adj2=-1;
		}
		hi1=xi1-xi1/g->xdim[b1][0]+adj1;
		hi2=xi2-xi2/g->xdim[b2][0]+adj2;
		mm->vxx[b1][xi1]=mm->vxx[b2][xi2]=
				0.5*(s->v[b1][hi1]+s->v[b2][hi2]);
	}
}
void fillxxBoundaryInterface(State* s,Momentum* mm,Grid* g,Interfaces* is){
	int i,side;
	for(i=0;i<is->totalinterfaces;i++){
		side=is->side[i];
		if(side%2==0){
			fillxxBoundaryInterface02(s,mm,g,is,i);
		} else {
			fillxxBoundaryInterface13(s,mm,g,is,i);
		}
	}
}
void fillxxBoundary(State*s,Momentum*mm,Grid*g,
		BoundaryConditions*bc,Interfaces*is){
	fillxxBoundaryZeroGradient(s,mm,g,bc);
	fillxxBoundaryFixed(mm,g,bc);
	fillxxBoundaryWall(mm,g,bc);
	fillxxBoundaryInterface(s,mm,g,is);
}
