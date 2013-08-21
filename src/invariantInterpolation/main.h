void fillcc(State* s,Momentum* mm,Grid* grid){
	int ci,vi,hi,b;
	int cd3,cd0,hd0;
	for(b=0;b<grid->totalblocks;b++){
		cd0=grid->cdim[b][0];
		cd3=grid->cdim[b][3];
		hd0=grid->hdim[b][0];
		for(ci=0;ci<cd3;ci++){
			vi=ci+ci/cd0;
			mm->ucc[b][ci]=0.5*(s->u[b][vi]+s->u[b][vi+1]);
		}
		for(ci=0;ci<cd3;ci++){
			hi=ci;
			mm->vcc[b][ci]=0.5*(s->v[b][hi]+s->v[b][hi+hd0]);
		}
	}
}
void fillxx(State* s,Momentum* mm,Grid* g,
		BoundaryConditions*bc,Interfaces* is){
	int xd1,xd0,vd0;
	int i,j,xi,vi,hi,b;
	for(b=0;b<g->totalblocks;b++){
		xd0=g->xdim[b][0];
		xd1=g->xdim[b][1];
		vd0=g->vdim[b][0];
		//Interior.
		for(j=0;j<xd1;j++){
			for(i=1;i<xd0-1;i++){
				xi=i+j*xd0;
				hi=xi-j;
				mm->vxx[b][xi]=0.5*(s->v[b][hi]+s->v[b][hi-1]);
			}
		}
		for(j=1;j<xd1-1;j++){
			for(i=0;i<xd0;i++){
				xi=i+j*xd0;
				vi=xi;
				mm->uxx[b][xi]=0.5*(s->u[b][vi]+s->u[b][vi-vd0]);
			}
		}
	}
	fillxxBoundary(s,mm,g,bc,is);
}
void filluhh(State* s,Grid* g,Momentum* mm,BoundaryConditions*bc,Interfaces* is){
	int b,i,j,hi,vi;
	int hd1,hd0,vd0;
	//interior.
	for(b=0;b<g->totalblocks;b++){
		hd1=g->hdim[b][1];
		hd0=g->hdim[b][0];
		vd0=g->vdim[b][0];
		for(j=1;j<hd1-1;j++){
			for(i=0;i<hd0;i++){
				hi=i+j*hd0;
				vi=hi+j;
				mm->uhh[b][hi]=0.25*(
						s->u[b][vi]+
						s->u[b][vi+1]+
						s->u[b][vi-vd0]+
						s->u[b][vi-vd0+1]);
			}
		}
	}
	filluhhBoundary(s,g,mm,bc,is);
}
void fillvvv(State* s,Grid* g,Momentum* mm,BoundaryConditions*bc,Interfaces* is){
	int b,i,j,hi,vi;
	int vd1,vd0,hd0;
	for(b=0;b<g->totalblocks;b++){
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		hd0=g->hdim[b][0];
		for(j=0;j<vd1;j++){
			for(i=1;i<vd0-1;i++){
				vi=i+j*vd0;
				hi=vi-j;
				mm->vvv[b][vi]=0.25*(
						s->v[b][hi]+
						s->v[b][hi+hd0]+
						s->v[b][hi-1]+
						s->v[b][hi-1+hd0]);
			}
		}
	}
	fillvvvBoundary(s,g,mm,bc,is);
}
