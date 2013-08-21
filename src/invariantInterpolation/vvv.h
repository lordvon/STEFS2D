void fillvvvInterface(State* s,Grid* g,Momentum* mm,Interfaces* is){
	int i,b1,b2,s1,s2,e;
	int vi1,vi2,hi1,hi2;
	int hi01,hi02,vi01,vi02,size;
	int hiv1,hiv2,viv1,viv2;
	for(i=0;i<is->totalinterfaces;i++){
		s1=is->side[i];
		if(s1%2==1){
			s2=(s1+2)%4;
			b1=is->block[i];
			b2=is->adjacentBlock[i];
			size=g->sidevv[b1].b[s1].size;
			hi01=g->sidehh[b1].b[s1].st;
			hi02=g->sidehh[b2].b[s2].st;
			vi01=g->sidevv[b1].b[s1].st;
			vi02=g->sidevv[b2].b[s2].st;
			hiv1=g->sidehh[b1].b[s1].iv;
			hiv2=g->sidehh[b2].b[s2].iv;
			viv1=g->sidevv[b1].b[s1].iv;
			viv2=g->sidevv[b2].b[s2].iv;
			for(e=0;e<size;e++){
				hi1=hi01+e*hiv1;
				hi2=hi02+e*hiv2;
				vi1=vi01+e*viv1;
				vi2=vi02+e*viv2;
				mm->vvv[b1][vi1]=mm->vvv[b2][vi2]=0.25*(
						s->v[b1][hi1]+
						s->v[b1][hi1+hiv1]+
						s->v[b2][hi2]+
						s->v[b2][hi2+hiv2]);
			}
		}
	}
}
void fillvvvBoundaryZeroGradientSub(State*s,Momentum*mm,Grid*g,int block,int side){
	int e,vi,hi;
	int hi0,vi0,hiv,viv,size;
	size=g->sidevv[block].b[side].size;
	vi0=g->sidevv[block].b[side].st;
	viv=g->sidevv[block].b[side].iv;
	hi0=g->sidehh[block].b[side].st;
	hiv=g->sidehh[block].b[side].iv;
	for(e=0;e<size;e++){
		vi=vi0+e*viv;
		hi=hi0+e*hiv;
		mm->vvv[block][vi]=0.5*(
				s->v[block][hi]+
				s->v[block][hi+hiv]);
	}
}
void fillvvvBoundaryZeroGradient(State*s,Momentum*mm,Grid*g,BoundaryConditions*bc){
	int b,side;
	char sbc;
	for(b=0;b<g->totalblocks;b++){
		for(side=1;side<=3;side+=2){
			sbc=bc->id[b][side];
			if(sbc=='z'){
				fillvvvBoundaryZeroGradientSub(s,mm,g,b,side);
			}
		}
	}
}
void fillvvvBoundaryFixedSub(Momentum*mm,Grid*g,int block,int side,double*cart){
	int vi;
	double a[2];
	Boundary*bo=&g->sidevv[block].b[side];
	for(vi=bo->st;vi<=bo->en;vi+=bo->iv){
		a[0]=g->a11vv[block][vi];
		a[1]=g->a12vv[block][vi];
		mm->vvv[block][vi]=-cross(cart,a);
	}
}
void fillvvvBoundaryFixed(Momentum*mm,Grid*g,BoundaryConditions*bc){
	int b,side;
	char sbc;
	for(b=0;b<g->totalblocks;b++){
		for(side=1;side<=3;side+=2){
			sbc=bc->id[b][side];
			if(sbc=='f'){
				fillvvvBoundaryFixedSub(mm,g,b,side,bc->fv[b]);
			}
		}
	}
}
void fillvvvBoundaryWallSub(Momentum*mm,Grid*g,int block,int side){
	int vi;
	Boundary*bo=&g->sidevv[block].b[side];
	for(vi=bo->st;vi<=bo->en;vi+=bo->iv){
		mm->vvv[block][vi]=0;
	}
}
void fillvvvBoundaryWall(Momentum*mm,Grid*g,BoundaryConditions*bc){
	int b,side;
	char sbc;
	for(b=0;b<g->totalblocks;b++){
		for(side=1;side<=3;side+=2){
			sbc=bc->id[b][side];
			if(sbc=='w'){
				fillvvvBoundaryWallSub(mm,g,b,side);
			}
		}
	}
}
void fillvvvBoundary(State*s,Grid*g,Momentum*mm,BoundaryConditions*bc,Interfaces*is){
	fillvvvBoundaryZeroGradient(s,mm,g,bc);
	fillvvvInterface(s,g,mm,is);
	fillvvvBoundaryFixed(mm,g,bc);
	fillvvvBoundaryWall(mm,g,bc);
}
