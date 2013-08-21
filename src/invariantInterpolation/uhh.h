void filluhhInterface(State* s,Grid* g,Momentum* mm,Interfaces* is){
	int i,b1,b2,s1,s2,e;
	int vi1,vi2,hi1,hi2;
	int hi01,hi02,vi01,vi02,size;
	for(i=0;i<is->totalinterfaces;i++){
		s1=is->side[i];
		if(s1%2==0){
			s2=(s1+2)%4;
			b1=is->block[i];
			b2=is->adjacentBlock[i];
			size=g->sidehh[b1].b[s1].size;
			hi01=g->sidehh[b1].b[s1].st;
			hi02=g->sidehh[b2].b[s2].st;
			vi01=g->sidevv[b1].b[s1].st;
			vi02=g->sidevv[b2].b[s2].st;
			for(e=0;e<size;e++){
				hi1=hi01+e;
				hi2=hi02+e;
				vi1=vi01+e;
				vi2=vi02+e;
				mm->uhh[b1][hi1]=mm->uhh[b2][hi2]=0.25*(
								s->u[b1][vi1]+
								s->u[b1][vi1+1]+
								s->u[b2][vi2]+
								s->u[b2][vi2+1]);
			}
		}
	}
}
void filluhhBoundaryZeroGradientSub(State*s,Momentum*mm,Grid*g,int block,int side){
	int e,vi,hi;
	int hi0,vi0,size;
	size=g->sidehh[block].b[side].size;
	vi0=g->sidevv[block].b[side].st;
	hi0=g->sidehh[block].b[side].st;
	for(e=0;e<size;e++){
		vi=vi0+e;
		hi=hi0+e;
		mm->uhh[block][hi]=0.5*(
				s->u[block][vi]+
				s->u[block][vi+1]);
	}
}
void filluhhBoundaryZeroGradient(State*s,Momentum*mm,Grid*g,BoundaryConditions*bc){
	int b,side;
	char sbc;
	for(b=0;b<g->totalblocks;b++){
		for(side=0;side<=2;side+=2){
			sbc=bc->id[b][side];
			if(sbc=='z'){
				filluhhBoundaryZeroGradientSub(s,mm,g,b,side);
			}
		}
	}
}
void filluhhBoundaryFixedSub(Momentum*mm,Grid*g,int block,int side,double*cart){
	int hi;
	double a[2];
	Boundary*bo=&g->sidehh[block].b[side];
	for(hi=bo->st;hi<=bo->en;hi+=bo->iv){
		a[0]=g->a21hh[block][hi];
		a[1]=g->a22hh[block][hi];
		mm->uhh[block][hi]=cross(cart,a);
	}
}
void filluhhBoundaryFixed(Momentum*mm,Grid*g,BoundaryConditions*bc){
	int b,side;
	char sbc;
	for(b=0;b<g->totalblocks;b++){
		for(side=0;side<=2;side+=2){
			sbc=bc->id[b][side];
			if(sbc=='f'){
				filluhhBoundaryFixedSub(mm,g,b,side,bc->fv[b]);
			}
		}
	}
}
void filluhhBoundaryWallSub(Momentum*mm,Grid*g,int block,int side){
	int hi;
	Boundary*bo=&g->sidehh[block].b[side];
	for(hi=bo->st;hi<=bo->en;hi+=bo->iv){
		mm->uhh[block][hi]=0;
	}
}
void filluhhBoundaryWall(Momentum*mm,Grid*g,BoundaryConditions*bc){
	int b,side;
	char sbc;
	for(b=0;b<g->totalblocks;b++){
		for(side=0;side<=2;side+=2){
			sbc=bc->id[b][side];
			if(sbc=='w'){
				filluhhBoundaryWallSub(mm,g,b,side);
			}
		}
	}
}
void filluhhBoundary(State*s,Grid*g,Momentum*mm,BoundaryConditions*bc,Interfaces*is){
	filluhhBoundaryZeroGradient(s,mm,g,bc);
	filluhhInterface(s,g,mm,is);
	filluhhBoundaryFixed(mm,g,bc);
	filluhhBoundaryWall(mm,g,bc);
}
