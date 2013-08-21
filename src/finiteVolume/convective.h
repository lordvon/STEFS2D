void fillI11Interface(Grid* g,Momentum* mm,Interfaces* is){
	int i,e,size;
	int b1,b2,s1,s2,vi1,vi2,ci1,ci2;
	int vi01,vi02,viv1,viv2,ci01,ci02,civ1,civ2;
	double term1,term2,sign;
	for(i=0;i<is->totalinterfaces;i++){
		s1=is->side[i];
		if(s1%2==1){
			s2=(s1+2)%4;
			b1=is->block[i];
			b2=is->adjacentBlock[i];
			size=g->sidevv[b1].b[s1].size;
			vi01=g->sidevv[b1].b[s1].st;
			vi02=g->sidevv[b2].b[s2].st;
			viv1=g->sidevv[b1].b[s1].iv;
			viv2=g->sidevv[b2].b[s2].iv;
			ci01=g->sidecc[b1].b[s1].st;
			ci02=g->sidecc[b2].b[s2].st;
			civ1=g->sidecc[b1].b[s1].iv;
			civ2=g->sidecc[b2].b[s2].iv;
			if(s1==1){ sign=1; } else { sign=-1; }
			for(e=0;e<size;e++){
				vi1=vi01+e*viv1;
				vi2=vi02+e*viv2;
				ci1=ci01+e*civ1;
				ci2=ci02+e*civ2;
				term1=mm->ucartcc[b1][ci1]*mm->ucc[b1][ci1];
				term2=mm->ucartcc[b2][ci2]*mm->ucc[b2][ci2];
				mm->I11_1[b1][vi1]=mm->I11_1[b2][vi2]=sign*(term2-term1);
			}
			for(e=0;e<size;e++){
				vi1=vi01+e*viv1;
				vi2=vi02+e*viv2;
				ci1=ci01+e*civ1;
				ci2=ci02+e*civ2;
				term1=mm->vcartcc[b1][ci1]*mm->ucc[b1][ci1];
				term2=mm->vcartcc[b2][ci2]*mm->ucc[b2][ci2];
				mm->I11_2[b1][vi1]=mm->I11_2[b2][vi2]=sign*(term2-term1);
			}
		}
	}
}
void fillI11Zero(Momentum*mm,Grid*g,int block,int side){
	int vi;
	Boundary*bo=&g->sidevv[block].b[side];
	for(vi=bo->st;vi<=bo->en;vi+=bo->iv){
		mm->I11_1[block][vi]=0;
	}
	for(vi=bo->st;vi<=bo->en;vi+=bo->iv){
		mm->I11_2[block][vi]=0;
	}
}
void fillI11ZeroGradientSub(Momentum*mm,Grid*g,int block,int side){
	int vi,vi0,viv,vif,in;
	vi0=g->sidevv[block].b[side].st;
	viv=g->sidevv[block].b[side].iv;
	vif=g->sidevv[block].b[side].en;
	in=g->sidevv[block].b[side].in;
	for(vi=vi0;vi<=vif;vi+=viv){
		mm->I11_1[block][vi]=mm->I11_1[block][vi+in];
	}
	for(vi=vi0;vi<=vif;vi+=viv){
		mm->I11_2[block][vi]=mm->I11_2[block][vi+in];
	}
}
void fillI11ZeroGradient(Grid*g,Momentum*mm,BoundaryConditions*bc){
	int b,side;
	char sbc;
	for(b=0;b<g->totalblocks;b++){
		for(side=1;side<=3;side+=2){
			sbc=bc->id[b][side];
			if(sbc=='z'){
				fillI11ZeroGradientSub(mm,g,b,side);
			}
		}
	}
}
void fillI11Boundary(Grid*g,Momentum*mm,BoundaryConditions*bc,Interfaces*is){
	int b,side;
	char sbc;
	//Zero Gradient
	fillI11ZeroGradient(g,mm,bc);
	//Interface
	fillI11Interface(g,mm,is);
	//Fixed and Wall
	for(b=0;b<g->totalblocks;b++){
		for(side=1;side<=3;side+=2){
			sbc=bc->id[b][side];
			if((sbc=='f') | (sbc=='w')){
				fillI11Zero(mm,g,b,side);
			}
		}
	}
}
void fillI11(Grid* g,Momentum* mm,BoundaryConditions*bc,Interfaces* is){
	int b,i,j,ci,vi;
	int vd1,vd0;
	for(b=0;b<g->totalblocks;b++){
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		for(j=0;j<vd1;j++){
			for(i=1;i<vd0-1;i++){
				vi=i+j*vd0;
				ci=vi-j;
				mm->I11_1[b][vi]=mm->ucartcc[b][ci]*mm->ucc[b][ci]-
						mm->ucartcc[b][ci-1]*mm->ucc[b][ci-1];
			}
		}
		for(j=0;j<vd1;j++){
			for(i=1;i<vd0-1;i++){
				vi=i+j*vd0;
				ci=vi-j;
				mm->I11_2[b][vi]=mm->vcartcc[b][ci]*mm->ucc[b][ci]-
						mm->vcartcc[b][ci-1]*mm->ucc[b][ci-1];
			}
		}
	}
	fillI11Boundary(g,mm,bc,is);
}
void fillI12(Grid* g,Momentum* mm){
	int b,xi,vi;
	int xd0,vd3;
	for(b=0;b<g->totalblocks;b++){
		vd3=g->vdim[b][3];
		xd0=g->xdim[b][0];
		for(vi=0;vi<vd3;vi++){
			xi=vi;
			mm->I12_1[b][vi]=mm->ucartxx[b][xi+xd0]*mm->vxx[b][xi+xd0]-
					mm->ucartxx[b][xi]*mm->vxx[b][xi];
		}
		for(vi=0;vi<vd3;vi++){
			xi=vi;
			mm->I12_2[b][vi]=mm->vcartxx[b][xi+xd0]*mm->vxx[b][xi+xd0]-
					mm->vcartxx[b][xi]*mm->vxx[b][xi];
		}
	}
}
void fillI21(Grid* g,Momentum* mm){
	int b,xi,hi;
	int hd0,hd3;
	for(b=0;b<g->totalblocks;b++){
		hd0=g->hdim[b][0];
		hd3=g->hdim[b][3];
		for(hi=0;hi<hd3;hi++){
			xi=hi+hi/hd0;
			mm->I21_1[b][hi]=mm->ucartxx[b][xi+1]*mm->uxx[b][xi+1]-
					mm->ucartxx[b][xi]*mm->uxx[b][xi];
		}
		for(hi=0;hi<hd3;hi++){
			xi=hi+hi/hd0;
			mm->I21_2[b][hi]=mm->vcartxx[b][xi+1]*mm->uxx[b][xi+1]-
					mm->vcartxx[b][xi]*mm->uxx[b][xi];
		}
	}
}
void fillI22Interface(Grid* g,Momentum* mm,Interfaces* is){
	int i,e,size;
	int b1,b2,s1,s2,hi1,hi2,ci1,ci2;
	int hi01,hi02,ci01,ci02;
	double term1,term2,sign;
	for(i=0;i<is->totalinterfaces;i++){
		s1=is->side[i];
		if(s1%2==0){
			s2=(s1+2)%4;
			b1=is->block[i];
			b2=is->adjacentBlock[i];
			size=g->sidehh[b1].b[s1].size;
			hi01=g->sidehh[b1].b[s1].st;
			hi02=g->sidehh[b2].b[s2].st;
			ci01=g->sidecc[b1].b[s1].st;
			ci02=g->sidecc[b2].b[s2].st;
			if(s1==0){ sign=1; } else { sign=-1; }
			for(e=0;e<size;e++){
				hi1=hi01+e;
				hi2=hi02+e;
				ci1=ci01+e;
				ci2=ci02+e;
				term1=mm->ucartcc[b1][ci1]*mm->vcc[b1][ci1];
				term2=mm->ucartcc[b2][ci2]*mm->vcc[b2][ci2];
				mm->I22_1[b1][hi1]=mm->I22_1[b2][hi2]=sign*(term1-term2);
			}
			for(e=0;e<size;e++){
				hi1=hi01+e;
				hi2=hi02+e;
				ci1=ci01+e;
				ci2=ci02+e;
				term1=mm->vcartcc[b1][ci1]*mm->vcc[b1][ci1];
				term2=mm->vcartcc[b2][ci2]*mm->vcc[b2][ci2];
				mm->I22_2[b1][hi1]=mm->I22_2[b2][hi2]=sign*(term1-term2);
			}
		}
	}
}
void fillI22Zero(Momentum*mm,Grid*g,int block,int side){
	int hi;
	Boundary*bo=&g->sidehh[block].b[side];
	for(hi=bo->st;hi<=bo->en;hi+=bo->iv){
		mm->I22_1[block][hi]=0;
	}
	for(hi=bo->st;hi<=bo->en;hi+=bo->iv){
		mm->I22_2[block][hi]=0;
	}
}
void fillI22ZeroGradientSub(Momentum*mm,Grid*g,int block,int side){
	int hi,hi0,hiv,hif,in;
	hi0=g->sidehh[block].b[side].st;
	hiv=g->sidehh[block].b[side].iv;
	hif=g->sidehh[block].b[side].en;
	in=g->sidehh[block].b[side].in;
	for(hi=hi0;hi<=hif;hi+=hiv){
		mm->I22_1[block][hi]=mm->I22_1[block][hi+in];
	}
	for(hi=hi0;hi<=hif;hi+=hiv){
		mm->I22_2[block][hi]=mm->I22_2[block][hi+in];
	}
}
void fillI22ZeroGradient(Grid*g,Momentum*mm,BoundaryConditions*bc){
	int b,side;
	char sbc;
	for(b=0;b<g->totalblocks;b++){
		for(side=0;side<=2;side+=2){
			sbc=bc->id[b][side];
			if(sbc=='z'){
				fillI22ZeroGradientSub(mm,g,b,side);
			}
		}
	}
}
void fillI22Boundary(Grid*g,Momentum*mm,BoundaryConditions*bc,Interfaces*is){
	int b,side;
	char sbc;
	//Order matters: highest precedence last.
	//Zero Gradient
	fillI22ZeroGradient(g,mm,bc);
	//Interface
	fillI22Interface(g,mm,is);
	//Boundary
	for(b=0;b<g->totalblocks;b++){
		for(side=0;side<=2;side+=2){
			sbc=bc->id[b][side];
			if((sbc=='f') | (sbc=='w')){
				fillI22Zero(mm,g,b,side);
			}
		}
	}
}
void fillI22(Grid* g,Momentum* mm,BoundaryConditions*bc,Interfaces* is){
	int b,i,j,ci,hi;
	int hd1,hd0,cd0;
	for(b=0;b<g->totalblocks;b++){
		hd0=g->hdim[b][0];
		hd1=g->hdim[b][1];
		cd0=g->cdim[b][0];
		for(j=1;j<hd1-1;j++){
			for(i=0;i<hd0;i++){
				hi=i+j*hd0;
				ci=hi;
				mm->I22_1[b][hi]=mm->ucartcc[b][ci]*mm->vcc[b][ci]-
						mm->ucartcc[b][ci-cd0]*mm->vcc[b][ci-cd0];
			}
		}
		for(j=1;j<hd1-1;j++){
			for(i=0;i<hd0;i++){
				hi=i+j*hd0;
				ci=hi;
				mm->I22_2[b][hi]=mm->vcartcc[b][ci]*mm->vcc[b][ci]-
						mm->vcartcc[b][ci-cd0]*mm->vcc[b][ci-cd0];
			}
		}
	}
	fillI22Boundary(g,mm,bc,is);
}
