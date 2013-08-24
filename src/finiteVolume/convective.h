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
void fillI11(Grid* g,Momentum* mm,BoundaryConditions*bc,Interfaces* is){
	fillI11InteriorCentralDifference(g,mm);
	fillI11ZeroGradient(g,mm,bc);
	fillI11InterfaceCentralDifference(g,mm,is);
	//Fixed and Wall
	int b,side;
	char sbc;
	for(b=0;b<g->totalblocks;b++){
		for(side=1;side<=3;side+=2){
			sbc=bc->id[b][side];
			if((sbc=='f') | (sbc=='w')){
				fillI11Zero(mm,g,b,side);
			}
		}
	}
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
void fillI22(Grid* g,Momentum* mm,BoundaryConditions*bc,Interfaces* is){
	fillI22InteriorCentralDifference(g,mm);
	fillI22ZeroGradient(g,mm,bc);
	fillI22InterfaceCentralDifference(g,mm,is);
	//Fixed and Walls
	int b,side;
	char sbc;
	for(b=0;b<g->totalblocks;b++){
		for(side=0;side<=2;side+=2){
			sbc=bc->id[b][side];
			if((sbc=='f') | (sbc=='w')){
				fillI22Zero(mm,g,b,side);
			}
		}
	}
}
