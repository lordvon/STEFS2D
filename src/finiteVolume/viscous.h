void fillterms14and25(Property*p,Grid*g,Momentum*mm){
	int b,ci;
	double e1[2],e2[2],nueff;
	for(b=0;b<g->totalblocks;b++){
		for(ci=0;ci<g->cdim[b][3];ci++){
			nueff=mm->nut[b][ci]+p->nu;
			e1[0]=2*mm->uxcc[b][ci];
			e1[1]=e2[0]=mm->uycc[b][ci]+mm->vxcc[b][ci];
			e2[1]=2*mm->vycc[b][ci];
			mm->term14_1[b][ci]=nueff*(
					g->a22cc[b][ci]*e1[0]-
					g->a21cc[b][ci]*e2[0]);
			mm->term14_2[b][ci]=nueff*(
					g->a22cc[b][ci]*e1[1]-
					g->a21cc[b][ci]*e2[1]);
			mm->term25_1[b][ci]=nueff*(
					-g->a12cc[b][ci]*e1[0]+
					g->a11cc[b][ci]*e2[0]);
			mm->term25_2[b][ci]=nueff*(
					-g->a12cc[b][ci]*e1[1]+
					g->a11cc[b][ci]*e2[1]);
		}
	}
}
void fillterms15and24(Property*p,Grid* g,Momentum* mm){
	int b,xi;
	double e1[2],e2[2],nueff;
	for(b=0;b<g->totalblocks;b++){
		for(xi=0;xi<g->xdim[b][3];xi++){
			nueff=mm->nutxx[b][xi]+p->nu;
			e1[0]=2*mm->uxxx[b][xi];
			e1[1]=e2[0]=mm->uyxx[b][xi]+mm->vxxx[b][xi];
			e2[1]=2*mm->vyxx[b][xi];
			mm->term15_1[b][xi]=nueff*(
					-g->a12xx[b][xi]*e1[0]+
					g->a11xx[b][xi]*e2[0]);
			mm->term15_2[b][xi]=nueff*(
					-g->a12xx[b][xi]*e1[1]+
					g->a11xx[b][xi]*e2[1]);
			mm->term24_1[b][xi]=nueff*(
					g->a22xx[b][xi]*e1[0]-
					g->a21xx[b][xi]*e2[0]);
			mm->term24_2[b][xi]=nueff*(
					g->a22xx[b][xi]*e1[1]-
					g->a21xx[b][xi]*e2[1]);
		}
	}
}
void fillI14Interface(Grid*g,Momentum*mm,Interfaces*is){
	int i,b1,b2,vi1,vi2,ci1,ci2;
	double term1,term2,sign;
	int vi01,vi02,ci01,ci02,e,s1,s2,size;
	int viv1,viv2,civ1,civ2;
	for(i=0;i<is->totalinterfaces;i++){
		s1=is->side[i];
		if(s1%2==1){
			s2=(s1+2)%4;
			b1=is->block[i];
			b2=is->adjacentBlock[i];
			size=g->sidevv[b1].b[s1].size;
			vi01=g->sidevv[b1].b[s1].st;
			vi02=g->sidevv[b2].b[s2].st;
			ci01=g->sidecc[b1].b[s1].st;
			ci02=g->sidecc[b2].b[s2].st;
			viv1=g->sidevv[b1].b[s1].iv;
			viv2=g->sidevv[b2].b[s2].iv;
			civ1=g->sidecc[b1].b[s1].iv;
			civ2=g->sidecc[b2].b[s2].iv;
			if(s1==1){ sign=1; } else { sign=-1; }
			for(e=0;e<size;e++){
				vi1=vi01+e*viv1;
				vi2=vi02+e*viv2;
				ci1=ci01+e*civ1;
				ci2=ci02+e*civ2;
				term1=mm->term14_1[b1][ci1];
				term2=mm->term14_1[b2][ci2];
				mm->I14_1[b1][vi1]=mm->I14_1[b2][vi2]=
						sign*(term2-term1);
			}
			for(e=0;e<size;e++){
				vi1=vi01+e*viv1;
				vi2=vi02+e*viv2;
				ci1=ci01+e*civ1;
				ci2=ci02+e*civ2;
				term1=mm->term14_2[b1][ci1];
				term2=mm->term14_2[b2][ci2];
				mm->I14_2[b1][vi1]=mm->I14_2[b2][vi2]=
						sign*(term2-term1);
			}
		}
	}
}
void fillI14Zero(Momentum*mm,Grid*g,int block,int side){
	int vi;
	Boundary*bo=&g->sidevv[block].b[side];
	for(vi=bo->st;vi<=bo->en;vi+=bo->iv){
		mm->I14_1[block][vi]=0;
	}
	for(vi=bo->st;vi<=bo->en;vi+=bo->iv){
		mm->I14_2[block][vi]=0;
	}
}
void fillI14ZeroGradientSub(Momentum*mm,Grid*g,int block,int side){
	int vi,vi0,viv,vif,in;
	vi0=g->sidevv[block].b[side].st;
	viv=g->sidevv[block].b[side].iv;
	vif=g->sidevv[block].b[side].en;
	in=g->sidevv[block].b[side].in;
	for(vi=vi0;vi<=vif;vi+=viv){
		mm->I14_1[block][vi]=mm->I14_1[block][vi+in];
	}
	for(vi=vi0;vi<=vif;vi+=viv){
		mm->I14_2[block][vi]=mm->I14_2[block][vi+in];
	}
}
void fillI14ZeroGradient(Grid*g,Momentum*mm,BoundaryConditions*bc){
	int b,side;
	char sbc;
	for(b=0;b<g->totalblocks;b++){
		for(side=1;side<=3;side+=2){
			sbc=bc->id[b][side];
			if(sbc=='z'){
				fillI14ZeroGradientSub(mm,g,b,side);
			}
		}
	}
}
void fillI14Boundary(Grid*g,Momentum*mm,BoundaryConditions*bc,Interfaces*is){
	int b,side;
	char sbc;
	//Order matters: highest precedence last.
	//Zero Gradient.
	fillI14ZeroGradient(g,mm,bc);
	//Interface.
	fillI14Interface(g,mm,is);
	//Wall, fixed, and zero gradient.
	for(b=0;b<g->totalblocks;b++){
		for(side=1;side<=3;side+=2){
			sbc=bc->id[b][side];
			if((sbc=='f') | (sbc=='w')){
				fillI14Zero(mm,g,b,side);
			}
		}
	}
}
void fillI14(Grid*g,Momentum*mm,BoundaryConditions*bc,Interfaces*is){
	int b,i,j,vi,ci;
	int vd1,vd0;
	for(b=0;b<g->totalblocks;b++){
		vd1=g->vdim[b][1];
		vd0=g->vdim[b][0];
		for(j=0;j<vd1;j++){
			for(i=1;i<vd0-1;i++){
				vi=i+j*vd0;
				ci=vi-j;
				mm->I14_1[b][vi]=mm->term14_1[b][ci]-mm->term14_1[b][ci-1];
				mm->I14_2[b][vi]=mm->term14_2[b][ci]-mm->term14_2[b][ci-1];
			}
		}
	}
	fillI14Boundary(g,mm,bc,is);
}
void fillI15(Grid* g,Momentum* mm){
	int b,vi,xi;
	int xd0,vd3;
	for(b=0;b<g->totalblocks;b++){
		xd0=g->xdim[b][0];
		vd3=g->vdim[b][3];
		for(vi=0;vi<vd3;vi++){
			xi=vi;
			mm->I15_1[b][vi]=mm->term15_1[b][xi+xd0]-mm->term15_1[b][xi];
		}
		for(vi=0;vi<vd3;vi++){
			xi=vi;
			mm->I15_2[b][vi]=mm->term15_2[b][xi+xd0]-mm->term15_2[b][xi];
		}
	}
}
void fillI24(Grid* g,Momentum* mm){
	int b,hi,xi;
	int hd0,hd3;
	for(b=0;b<g->totalblocks;b++){
		hd0=g->hdim[b][0];
		hd3=g->hdim[b][3];
		for(hi=0;hi<hd3;hi++){
			xi=hi+hi/hd0;
			mm->I24_1[b][hi]=mm->term24_1[b][xi+1]-mm->term24_1[b][xi];
		}
		for(hi=0;hi<hd3;hi++){
			xi=hi+hi/hd0;
			mm->I24_2[b][hi]=mm->term24_2[b][xi+1]-mm->term24_2[b][xi];
		}
	}
}
void fillI25Interface(Grid* g,Momentum* mm,Interfaces* is){
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
			hi01=g->sidehh[b1].b[s1].st;
			hi02=g->sidehh[b2].b[s2].st;
			ci01=g->sidecc[b1].b[s1].st;
			ci02=g->sidecc[b2].b[s2].st;
			size=g->sidehh[b1].b[s1].size;
			if(s1==0){ sign=1; } else { sign=-1; }
			for(e=0;e<size;e++){
				hi1=hi01+e;
				hi2=hi02+e;
				ci1=ci01+e;
				ci2=ci02+e;
				term1=mm->term25_1[b1][ci1];
				term2=mm->term25_1[b2][ci2];
				mm->I25_1[b1][hi1]=mm->I25_1[b2][hi2]=
						sign*(term1-term2);
			}
			for(e=0;e<size;e++){
				hi1=hi01+e;
				hi2=hi02+e;
				ci1=ci01+e;
				ci2=ci02+e;
				term1=mm->term25_2[b1][ci1];
				term2=mm->term25_2[b2][ci2];
				mm->I25_2[b1][hi1]=mm->I25_2[b2][hi2]=
						sign*(term1-term2);
			}
		}
	}
}
void fillI25Zero(Momentum*mm,Grid*g,int block,int side){
	int hi;
	Boundary*bo=&g->sidehh[block].b[side];
	for(hi=bo->st;hi<=bo->en;hi+=bo->iv){
		mm->I25_1[block][hi]=0;
	}
	for(hi=bo->st;hi<=bo->en;hi+=bo->iv){
		mm->I25_2[block][hi]=0;
	}
}
void fillI25ZeroGradientSub(Momentum*mm,Grid*g,int block,int side){
	int hi,hi0,hiv,hif,in;
	hi0=g->sidehh[block].b[side].st;
	hiv=g->sidehh[block].b[side].iv;
	hif=g->sidehh[block].b[side].en;
	in=g->sidehh[block].b[side].in;
	for(hi=hi0;hi<=hif;hi+=hiv){
		mm->I25_1[block][hi]=mm->I25_1[block][hi+in];
	}
	for(hi=hi0;hi<=hif;hi+=hiv){
		mm->I25_2[block][hi]=mm->I25_2[block][hi+in];
	}
}
void fillI25ZeroGradient(Grid*g,Momentum*mm,BoundaryConditions*bc){
	int b,side;
	char sbc;
	for(b=0;b<g->totalblocks;b++){
		for(side=0;side<=2;side+=2){
			sbc=bc->id[b][side];
			if(sbc=='z'){
				fillI25ZeroGradientSub(mm,g,b,side);
			}
		}
	}
}
void fillI25Boundary(Grid*g,Momentum*mm,BoundaryConditions*bc,Interfaces*is){
	int b,side;
	char sbc;
	//Order matters: highest precedence last.
	//Zero Gradient.
	fillI25ZeroGradient(g,mm,bc);
	//Interface.
	fillI25Interface(g,mm,is);
	//Wall and Fixed.
	for(b=0;b<g->totalblocks;b++){
		for(side=0;side<=2;side+=2){
			sbc=bc->id[b][side];
			if((sbc=='f') | (sbc=='w')){
				fillI25Zero(mm,g,b,side);
			}
		}
	}
}
void fillI25(Grid* g,Momentum* mm,BoundaryConditions*bc,Interfaces*is){
	int b,i,j,hi,ci;
	int hd1,hd0,cd0;
	for(b=0;b<g->totalblocks;b++){
		hd1=g->hdim[b][1];
		hd0=g->hdim[b][0];
		cd0=g->cdim[b][0];
		for(j=1;j<hd1-1;j++){
			for(i=0;i<hd0;i++){
				hi=i+j*hd0;
				ci=hi;
				mm->I25_1[b][hi]=mm->term25_1[b][ci]-mm->term25_1[b][ci-cd0];
			}
		}
		for(j=1;j<hd1-1;j++){
			for(i=0;i<hd0;i++){
				hi=i+j*hd0;
				ci=hi;
				mm->I25_2[b][hi]=mm->term25_2[b][ci]-mm->term25_2[b][ci-cd0];
			}
		}
	}
	fillI25Boundary(g,mm,bc,is);
}
