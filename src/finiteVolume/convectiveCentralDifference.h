void fillI11InterfaceCentralDifference(Grid*g,Momentum*mm,Interfaces*is){
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
void fillI11InteriorCentralDifference(Grid* g,Momentum* mm){
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
}
void fillI22InterfaceCentralDifference(Grid*g,Momentum*mm,Interfaces*is){
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
void fillI22InteriorCentralDifference(Grid* g,Momentum* mm){
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
}
