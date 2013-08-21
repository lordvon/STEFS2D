void fillCartXixx(Grid* g,Momentum* mm){
	int b,i,j,xi,hi;
	int xd1,xd0;
	for(b=0;b<g->totalblocks;b++){
		xd1=g->xdim[b][1];
		xd0=g->xdim[b][0];
		for(j=0;j<xd1;j++){
			for(i=0;i<xd0;i++){
				xi=i+j*xd0;
				hi=xi-j;
				if(i==0){
					mm->ucartxi[b][xi]=mm->ucartxx[b][xi+1]-mm->ucartxx[b][xi];
				} else if(i==xd0-1){
					mm->ucartxi[b][xi]=mm->ucartxx[b][xi]-mm->ucartxx[b][xi-1];
				} else {
					mm->ucartxi[b][xi]=mm->ucarthh[b][hi]-mm->ucarthh[b][hi-1];
				}
			}
		}
		for(j=0;j<xd1;j++){
			for(i=0;i<xd0;i++){
				xi=i+j*xd0;
				hi=xi-j;
				if(i==0){
					mm->vcartxi[b][xi]=mm->vcartxx[b][xi+1]-mm->vcartxx[b][xi];
				} else if(i==xd0-1){
					mm->vcartxi[b][xi]=mm->vcartxx[b][xi]-mm->vcartxx[b][xi-1];
				} else {
					mm->vcartxi[b][xi]=mm->vcarthh[b][hi]-mm->vcarthh[b][hi-1];
				}
			}
		}
	}
}
void fillCartEtaxx(Grid* g,Momentum* mm){
	int b,i,j,xi,vi;
	int xd1,xd0,vd0;
	for(b=0;b<g->totalblocks;b++){
		xd1=g->xdim[b][1];
		xd0=g->xdim[b][0];
		vd0=g->vdim[b][0];
		for(j=0;j<xd1;j++){
			for(i=0;i<xd0;i++){
				xi=i+j*xd0;
				vi=xi;
				if(j==0){
					mm->ucarteta[b][xi]=mm->ucartxx[b][xi+xd0]-mm->ucartxx[b][xi];
					mm->vcarteta[b][xi]=mm->vcartxx[b][xi+xd0]-mm->vcartxx[b][xi];
				} else if(j==xd1-1){
					mm->ucarteta[b][xi]=mm->ucartxx[b][xi]-mm->ucartxx[b][xi-xd0];
					mm->vcarteta[b][xi]=mm->vcartxx[b][xi]-mm->vcartxx[b][xi-xd0];
				} else {
					mm->ucarteta[b][xi]=mm->ucartvv[b][vi]-mm->ucartvv[b][vi-vd0];
					mm->vcarteta[b][xi]=mm->vcartvv[b][vi]-mm->vcartvv[b][vi-vd0];
				}
			}
		}
	}
}
void fillCartXixxInterface(Grid* g,Momentum* mm,Interfaces* is){
	int i,xi1,xi2,hi1,hi2;
	int b1,b2,s1,s2;
	int e,hi01,hi02,xi01,xi02,size;
	int hiv1,hiv2,xiv1,xiv2;
	double sign;
	for(i=0;i<is->totalinterfaces;i++){
		s1=is->side[i];
		if(s1%2==1){
			b1=is->block[i];
			b2=is->adjacentBlock[i];
			s2=(s1+2)%4;
			size=g->sidexx[b1].b[s1].size;
			hi01=g->sidehh[b1].b[s1].st;
			hi02=g->sidehh[b2].b[s2].st;
			xi01=g->sidexx[b1].b[s1].st;
			xi02=g->sidexx[b2].b[s2].st;
			hiv1=g->sidehh[b1].b[s1].iv;
			hiv2=g->sidehh[b2].b[s2].iv;
			xiv1=g->sidexx[b1].b[s1].iv;
			xiv2=g->sidexx[b2].b[s2].iv;
			if(s1==1){ sign=1; } else { sign=-1; }
			for(e=0;e<size;e++){
				hi1=hi01+e*hiv1;
				hi2=hi02+e*hiv2;
				xi1=xi01+e*xiv1;
				xi2=xi02+e*xiv2;
				mm->ucartxi[b1][xi1]=mm->ucartxi[b2][xi2]=sign*(
						mm->ucarthh[b2][hi2]-mm->ucarthh[b1][hi1]);
			}
			for(e=0;e<size;e++){
				hi1=hi01+e*hiv1;
				hi2=hi02+e*hiv2;
				xi1=xi01+e*xiv1;
				xi2=xi02+e*xiv2;
				mm->vcartxi[b1][xi1]=mm->vcartxi[b2][xi2]=sign*(
						mm->vcarthh[b2][hi2]-mm->vcarthh[b1][hi1]);
			}
		}
	}
}
void fillCartEtaxxInterface(Grid* g,Momentum* mm,Interfaces* is){
	int i,xi1,xi2,vi1,vi2,b1,b2,s1,s2;
	int e,vi01,vi02,xi01,xi02,size;
	double sign;
	for(i=0;i<is->totalinterfaces;i++){
		s1=is->side[i];
		if(s1%2==0){
			b1=is->block[i];
			b2=is->adjacentBlock[i];
			s2=(s1+2)%4;
			size=g->sidexx[b1].b[s1].size;
			vi01=g->sidevv[b1].b[s1].st;
			vi02=g->sidevv[b2].b[s2].st;
			xi01=g->sidexx[b1].b[s1].st;
			xi02=g->sidexx[b2].b[s2].st;
			if(s1==0){ sign=1; } else { sign=-1; }
			for(e=0;e<size;e++){
				vi1=vi01+e;
				vi2=vi02+e;
				xi1=xi01+e;
				xi2=xi02+e;
				mm->ucarteta[b1][xi1]=mm->ucarteta[b2][xi2]=sign*(
						mm->ucartvv[b1][vi1]-mm->ucartvv[b2][vi2]);
			}
			for(e=0;e<size;e++){
				vi1=vi01+e;
				vi2=vi02+e;
				xi1=xi01+e;
				xi2=xi02+e;
				mm->vcarteta[b1][xi1]=mm->vcarteta[b2][xi2]=sign*(
						mm->vcartvv[b1][vi1]-mm->vcartvv[b2][vi2]);
			}
		}
	}
}
