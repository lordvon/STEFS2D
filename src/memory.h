void freeGrid(Grid*g){
	int b;
	for(b=0;b<g->totalblocks;b++){
		free(g->x[b]);
		free(g->y[b]);
		free(g->xdim[b]);
		free(g->cdim[b]);
		free(g->vdim[b]);
		free(g->hdim[b]);
		free(g->blockNames[b]);
	}
	free(g->x);
	free(g->y);
	free(g->xdim);
	free(g->cdim);
	free(g->vdim);
	free(g->hdim);
	free(g->xstart);
	free(g->cstart);
	free(g->estart);
	free(g->blockNames);
	free(g->sidexx);
	free(g->sidevv);
	free(g->sidehh);
	free(g->sidecc);
}
void mallocMomentum(State* s,Grid* grid,Momentum* mm,Interfaces*is){
	int b;
	int vd3,hd3,cd3,xd3;
	mm->out=fopen("out/momentum.out","w");
	int uvecsize=grid->totaledges+is->totalinterfacepoints;
	s->uvec=malloc(uvecsize*sizeof(double));
	s->u=malloc(grid->totalblocks*sizeof(double*));
	s->v=malloc(grid->totalblocks*sizeof(double*));
	s->sanu=malloc(grid->totalblocks*sizeof(double*));
	s->ut=malloc(grid->totalblocks*sizeof(double*));
	s->vt=malloc(grid->totalblocks*sizeof(double*));
	s->sanut=malloc(grid->totalblocks*sizeof(double*));
	mm->nut=malloc(grid->totalblocks*sizeof(double*));
	mm->ucc=malloc(grid->totalblocks*sizeof(double*));
	mm->vcc=malloc(grid->totalblocks*sizeof(double*));
	mm->uxx=malloc(grid->totalblocks*sizeof(double*));
	mm->vxx=malloc(grid->totalblocks*sizeof(double*));
	mm->uhh=malloc(grid->totalblocks*sizeof(double*));
	mm->vvv=malloc(grid->totalblocks*sizeof(double*));
	mm->uxcc=malloc(grid->totalblocks*sizeof(double*));
	mm->uycc=malloc(grid->totalblocks*sizeof(double*));
	mm->vxcc=malloc(grid->totalblocks*sizeof(double*));
	mm->vycc=malloc(grid->totalblocks*sizeof(double*));
	mm->ucartxi=malloc(grid->totalblocks*sizeof(double*));
	mm->ucarteta=malloc(grid->totalblocks*sizeof(double*));
	mm->vcartxi=malloc(grid->totalblocks*sizeof(double*));
	mm->vcarteta=malloc(grid->totalblocks*sizeof(double*));
	mm->uxxx=malloc(grid->totalblocks*sizeof(double*));
	mm->vxxx=malloc(grid->totalblocks*sizeof(double*));
	mm->uyxx=malloc(grid->totalblocks*sizeof(double*));
	mm->vyxx=malloc(grid->totalblocks*sizeof(double*));
	mm->ucartcc=malloc(grid->totalblocks*sizeof(double*));
	mm->vcartcc=malloc(grid->totalblocks*sizeof(double*));
	mm->ucartxx=malloc(grid->totalblocks*sizeof(double*));
	mm->vcartxx=malloc(grid->totalblocks*sizeof(double*));
	mm->ucartvv=malloc(grid->totalblocks*sizeof(double*));
	mm->vcartvv=malloc(grid->totalblocks*sizeof(double*));
	mm->ucarthh=malloc(grid->totalblocks*sizeof(double*));
	mm->vcarthh=malloc(grid->totalblocks*sizeof(double*));
	mm->term14_1=malloc(grid->totalblocks*sizeof(double*));
	mm->term14_2=malloc(grid->totalblocks*sizeof(double*));
	mm->term15_1=malloc(grid->totalblocks*sizeof(double*));
	mm->term15_2=malloc(grid->totalblocks*sizeof(double*));
	mm->term24_1=malloc(grid->totalblocks*sizeof(double*));
	mm->term24_2=malloc(grid->totalblocks*sizeof(double*));
	mm->term25_1=malloc(grid->totalblocks*sizeof(double*));
	mm->term25_2=malloc(grid->totalblocks*sizeof(double*));
	mm->I11_1=malloc(grid->totalblocks*sizeof(double*));
	mm->I11_2=malloc(grid->totalblocks*sizeof(double*));
	mm->I12_1=malloc(grid->totalblocks*sizeof(double*));
	mm->I12_2=malloc(grid->totalblocks*sizeof(double*));
	mm->I14_1=malloc(grid->totalblocks*sizeof(double*));
	mm->I14_2=malloc(grid->totalblocks*sizeof(double*));
	mm->I15_1=malloc(grid->totalblocks*sizeof(double*));
	mm->I15_2=malloc(grid->totalblocks*sizeof(double*));
	mm->I21_1=malloc(grid->totalblocks*sizeof(double*));
	mm->I21_2=malloc(grid->totalblocks*sizeof(double*));
	mm->I22_1=malloc(grid->totalblocks*sizeof(double*));
	mm->I22_2=malloc(grid->totalblocks*sizeof(double*));
	mm->I24_1=malloc(grid->totalblocks*sizeof(double*));
	mm->I24_2=malloc(grid->totalblocks*sizeof(double*));
	mm->I25_1=malloc(grid->totalblocks*sizeof(double*));
	mm->I25_2=malloc(grid->totalblocks*sizeof(double*));
	mm->nutxx=malloc(grid->totalblocks*sizeof(double*));
	for(b=0;b<grid->totalblocks;b++){
		vd3=grid->vdim[b][3];
		hd3=grid->hdim[b][3];
		cd3=grid->cdim[b][3];
		xd3=grid->xdim[b][3];
		//s->u[b]=malloc(vd3*sizeof(double));
		//s->v[b]=malloc(hd3*sizeof(double));
		s->u[b]=&s->uvec[grid->estart[b]];
		s->v[b]=&s->uvec[grid->estart[b]+vd3];
		s->sanu[b]=malloc(cd3*sizeof(double));
		s->ut[b]=malloc(vd3*sizeof(double));
		s->vt[b]=malloc(hd3*sizeof(double));
		s->sanut[b]=malloc(cd3*sizeof(double));
		mm->nut[b]=malloc(cd3*sizeof(double));
		mm->ucc[b]=malloc(cd3*sizeof(double));
		mm->vcc[b]=malloc(cd3*sizeof(double));
		mm->uxx[b]=malloc(xd3*sizeof(double));
		mm->vxx[b]=malloc(xd3*sizeof(double));
		mm->uhh[b]=malloc(hd3*sizeof(double));
		mm->vvv[b]=malloc(vd3*sizeof(double));
		mm->uxcc[b]=malloc(cd3*sizeof(double));
		mm->uycc[b]=malloc(cd3*sizeof(double));
		mm->vxcc[b]=malloc(cd3*sizeof(double));
		mm->vycc[b]=malloc(cd3*sizeof(double));
		mm->ucartxi[b]=malloc(xd3*sizeof(double));
		mm->ucarteta[b]=malloc(xd3*sizeof(double));
		mm->vcartxi[b]=malloc(xd3*sizeof(double));
		mm->vcarteta[b]=malloc(xd3*sizeof(double));
		mm->uxxx[b]=malloc(xd3*sizeof(double));
		mm->vxxx[b]=malloc(xd3*sizeof(double));
		mm->uyxx[b]=malloc(xd3*sizeof(double));
		mm->vyxx[b]=malloc(xd3*sizeof(double));
		mm->ucartcc[b]=malloc(cd3*sizeof(double));
		mm->vcartcc[b]=malloc(cd3*sizeof(double));
		mm->ucartxx[b]=malloc(xd3*sizeof(double));
		mm->vcartxx[b]=malloc(xd3*sizeof(double));
		mm->ucartvv[b]=malloc(vd3*sizeof(double));
		mm->vcartvv[b]=malloc(vd3*sizeof(double));
		mm->ucarthh[b]=malloc(hd3*sizeof(double));
		mm->vcarthh[b]=malloc(hd3*sizeof(double));
		mm->term14_1[b]=malloc(cd3*sizeof(double));
		mm->term14_2[b]=malloc(cd3*sizeof(double));
		mm->term15_1[b]=malloc(xd3*sizeof(double));
		mm->term15_2[b]=malloc(xd3*sizeof(double));
		mm->term24_1[b]=malloc(xd3*sizeof(double));
		mm->term24_2[b]=malloc(xd3*sizeof(double));
		mm->term25_1[b]=malloc(cd3*sizeof(double));
		mm->term25_2[b]=malloc(cd3*sizeof(double));
		mm->I11_1[b]=malloc(vd3*sizeof(double));
		mm->I11_2[b]=malloc(vd3*sizeof(double));
		mm->I12_1[b]=malloc(vd3*sizeof(double));
		mm->I12_2[b]=malloc(vd3*sizeof(double));
		mm->I14_1[b]=malloc(vd3*sizeof(double));
		mm->I14_2[b]=malloc(vd3*sizeof(double));
		mm->I15_1[b]=malloc(vd3*sizeof(double));
		mm->I15_2[b]=malloc(vd3*sizeof(double));
		mm->I21_1[b]=malloc(hd3*sizeof(double));
		mm->I21_2[b]=malloc(hd3*sizeof(double));
		mm->I22_1[b]=malloc(hd3*sizeof(double));
		mm->I22_2[b]=malloc(hd3*sizeof(double));
		mm->I24_1[b]=malloc(hd3*sizeof(double));
		mm->I24_2[b]=malloc(hd3*sizeof(double));
		mm->I25_1[b]=malloc(hd3*sizeof(double));
		mm->I25_2[b]=malloc(hd3*sizeof(double));
		mm->nutxx[b]=malloc(xd3*sizeof(double));
	}
}
void freeState(State*s,Grid*g){
	int b;
	for(b=0;b<g->totalblocks;b++){
		//free(s->u[b]);
		//free(s->v[b]);
		free(s->sanu[b]);
		free(s->ut[b]);
		free(s->vt[b]);
		free(s->sanut[b]);
	}
	free(s->uvec);
	free(s->u);
	free(s->v);
	free(s->sanu);
	free(s->ut);
	free(s->vt);
	free(s->sanut);
}
void freeMomentum(Momentum*mm,Grid*g){
	int b;
	for(b=0;b<g->totalblocks;b++){
		free(mm->nut[b]);
		free(mm->ucc[b]);
		free(mm->vcc[b]);
		free(mm->uxx[b]);
		free(mm->vxx[b]);
		free(mm->uhh[b]);
		free(mm->vvv[b]);
		free(mm->uxcc[b]);
		free(mm->uycc[b]);
		free(mm->vxcc[b]);
		free(mm->vycc[b]);
		free(mm->ucartxi[b]);
		free(mm->ucarteta[b]);
		free(mm->vcartxi[b]);
		free(mm->vcarteta[b]);
		free(mm->uxxx[b]);
		free(mm->vxxx[b]);
		free(mm->uyxx[b]);
		free(mm->vyxx[b]);
		free(mm->ucartcc[b]);
		free(mm->vcartcc[b]);
		free(mm->ucartxx[b]);
		free(mm->vcartxx[b]);
		free(mm->ucartvv[b]);
		free(mm->vcartvv[b]);
		free(mm->ucarthh[b]);
		free(mm->vcarthh[b]);
		free(mm->term14_1[b]);
		free(mm->term14_2[b]);
		free(mm->term15_1[b]);
		free(mm->term15_2[b]);
		free(mm->term24_1[b]);
		free(mm->term24_2[b]);
		free(mm->term25_1[b]);
		free(mm->term25_2[b]);
		free(mm->I11_1[b]);
		free(mm->I11_2[b]);
		free(mm->I12_1[b]);
		free(mm->I12_2[b]);
		free(mm->I14_1[b]);
		free(mm->I14_2[b]);
		free(mm->I15_1[b]);
		free(mm->I15_2[b]);
		free(mm->I21_1[b]);
		free(mm->I21_2[b]);
		free(mm->I22_1[b]);
		free(mm->I22_2[b]);
		free(mm->I24_1[b]);
		free(mm->I24_2[b]);
		free(mm->I25_1[b]);
		free(mm->I25_2[b]);
		free(mm->nutxx[b]);
	}
	free(mm->nut);
	free(mm->ucc);
	free(mm->vcc);
	free(mm->uxx);
	free(mm->vxx);
	free(mm->uhh);
	free(mm->vvv);
	free(mm->uxcc);
	free(mm->uycc);
	free(mm->vxcc);
	free(mm->vycc);
	free(mm->ucartxi);
	free(mm->ucarteta);
	free(mm->vcartxi);
	free(mm->vcarteta);
	free(mm->uxxx);
	free(mm->vxxx);
	free(mm->uyxx);
	free(mm->vyxx);
	free(mm->ucartcc);
	free(mm->vcartcc);
	free(mm->ucartxx);
	free(mm->vcartxx);
	free(mm->ucartvv);
	free(mm->vcartvv);
	free(mm->ucarthh);
	free(mm->vcarthh);
	free(mm->term14_1);
	free(mm->term14_2);
	free(mm->term15_1);
	free(mm->term15_2);
	free(mm->term24_1);
	free(mm->term24_2);
	free(mm->term25_1);
	free(mm->term25_2);
	free(mm->I11_1);
	free(mm->I11_2);
	free(mm->I12_1);
	free(mm->I12_2);
	free(mm->I14_1);
	free(mm->I14_2);
	free(mm->I15_1);
	free(mm->I15_2);
	free(mm->I21_1);
	free(mm->I21_2);
	free(mm->I22_1);
	free(mm->I22_2);
	free(mm->I24_1);
	free(mm->I24_2);
	free(mm->I25_1);
	free(mm->I25_2);
	free(mm->nutxx);
	fclose(mm->out);
}
void mallocTransformation(Grid* g){
	int b;
	int vd3,hd3,cd3,xd3;
	g->a11=malloc(g->totalblocks*sizeof(double *));
	g->a12=malloc(g->totalblocks*sizeof(double *));
	g->a21=malloc(g->totalblocks*sizeof(double *));
	g->a22=malloc(g->totalblocks*sizeof(double *));
	g->a11cc=malloc(g->totalblocks*sizeof(double *));
	g->a12cc=malloc(g->totalblocks*sizeof(double *));
	g->a21cc=malloc(g->totalblocks*sizeof(double *));
	g->a22cc=malloc(g->totalblocks*sizeof(double *));
	g->a11xx=malloc(g->totalblocks*sizeof(double *));
	g->a12xx=malloc(g->totalblocks*sizeof(double *));
	g->a21xx=malloc(g->totalblocks*sizeof(double *));
	g->a22xx=malloc(g->totalblocks*sizeof(double *));
	g->a11vv=malloc(g->totalblocks*sizeof(double *));
	g->a12vv=malloc(g->totalblocks*sizeof(double *));
	g->a21hh=malloc(g->totalblocks*sizeof(double *));
	g->a22hh=malloc(g->totalblocks*sizeof(double *));
	g->Jcc=malloc(g->totalblocks*sizeof(double *));
	g->Jxx=malloc(g->totalblocks*sizeof(double *));
	g->Jvv=malloc(g->totalblocks*sizeof(double *));
	g->Jhh=malloc(g->totalblocks*sizeof(double *));
	g->vareas=malloc(g->totalblocks*sizeof(double *));
	g->hareas=malloc(g->totalblocks*sizeof(double *));
	g->covc11=malloc(g->totalblocks*sizeof(double *));
	g->covc12=malloc(g->totalblocks*sizeof(double *));
	g->covc21=malloc(g->totalblocks*sizeof(double *));
	g->covc22=malloc(g->totalblocks*sizeof(double *));
	g->contrac11=malloc(g->totalblocks*sizeof(double *));
	g->contrac12=malloc(g->totalblocks*sizeof(double *));
	g->contrac21=malloc(g->totalblocks*sizeof(double *));
	g->contrac22=malloc(g->totalblocks*sizeof(double *));
	g->xixxx=malloc(g->totalblocks*sizeof(double *));
	g->xiyxx=malloc(g->totalblocks*sizeof(double *));
	g->etaxxx=malloc(g->totalblocks*sizeof(double *));
	g->etayxx=malloc(g->totalblocks*sizeof(double *));
	g->xixcc=malloc(g->totalblocks*sizeof(double *));
	g->xiycc=malloc(g->totalblocks*sizeof(double *));
	g->etaxcc=malloc(g->totalblocks*sizeof(double *));
	g->etaycc=malloc(g->totalblocks*sizeof(double *));
	g->B11=malloc(g->totalblocks*sizeof(double *));
	g->B22=malloc(g->totalblocks*sizeof(double *));
	g->B12=malloc(g->totalblocks*sizeof(double *));
	for(b=0;b<g->totalblocks;b++){
		vd3=g->vdim[b][3];
		hd3=g->hdim[b][3];
		cd3=g->cdim[b][3];
		xd3=g->xdim[b][3];
		g->a11[b]=malloc(hd3*sizeof(double));
		g->a12[b]=malloc(hd3*sizeof(double));
		g->a21[b]=malloc(vd3*sizeof(double));
		g->a22[b]=malloc(vd3*sizeof(double));
		g->a11cc[b]=malloc(cd3*sizeof(double));
		g->a12cc[b]=malloc(cd3*sizeof(double));
		g->a21cc[b]=malloc(cd3*sizeof(double));
		g->a22cc[b]=malloc(cd3*sizeof(double));
		g->a11xx[b]=malloc(xd3*sizeof(double));
		g->a12xx[b]=malloc(xd3*sizeof(double));
		g->a21xx[b]=malloc(xd3*sizeof(double));
		g->a22xx[b]=malloc(xd3*sizeof(double));
		g->a11vv[b]=malloc(vd3*sizeof(double));
		g->a12vv[b]=malloc(vd3*sizeof(double));
		g->a21hh[b]=malloc(hd3*sizeof(double));
		g->a22hh[b]=malloc(hd3*sizeof(double));
		g->Jcc[b]=malloc(cd3*sizeof(double));
		g->Jxx[b]=malloc(xd3*sizeof(double));
		g->Jvv[b]=malloc(vd3*sizeof(double));
		g->Jhh[b]=malloc(hd3*sizeof(double));
		g->vareas[b]=malloc(vd3*sizeof(double));
		g->hareas[b]=malloc(hd3*sizeof(double));
		g->covc11[b]=malloc(xd3*sizeof(double));
		g->covc12[b]=malloc(xd3*sizeof(double));
		g->covc21[b]=malloc(xd3*sizeof(double));
		g->covc22[b]=malloc(xd3*sizeof(double));
		g->contrac11[b]=malloc(xd3*sizeof(double));
		g->contrac12[b]=malloc(xd3*sizeof(double));
		g->contrac21[b]=malloc(xd3*sizeof(double));
		g->contrac22[b]=malloc(xd3*sizeof(double));
		g->xixxx[b]=malloc(xd3*sizeof(double));
		g->xiyxx[b]=malloc(xd3*sizeof(double));
		g->etaxxx[b]=malloc(xd3*sizeof(double));
		g->etayxx[b]=malloc(xd3*sizeof(double));
		g->xixcc[b]=malloc(cd3*sizeof(double));
		g->xiycc[b]=malloc(cd3*sizeof(double));
		g->etaxcc[b]=malloc(cd3*sizeof(double));
		g->etaycc[b]=malloc(cd3*sizeof(double));
		g->B11[b]=malloc(vd3*sizeof(double));
		g->B22[b]=malloc(hd3*sizeof(double));
		g->B12[b]=malloc(xd3*sizeof(double));
	}
}
void freeTransformation(Grid*g){
	int b;
	for(b=0;b<g->totalblocks;b++){
		free(g->a11[b]);
		free(g->a12[b]);
		free(g->a21[b]);
		free(g->a22[b]);
		free(g->a11cc[b]);
		free(g->a12cc[b]);
		free(g->a21cc[b]);
		free(g->a22cc[b]);
		free(g->a11xx[b]);
		free(g->a12xx[b]);
		free(g->a21xx[b]);
		free(g->a22xx[b]);
		free(g->a11vv[b]);
		free(g->a12vv[b]);
		free(g->a21hh[b]);
		free(g->a22hh[b]);
		free(g->Jcc[b]);
		free(g->Jxx[b]);
		free(g->Jvv[b]);
		free(g->Jhh[b]);
		free(g->vareas[b]);
		free(g->hareas[b]);
		free(g->covc11[b]);
		free(g->covc12[b]);
		free(g->covc21[b]);
		free(g->covc22[b]);
		free(g->contrac11[b]);
		free(g->contrac12[b]);
		free(g->contrac21[b]);
		free(g->contrac22[b]);
		free(g->xixxx[b]);
		free(g->xiyxx[b]);
		free(g->etaxxx[b]);
		free(g->etayxx[b]);
		free(g->xixcc[b]);
		free(g->xiycc[b]);
		free(g->etaxcc[b]);
		free(g->etaycc[b]);
		free(g->B11[b]);
		free(g->B22[b]);
		free(g->B12[b]);
	}
	free(g->a11);
	free(g->a12);
	free(g->a21);
	free(g->a22);
	free(g->a11cc);
	free(g->a12cc);
	free(g->a21cc);
	free(g->a22cc);
	free(g->a11xx);
	free(g->a12xx);
	free(g->a21xx);
	free(g->a22xx);
	free(g->a11vv);
	free(g->a12vv);
	free(g->a21hh);
	free(g->a22hh);
	free(g->Jcc);
	free(g->Jxx);
	free(g->Jvv);
	free(g->Jhh);
	free(g->vareas);
	free(g->hareas);
	free(g->covc11);
	free(g->covc12);
	free(g->covc21);
	free(g->covc22);
	free(g->contrac11);
	free(g->contrac12);
	free(g->contrac21);
	free(g->contrac22);
	free(g->xixxx);
	free(g->xiyxx);
	free(g->etaxxx);
	free(g->etayxx);
	free(g->xixcc);
	free(g->xiycc);
	free(g->etaxcc);
	free(g->etaycc);
	free(g->B11);
	free(g->B22);
	free(g->B12);
}
void mallocSymlap(SymLap*s,Grid*g,Interfaces*is){
	int tb=g->totalblocks;
	int tc=g->totalcells;
	int b,cd3,xd3;
	//
	s->dim=tc;
	s->b=malloc(sizeof(double)*tc);
	s->x=malloc(sizeof(double)*tc);
	s->d0=malloc(sizeof(double *)*tb);
	s->d1=malloc(sizeof(double *)*tb);
	s->d2=malloc(sizeof(double *)*tb);
	s->d3=malloc(sizeof(double *)*tb);
	s->d4=malloc(sizeof(double *)*tb);
	s->cc=malloc(sizeof(double *)*tb);
	s->xx=malloc(sizeof(double *)*tb);
	for(b=0;b<tb;b++){
		cd3=g->cdim[b][3];
		xd3=g->xdim[b][3];
		s->d0[b]=malloc(sizeof(double)*cd3);
		s->d1[b]=malloc(sizeof(double)*cd3);//malloc(sizeof(double)*(cd3-1));
		s->d2[b]=malloc(sizeof(double)*cd3);//malloc(sizeof(double)*(cd3-cd0+1));
		s->d3[b]=malloc(sizeof(double)*cd3);//malloc(sizeof(double)*(cd3-cd0));
		s->d4[b]=malloc(sizeof(double)*cd3);//malloc(sizeof(double)*(cd3-cd0-1));
		s->cc[b]=malloc(sizeof(double)*cd3);
		s->xx[b]=malloc(sizeof(double)*xd3);
	}
	//csr
	//s->csr.nr=tc;
	s->csr.nn=is->totalinterfacepoints*3;//(is->totalinterfacepoints-is->totalinterfaces)*3-2*is->totalinterfaces;
	s->csr.v=malloc(sizeof(double)*s->csr.nn);
	s->csr.ci=malloc(sizeof(int)*s->csr.nn);
	s->csr.ri=malloc(sizeof(int)*s->csr.nn);
}
void mallocCSR2(SymLap*s,int diagcorners){
	s->csr2.nn=diagcorners;
	s->csr2.ri=malloc(sizeof(int)*s->csr2.nn);
	s->csr2.v=malloc(sizeof(double)*s->csr2.nn);
	s->csr2.ci=malloc(sizeof(int)*s->csr2.nn);
}
void freeSymlap(SymLap*s,Grid*g){
	int b;
	for(b=0;b<g->totalblocks;b++){
		free(s->d0[b]);
		free(s->d1[b]);
		free(s->d2[b]);
		free(s->d3[b]);
		free(s->d4[b]);
		free(s->cc[b]);
		free(s->xx[b]);
	}
	free(s->b);
	free(s->x);
	free(s->d0);
	free(s->d1);
	free(s->d2);
	free(s->d3);
	free(s->d4);
	free(s->cc);
	free(s->xx);
	//csr
	free(s->csr.ci);
	free(s->csr.ri);
	free(s->csr.v);
	free(s->csr2.ci);
	free(s->csr2.ri);
	free(s->csr2.v);
}
void prepls(Grid* g,Interfaces* is,LinearSystem* ls){
	ls->crows=g->totaledges+is->totalinterfacepoints;
	ls->ccols=g->totalnodes;
	ls->ctcmain	=malloc(ls->ccols*sizeof(double));
	ls->s		=malloc(ls->ccols*sizeof(double));
	ls->ctrhs	=malloc(ls->ccols*sizeof(double));
	ls->cs		=malloc(ls->crows*sizeof(double));
	ls->rhs		=malloc(ls->crows*sizeof(double));
	ls->bp		=malloc(ls->ccols*sizeof(double));
	ls->r		=malloc(ls->ccols*sizeof(double));
	ls->d		=malloc(ls->ccols*sizeof(double));
	ls->q		=malloc(ls->ccols*sizeof(double));
}
void freels(LinearSystem*ls){
	free(ls->ctcmain);
	free(ls->s);
	free(ls->ctrhs);
	free(ls->cs);
	free(ls->rhs);
	free(ls->bp);
	free(ls->r);
	free(ls->d);
	free(ls->q);
}
void mallocSA(SpalartAllmaras*sa,Grid*g){
	int b,cd3;
	sa->bcid=malloc(sizeof(char*)*g->totalblocks);
	sa->fixed=malloc(sizeof(double)*g->totalblocks);
	sa->sanuxi=malloc(sizeof(double*)*g->totalblocks);
	sa->sanueta=malloc(sizeof(double*)*g->totalblocks);
	sa->sanux=malloc(sizeof(double*)*g->totalblocks);
	sa->sanuy=malloc(sizeof(double*)*g->totalblocks);
	sa->dissTermx=malloc(sizeof(double*)*g->totalblocks);
	sa->dissTermy=malloc(sizeof(double*)*g->totalblocks);
	sa->dissTermxxi=malloc(sizeof(double*)*g->totalblocks);
	sa->dissTermyxi=malloc(sizeof(double*)*g->totalblocks);
	sa->dissTermxeta=malloc(sizeof(double*)*g->totalblocks);
	sa->dissTermyeta=malloc(sizeof(double*)*g->totalblocks);
	sa->dissTermxx=malloc(sizeof(double*)*g->totalblocks);
	sa->dissTermyy=malloc(sizeof(double*)*g->totalblocks);
	sa->advection=malloc(sizeof(double*)*g->totalblocks);
	sa->production=malloc(sizeof(double*)*g->totalblocks);
	sa->wallDestruction=malloc(sizeof(double*)*g->totalblocks);
	sa->dissipation=malloc(sizeof(double*)*g->totalblocks);
	sa->chi=malloc(sizeof(double*)*g->totalblocks);
	sa->fv1=malloc(sizeof(double*)*g->totalblocks);
	sa->fv2=malloc(sizeof(double*)*g->totalblocks);
	sa->ft2=malloc(sizeof(double*)*g->totalblocks);
	sa->S2M=malloc(sizeof(double*)*g->totalblocks);
	sa->SM=malloc(sizeof(double*)*g->totalblocks);
	sa->S2=malloc(sizeof(double*)*g->totalblocks);
	sa->S=malloc(sizeof(double*)*g->totalblocks);
	sa->r=malloc(sizeof(double*)*g->totalblocks);
	sa->g=malloc(sizeof(double*)*g->totalblocks);
	sa->fw=malloc(sizeof(double*)*g->totalblocks);
	sa->fn=malloc(sizeof(double*)*g->totalblocks);
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		sa->sanuxi[b]=malloc(sizeof(double)*cd3);
		sa->sanueta[b]=malloc(sizeof(double)*cd3);
		sa->sanux[b]=malloc(sizeof(double)*cd3);
		sa->sanuy[b]=malloc(sizeof(double)*cd3);
		sa->dissTermx[b]=malloc(sizeof(double)*cd3);
		sa->dissTermy[b]=malloc(sizeof(double)*cd3);
		sa->dissTermxxi[b]=malloc(sizeof(double)*cd3);
		sa->dissTermyxi[b]=malloc(sizeof(double)*cd3);
		sa->dissTermxeta[b]=malloc(sizeof(double)*cd3);
		sa->dissTermyeta[b]=malloc(sizeof(double)*cd3);
		sa->dissTermxx[b]=malloc(sizeof(double)*cd3);
		sa->dissTermyy[b]=malloc(sizeof(double)*cd3);
		sa->advection[b]=malloc(sizeof(double)*cd3);
		sa->production[b]=malloc(sizeof(double)*cd3);
		sa->wallDestruction[b]=malloc(sizeof(double)*cd3);
		sa->dissipation[b]=malloc(sizeof(double)*cd3);
		sa->chi[b]=malloc(sizeof(double)*cd3);
		sa->fv1[b]=malloc(sizeof(double)*cd3);
		sa->fv2[b]=malloc(sizeof(double)*cd3);
		sa->ft2[b]=malloc(sizeof(double)*cd3);
		sa->S2M[b]=malloc(sizeof(double)*cd3);
		sa->SM[b]=malloc(sizeof(double)*cd3);
		sa->S2[b]=malloc(sizeof(double)*cd3);
		sa->S[b]=malloc(sizeof(double)*cd3);
		sa->r[b]=malloc(sizeof(double)*cd3);
		sa->g[b]=malloc(sizeof(double)*cd3);
		sa->fw[b]=malloc(sizeof(double)*cd3);
		sa->fn[b]=malloc(sizeof(double)*cd3);
		sa->fixed[b]=sa->fixedSanu;
		sa->bcid[b]=malloc(sizeof(char)*4);
	}
}
void freeSA(SpalartAllmaras*sa,Grid*g){
	int b;
	for(b=0;b<g->totalblocks;b++){
		free(sa->bcid[b]);
		free(sa->sanuxi[b]);
		free(sa->sanueta[b]);
		free(sa->sanux[b]);
		free(sa->sanuy[b]);
		free(sa->dissTermx[b]);
		free(sa->dissTermy[b]);
		free(sa->dissTermxxi[b]);
		free(sa->dissTermyxi[b]);
		free(sa->dissTermxeta[b]);
		free(sa->dissTermyeta[b]);
		free(sa->dissTermxx[b]);
		free(sa->dissTermyy[b]);
		free(sa->advection[b]);
		free(sa->production[b]);
		free(sa->wallDestruction[b]);
		free(sa->dissipation[b]);
		free(sa->chi[b]);
		free(sa->fv1[b]);
		free(sa->fv2[b]);
		free(sa->ft2[b]);
		free(sa->S2M[b]);
		free(sa->SM[b]);
		free(sa->S2[b]);
		free(sa->S[b]);
		free(sa->r[b]);
		free(sa->g[b]);
		free(sa->fw[b]);
		free(sa->fn[b]);
	}
	free(sa->bcid);
	free(sa->fixed);
	free(sa->sanuxi);
	free(sa->sanueta);
	free(sa->sanux);
	free(sa->sanuy);
	free(sa->dissTermx);
	free(sa->dissTermy);
	free(sa->dissTermxxi);
	free(sa->dissTermyxi);
	free(sa->dissTermxeta);
	free(sa->dissTermyeta);
	free(sa->dissTermxx);
	free(sa->dissTermyy);
	free(sa->advection);
	free(sa->production);
	free(sa->wallDestruction);
	free(sa->dissipation);
	free(sa->chi);
	free(sa->fv1);
	free(sa->fv2);
	free(sa->ft2);
	free(sa->S2M);
	free(sa->SM);
	free(sa->S2);
	free(sa->S);
	free(sa->r);
	free(sa->g);
	free(sa->fw);
	free(sa->fn);
}
void mallocWd(WallDistances*wd,Grid*g){
	int b,tb=g->totalblocks;
	wd->bcid=malloc(sizeof(char*)*tb);
	wd->fixed=malloc(tb*sizeof(double));
	wd->phi=malloc(tb*sizeof(double*));
	wd->phixi=malloc(tb*sizeof(double*));
	wd->phieta=malloc(tb*sizeof(double*));
	wd->cc=malloc(tb*sizeof(double*));
	wd->xx=malloc(tb*sizeof(double*));
	for(b=0;b<tb;b++){
		wd->phi[b]=malloc(g->cdim[b][3]*sizeof(double));
		wd->phixi[b]=malloc(g->cdim[b][3]*sizeof(double));
		wd->phieta[b]=malloc(g->cdim[b][3]*sizeof(double));
		wd->cc[b]=malloc(g->cdim[b][3]*sizeof(double));
		wd->xx[b]=malloc(g->xdim[b][3]*sizeof(double));
		wd->fixed[b]=0;
		wd->bcid[b]=malloc(sizeof(char)*4);
	}
}
void freeWd(WallDistances*wd,Grid*g){
	int b;
	for(b=0;b<g->totalblocks;b++){
		free(wd->bcid[b]);
		free(wd->phi[b]);
		free(wd->phixi[b]);
		free(wd->phieta[b]);
		free(wd->cc[b]);
		free(wd->xx[b]);
	}
	free(wd->fixed);
	free(wd->bcid);
	free(wd->phi);
	free(wd->phixi);
	free(wd->phieta);
	free(wd->cc);
	free(wd->xx);
}
void mallocInterfaces(Interfaces*is,int ti){
	is->totalinterfaces=ti;
	is->block=malloc(sizeof(int)*ti);
	is->side=malloc(sizeof(int)*ti);
	is->adjacentBlock=malloc(sizeof(int)*ti);
	is->xsize=malloc(sizeof(int)*ti);
	is->rowstart=malloc(sizeof(int)*ti);
	is->blockx=malloc(sizeof(int *)*ti);
	is->adjacentBlockx=malloc(sizeof(int *)*ti);
	is->e=malloc(sizeof(int *)*ti);
	is->ae=malloc(sizeof(int *)*ti);
	is->out=fopen("out/interfaces.out","w");
}
void freeInterfaces(Interfaces*is,Grid*g){
	int i,b;
	for(i=0;i<is->totalinterfaces;i++){
		free(is->blockx[i]);
		free(is->adjacentBlockx[i]);
		free(is->e[i]);
		free(is->ae[i]);
	}
	for(i=0;i<is->totalicorners;i++){
		free(is->icorners[i]);
	}
	for(b=0;b<g->totalblocks;b++){
		free(is->conn[b]);
	}
	free(is->block);
	free(is->side);
	free(is->adjacentBlock);
	free(is->xsize);
	free(is->rowstart);
	free(is->blockx);
	free(is->adjacentBlockx);
	free(is->e);
	free(is->ae);
	if(is->totalicorners>0){
		free(is->icorners);
	}
	free(is->conn);
	fclose(is->out);
}
void freeBc(BoundaryConditions*bc,Grid*g){
	int b,s;
	for(b=0;b<g->totalblocks;b++){
		free(bc->id[b]);
		free(bc->fv[b]);
		for(s=0;s<4;s++){
			free(bc->side[b][s]);
		}
		free(bc->side[b]);
	}
	free(bc->side);
	free(bc->id);
	free(bc->fv);
}
void mallocRK(Grid*g,RungeKutta*rk){
	int b,s,cd3,vd3,hd3,tb;
	rk->f=malloc(sizeof(State)*rk->stages);
	tb=g->totalblocks;
	for(s=0;s<rk->stages;s++){
		rk->f[s].u=malloc(sizeof(double*)*tb);
		rk->f[s].v=malloc(sizeof(double*)*tb);
		rk->f[s].sanu=malloc(sizeof(double*)*tb);
		rk->f[s].ut=malloc(sizeof(double*)*tb);
		rk->f[s].vt=malloc(sizeof(double*)*tb);
		rk->f[s].sanut=malloc(sizeof(double*)*tb);
		for(b=0;b<tb;b++){
			cd3=g->cdim[b][3];
			vd3=g->vdim[b][3];
			hd3=g->hdim[b][3];
			rk->f[s].u[b]=malloc(sizeof(double)*vd3);
			rk->f[s].v[b]=malloc(sizeof(double)*hd3);
			rk->f[s].sanu[b]=malloc(sizeof(double)*cd3);
			rk->f[s].ut[b]=malloc(sizeof(double)*vd3);
			rk->f[s].vt[b]=malloc(sizeof(double)*hd3);
			rk->f[s].sanut[b]=malloc(sizeof(double)*cd3);
		}
	}
}
void freeRK(Grid*g,RungeKutta*rk){
	int b,s;
	for(s=0;s<rk->stages;s++){
		for(b=0;b<g->totalblocks;b++){
			free(rk->f[s].u[b]);
			free(rk->f[s].v[b]);
			free(rk->f[s].sanu[b]);
			free(rk->f[s].ut[b]);
			free(rk->f[s].vt[b]);
			free(rk->f[s].sanut[b]);
		}
		free(rk->f[s].u);
		free(rk->f[s].v);
		free(rk->f[s].sanu);
		free(rk->f[s].ut);
		free(rk->f[s].vt);
		free(rk->f[s].sanut);
	}
	free(rk->f);
}




void freeAll(State*st,Numerics*n,Switches*sw,Grid*g,Initial*in,
		Interfaces*is,LinearSystem*ls,Momentum*mm,BoundaryConditions*bc,
		SymLap*sl,SpalartAllmaras*sa,WallDistances*wd,RungeKutta*rk){

	freels(ls);

	freeState(st,g);

	freeTransformation(g);

	freeSymlap(sl,g);

	freeSA(sa,g);

	freeWd(wd,g);

	freeInterfaces(is,g);

	freeBc(bc,g);

	freeMomentum(mm,g);

	freeRK(g,rk);

	freeGrid(g);

}
