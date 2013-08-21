void fillDerivativecc(Grid* g,Momentum* mm){
	int b,ci,vi,hi;
	int cd3,cd0,hd0;
	double ucarteta,vcarteta,ucartxi,vcartxi;
	double Jinv;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		cd0=g->cdim[b][0];
		hd0=g->hdim[b][0];
		for(ci=0;ci<cd3;ci++){
			vi=ci+ci/cd0;
			hi=ci;
			ucarteta=mm->ucarthh[b][hi+hd0]-mm->ucarthh[b][hi];
			vcarteta=mm->vcarthh[b][hi+hd0]-mm->vcarthh[b][hi];
			ucartxi=mm->ucartvv[b][vi+1]-mm->ucartvv[b][vi];
			vcartxi=mm->vcartvv[b][vi+1]-mm->vcartvv[b][vi];
			Jinv=1/g->Jcc[b][ci];
			mm->uxcc[b][ci]=Jinv*(
					g->a22cc[b][ci]*ucartxi-
					g->a12cc[b][ci]*ucarteta);
			mm->vxcc[b][ci]=Jinv*(
					g->a22cc[b][ci]*vcartxi-
					g->a12cc[b][ci]*vcarteta);
			mm->uycc[b][ci]=Jinv*(
					-g->a21cc[b][ci]*ucartxi+
					g->a11cc[b][ci]*ucarteta);
			mm->vycc[b][ci]=Jinv*(
					-g->a21cc[b][ci]*vcartxi+
					g->a11cc[b][ci]*vcarteta);
		}
	}
}
void fillDerivativexx(Grid* g,Momentum* mm,Interfaces*is){
	//STILL NEED TO ENFORCE WALL CORNERS!
	fillCartXixx(g,mm);
	fillCartEtaxx(g,mm);
	fillCartXixxInterface(g,mm,is);
	fillCartEtaxxInterface(g,mm,is);
	int b,xi;
	for(b=0;b<g->totalblocks;b++){
		for(xi=0;xi<g->xdim[b][3];xi++){
			mm->uxxx[b][xi]=mm->ucartxi[b][xi]*g->contrac11[b][xi]+
					mm->ucarteta[b][xi]*g->contrac21[b][xi];
		}
	}
	for(b=0;b<g->totalblocks;b++){
		for(xi=0;xi<g->xdim[b][3];xi++){
			mm->uyxx[b][xi]=mm->ucartxi[b][xi]*g->contrac12[b][xi]+
					mm->ucarteta[b][xi]*g->contrac22[b][xi];
		}
	}
	for(b=0;b<g->totalblocks;b++){
		for(xi=0;xi<g->xdim[b][3];xi++){
			mm->vxxx[b][xi]=mm->vcartxi[b][xi]*g->contrac11[b][xi]+
					mm->vcarteta[b][xi]*g->contrac21[b][xi];
		}
	}
	for(b=0;b<g->totalblocks;b++){
		for(xi=0;xi<g->xdim[b][3];xi++){
			mm->vyxx[b][xi]=mm->vcartxi[b][xi]*g->contrac12[b][xi]+
					mm->vcarteta[b][xi]*g->contrac22[b][xi];
		}
	}
}
