void fillCartcc(Grid* g,Momentum* mm){
	int b,ci,cd3;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			mm->ucartcc[b][ci]=(g->a11cc[b][ci]*mm->ucc[b][ci]+
					g->a21cc[b][ci]*mm->vcc[b][ci])/g->Jcc[b][ci];
		}
		for(ci=0;ci<cd3;ci++){
			mm->vcartcc[b][ci]=(g->a12cc[b][ci]*mm->ucc[b][ci]+
					g->a22cc[b][ci]*mm->vcc[b][ci])/g->Jcc[b][ci];
		}
	}
}
void fillCartxx(Grid* g,Momentum* mm){
	//Assumes all BCs have been applied already to the invariant velocity field.
	int b,xi,xd3;
	for(b=0;b<g->totalblocks;b++){
		xd3=g->xdim[b][3];
		for(xi=0;xi<xd3;xi++){
			mm->ucartxx[b][xi]=(g->a11xx[b][xi]*mm->uxx[b][xi]+
					g->a21xx[b][xi]*mm->vxx[b][xi])/g->Jxx[b][xi];
		}
		for(xi=0;xi<xd3;xi++){
			mm->vcartxx[b][xi]=(g->a12xx[b][xi]*mm->uxx[b][xi]+
					g->a22xx[b][xi]*mm->vxx[b][xi])/g->Jxx[b][xi];
		}
	}
}
void fillCarthh(State* s,Grid* g,Momentum* mm){
	//Assumes all BCs have been applied already to the invariant velocity field.
	int b,hi,hd3;
	for(b=0;b<g->totalblocks;b++){
		hd3=g->hdim[b][3];
		for(hi=0;hi<hd3;hi++){
			mm->ucarthh[b][hi]=(g->a11[b][hi]*mm->uhh[b][hi]+
					g->a21hh[b][hi]*s->v[b][hi])/g->Jhh[b][hi];
		}
		for(hi=0;hi<hd3;hi++){
			mm->vcarthh[b][hi]=(g->a12[b][hi]*mm->uhh[b][hi]+
					g->a22hh[b][hi]*s->v[b][hi])/g->Jhh[b][hi];
		}
	}
}
void fillCartvv(State* s,Grid* g,Momentum* mm){
	//Assumes all BCs have been applied already to the invariant velocity field.
	int b,vi,vd3;
	for(b=0;b<g->totalblocks;b++){
		vd3=g->vdim[b][3];
		for(vi=0;vi<vd3;vi++){
			mm->ucartvv[b][vi]=(g->a11vv[b][vi]*s->u[b][vi]+
					g->a21[b][vi]*mm->vvv[b][vi])/g->Jvv[b][vi];
		}
		for(vi=0;vi<vd3;vi++){
			mm->vcartvv[b][vi]=(g->a12vv[b][vi]*s->u[b][vi]+
					g->a22[b][vi]*mm->vvv[b][vi])/g->Jvv[b][vi];
		}
	}
}
