void RK3(State*s0,RungeKutta*rk,
		Grid*g,BoundaryConditions*bc,Interfaces*is,WallDistances*wd,Momentum*mm,
		SpalartAllmaras*sa,Property*p,Switches*sw,Numerics*n){
	/* Takes u,v,sanu in s0 and fills the time derivatives using RK3.
	Butcher tableau:
	0		|	0
	0.5		|	0.5
	0.5		|	0.25	0.25
			|	0		-1		2
	*/
	fillTimeDerivative(s0,g,bc,is,wd,mm,sa,p,sw,n);
	int b,vi,hi,ci,vd3,cd3,hd3;
	double a1,a2;
	a1=0.5;a2=0.25;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		vd3=g->vdim[b][3];
		hd3=g->hdim[b][3];
		for(vi=0;vi<vd3;vi++){
			rk->f[0].u[b][vi]=s0->u[b][vi]+
					a1*n->dt*s0->ut[b][vi];
		}
		for(hi=0;hi<hd3;hi++){
			rk->f[0].v[b][hi]=s0->v[b][hi]+
					a1*n->dt*s0->vt[b][hi];
		}
		for(ci=0;ci<cd3;ci++){
			rk->f[0].sanu[b][ci]=s0->sanu[b][ci]+
					a1*n->dt*s0->sanut[b][ci];
		}
	}
	fillTimeDerivative(&rk->f[0],g,bc,is,wd,mm,sa,p,sw,n);
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		vd3=g->vdim[b][3];
		hd3=g->hdim[b][3];
		for(vi=0;vi<vd3;vi++){
			rk->f[1].u[b][vi]=s0->u[b][vi]+
					a2*n->dt*(s0->ut[b][vi]+
								rk->f[0].ut[b][vi]);
		}
		for(hi=0;hi<hd3;hi++){
			rk->f[1].v[b][hi]=s0->v[b][hi]+
					a2*n->dt*(s0->vt[b][hi]+
								rk->f[0].vt[b][hi]);
		}
		for(ci=0;ci<cd3;ci++){
			rk->f[1].sanu[b][ci]=s0->sanu[b][ci]+
					a2*n->dt*(s0->sanut[b][ci]+
								rk->f[0].sanut[b][ci]);
		}
	}
	fillTimeDerivative(&rk->f[1],g,bc,is,wd,mm,sa,p,sw,n);
	//Replace s0 derivatives with final values.
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		vd3=g->vdim[b][3];
		hd3=g->hdim[b][3];
		for(vi=0;vi<vd3;vi++){
			s0->ut[b][vi]=
					-rk->f[0].ut[b][vi]+
					2*rk->f[1].ut[b][vi];
		}
		for(hi=0;hi<hd3;hi++){
			s0->vt[b][hi]=
					-rk->f[0].vt[b][hi]+
					2*rk->f[1].vt[b][hi];
		}
		for(ci=0;ci<cd3;ci++){
			s0->sanut[b][ci]=
					-rk->f[0].sanut[b][ci]+
					2*rk->f[1].sanut[b][ci];
		}
	}
}
