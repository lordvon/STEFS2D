void invariantTimeDerivative(Grid* g,Momentum* mm,State* s){
	int b,vi,hi,vd3,hd3;
	double umom[2],vmom[2],a2[2],a1[2];
	//Explicit Euler (first-order) time derivatives.
	//Finite-volume formulation.
	for(b=0;b<g->totalblocks;b++){
		vd3=g->vdim[b][3];
		hd3=g->hdim[b][3];
		//vertical edges
		for(vi=0;vi<vd3;vi++){
			a2[0]=g->a21[b][vi];
			a2[1]=g->a22[b][vi];
			umom[0]=-mm->I11_1[b][vi]-mm->I12_1[b][vi]+
					mm->I14_1[b][vi]+mm->I15_1[b][vi];
			umom[1]=-mm->I11_2[b][vi]-mm->I12_2[b][vi]+
					mm->I14_2[b][vi]+mm->I15_2[b][vi];
			s->ut[b][vi]=-cross(a2,umom)/g->vareas[b][vi];
		}
		//horizontal edges
		for(hi=0;hi<hd3;hi++){
			a1[0]=g->a11[b][hi];
			a1[1]=g->a12[b][hi];
			vmom[0]=-mm->I21_1[b][hi]-mm->I22_1[b][hi]+
					mm->I24_1[b][hi]+mm->I25_1[b][hi];
			vmom[1]=-mm->I21_2[b][hi]-mm->I22_2[b][hi]+
					mm->I24_2[b][hi]+mm->I25_2[b][hi];
			s->vt[b][hi]=cross(a1,vmom)/g->hareas[b][hi];
		}
	}
}
void fillTimeDerivative(State*st,
		Grid*g,BoundaryConditions*bc,Interfaces*is,WallDistances*wd,
		Momentum*mm,SpalartAllmaras*sa,
		Property*p,Switches*sw,Numerics*n){
	//Calculates ut,vt,sanut based on the u,v,sanu in State*st.
	interpolateState(st,mm,g,is,bc);
	cartesianConvert(st,g,mm);
	cartesianDerivatives(g,mm,is);
	convection(g,mm,bc,is);
	if(sw->turbmod>0){
		computeSanut(st,sa,p,g,is,mm,wd->cc);
		updateNut(st,sa,g,p,mm,bc,is);
	}
	viscosity(p,g,mm,bc,is);
	invariantTimeDerivative(g,mm,st);
}
void fillRhs(Grid*g,Interfaces*is,State*s,Numerics*n,LinearSystem*ls){
	//Adds time derivative term to rhs.
	int b,vi,hi,vd3,hd3,ip,es;
	//Assign to RHS.
	for(b=0;b<g->totalblocks;b++){
		vd3=g->vdim[b][3];
		hd3=g->hdim[b][3];
		es=g->estart[b];
		for(vi=0;vi<vd3;vi++){
			ls->rhs[es+vi]=s->u[b][vi]+n->dt*s->ut[b][vi];
		}
		for(hi=0;hi<hd3;hi++){
			ls->rhs[es+vd3+hi]=s->v[b][hi]+n->dt*s->vt[b][hi];
		}
	}
	//printf("Last edge counted: %d\nTotal edges: %d\n",es+vd3+hi,g->totaledges);
	for(ip=0;ip<is->totalinterfacepoints;ip++){
		ls->rhs[g->totaledges+ip]=0;
	}
	//printf("Last interface row: %d\nTotal rows: %d\n",g->totaledges+ip,ls->crows);
}
