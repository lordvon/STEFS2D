void invariantTimeDerivativeCentralDifference(Grid* g,Momentum* mm,State* s){
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
void invariantTimeDerivativeUpwindCC(Dimension*d,Grid*g,Momentum*mm,State*s,
		BoundaryConditions*bc){
	int b,i,j,vi,ci,hi,adj;
	double umom[2],vmom[2],a2[2],a1[2];
	//Explicit Euler (first-order) time derivatives.
	//Finite-volume formulation.
	//Interior
	for(b=0;b<d->tb;b++){
		//vertical edges
		for(j=0;j<d->v[b].j;j++){
			for(i=1;i<d->v[b].i-1;i++){
				vi=i+j*d->v[b].i;
				ci=vi-j;
				a2[0]=g->a21[b][vi];
				a2[1]=g->a22[b][vi];
				umom[0]=-mm->I11_1[b][vi]-mm->I12_1[b][vi]+
						mm->I14_1[b][vi]+mm->I15_1[b][vi];
				umom[1]=-mm->I11_2[b][vi]-mm->I12_2[b][vi]+
						mm->I14_2[b][vi]+mm->I15_2[b][vi];
				if(s->u[b][vi]>=0){ adj=-1; }
				else { adj=0; }
				s->ut[b][vi]=-cross(a2,umom)/g->Jcc[b][ci+adj];
			}
		}
		//horizontal edges
		for(j=1;j<d->h[b].j-1;j++){
			for(i=0;i<d->h[b].i;i++){
				hi=i+j*d->h[b].i;
				ci=hi;
				a1[0]=g->a11[b][hi];
				a1[1]=g->a12[b][hi];
				vmom[0]=-mm->I21_1[b][hi]-mm->I22_1[b][hi]+
						mm->I24_1[b][hi]+mm->I25_1[b][hi];
				vmom[1]=-mm->I21_2[b][hi]-mm->I22_2[b][hi]+
						mm->I24_2[b][hi]+mm->I25_2[b][hi];
				if(s->v[b][hi]>=0){ adj=-d->c[b].i; }
				else { adj=0; }
				s->vt[b][hi]=cross(a1,vmom)/g->Jcc[b][ci+adj];
			}
		}
	}
	//Boundaries
	int w,z,f,side;
	Boundary*bo;
	for(w=0;w<bc->wn;w++){
		side=bc->w[w].s;
		b=bc->w[w].b;
		if(side%2==1){
			bo=&d->vside[b].b[side];
			for(vi=bo->st;vi<=bo->en;vi+=bo->iv){
				s->ut[b][vi]=0;
			}
		} else {
			bo=&d->hside[b].b[side];
			for(hi=bo->st;hi<=bo->en;hi+=bo->iv){
				s->vt[b][hi]=0;
			}
		}
	}
	for(f=0;f<bc->fn;f++){
		side=bc->f[f].s;
		b=bc->f[f].b;
		if(side%2==1){
			bo=&d->vside[b].b[side];
			for(vi=bo->st;vi<=bo->en;vi+=bo->iv){
				s->ut[b][vi]=0;
			}
		} else {
			bo=&d->hside[b].b[side];
			for(hi=bo->st;hi<=bo->en;hi+=bo->iv){
				s->vt[b][hi]=0;
			}
		}
	}
	for(z=0;z<bc->zn;z++){
		side=bc->z[z].s;
		b=bc->z[z].b;
		if(side%2==1){
			bo=&d->vside[b].b[side];
			for(vi=bo->st;vi<=bo->en;vi+=bo->iv){
				s->ut[b][vi]=s->ut[b][vi+bo->in];
			}
		} else {
			bo=&d->hside[b].b[side];
			for(hi=bo->st;hi<=bo->en;hi+=bo->iv){
				s->vt[b][hi]=s->vt[b][hi+bo->in];
			}
		}
	}
	//Interfaces
	int vi1,vi2,hi1,hi2;
	int sign;
	int ic,b1,b2,s1,s2;
	Boundary*bo1,*bo2,*boc1,*boc2;
	for(i=0;i<bc->in;i++){
		s1=bc->i[i].s;
		b1=bc->i[i].b;
		b2=bc->bid[b1].id[s1];
		boc1=&d->cside[b1].b[s1];
		boc2=&d->cside[b2].b[s2];
		s2=(s1+2)%4;
		switch(s1){
		case 0: sign=-1; break;
		case 1: sign=1; break;
		case 2: sign=1; break;
		case 3: sign=-1; break;
		}
		if(s1%2==1){
			bo1=&d->vside[b1].b[s1];
			bo2=&d->vside[b2].b[s2];
			for(ic=0;ic<bo1->size;ic++){
				if(sign*s->u[b1][vi1]>=0){
					//inside b1
					b=b1;
					ci=boc1->st+ic*boc1->iv;
				} else {
					//inside b2
					b=b2;
					ci=boc2->st+ic*boc2->iv;
				}
				vi1=bo1->st+ic*bo1->iv;
				a2[0]=g->a21[b1][vi1];
				a2[1]=g->a22[b1][vi1];
				umom[0]=-mm->I11_1[b1][vi1]-mm->I12_1[b1][vi1]+
						mm->I14_1[b1][vi1]+mm->I15_1[b1][vi1];
				umom[1]=-mm->I11_2[b1][vi1]-mm->I12_2[b1][vi1]+
						mm->I14_2[b1][vi1]+mm->I15_2[b1][vi1];
				vi2=bo2->st+ic*bo2->iv;
				s->ut[b1][vi1]=s->ut[b2][vi2]=
						-cross(a2,umom)/g->Jcc[b][ci];
			}
		} else {
			bo1=&d->hside[b1].b[s1];
			bo2=&d->hside[b2].b[s2];
			for(ic=0;ic<bo1->size;ic++){
				if(sign*s->v[b1][hi1]>=0){
					//inside b1
					b=b1;
					ci=boc1->st+ic*boc1->iv;
				} else {
					//inside b2
					b=b2;
					ci=boc2->st+ic*boc2->iv;
				}
				hi1=bo1->st+ic*bo1->iv;
				a1[0]=g->a11[b1][hi1];
				a1[1]=g->a12[b1][hi1];
				vmom[0]=-mm->I21_1[b1][hi1]-mm->I22_1[b1][hi1]+
						mm->I24_1[b1][hi1]+mm->I25_1[b1][hi1];
				vmom[1]=-mm->I21_2[b1][hi1]-mm->I22_2[b1][hi1]+
						mm->I24_2[b1][hi1]+mm->I25_2[b1][hi1];
				hi2=bo2->st+ic*bo2->iv;
				s->vt[b1][hi1]=s->vt[b2][hi2]=
						cross(a1,vmom)/g->Jcc[b][ci];
			}
		}
	}
}
void fillTimeDerivative(State*st,
		Grid*g,BoundaryConditions*bc,Interfaces*is,WallDistances*wd,
		Momentum*mm,SpalartAllmaras*sa,
		Property*p,Switches*sw,Numerics*n,Dimension*d){
	//Calculates ut,vt,sanut based on the u,v,sanu in State*st.
	interpolateState(st,mm,g,is,bc,d);
	cartesianConvert(st,g,mm);
	cartesianDerivatives(g,mm,is);
	//convectionCentralDifference(g,mm,bc,is);
	convectionUpwindCC(d,g,st,mm,bc);
	if(sw->turbmod>0){
		computeSanut(st,sa,p,g,is,mm,d,bc,wd->cc);
		updateNut(st,sa,g,p,mm,bc,is);
	}
	viscosity(p,g,mm,bc,is);
	invariantTimeDerivativeCentralDifference(g,mm,st);
}
void fillRhs(Dimension*d,State*s,EFS*efs,double dt){
	//Adds time derivative term to rhs.
	int b,vi,hi,es;
	//Assign to RHS.
	for(b=0;b<d->tb;b++){
		es=d->e0[b];
		for(vi=0;vi<d->v[b].n;vi++){
			efs->rhs[es+vi]=s->u[b][vi]+dt*s->ut[b][vi];
		}
		for(hi=0;hi<d->h[b].n;hi++){
			efs->rhs[es+d->v[b].n+hi]=s->v[b][hi]+dt*s->vt[b][hi];
		}
	}
}
