void fillI11InterfaceUpwindCC(Dimension*d,State*st,Momentum*mm,
		BoundaryConditions*bc){
	//Sides are 1 or 3.
	int i,vi1,vi2;
	int sign;
	int ic,b1,b2,s1,s2;
	int b,interface,inner;
	Boundary*bo1,*bo2;
	for(i=0;i<bc->in;i++){
		s1=bc->i[i].s;
		if(s1%2==1){
			b1=bc->i[i].b;
			b2=bc->bid[b1].id[s1];
			s2=(s1+2)%4;
			bo1=&d->vside[b1].b[s1];
			bo2=&d->vside[b2].b[s2];
			if(s1==1){ sign=1; }
			else { sign=-1; }
			for(ic=0;ic<bo1->size;ic++){
				vi1=bo1->st+ic*bo1->iv;
				vi2=bo2->st+ic*bo2->iv;
				if(sign*st->u[b1][vi1]>=0){
					//inside b1
					b=b1;
					interface=vi1;
					inner=vi1+bo1->in;
				} else {
					//inside b2
					b=b2;
					interface=vi2;
					inner=vi2+bo2->in;
				}
				mm->I11_1[b1][vi1]=mm->I11_1[b2][vi2]=
						sign*(	mm->ucartvv[b][interface]*st->u[b][interface]-
								mm->ucartvv[b][inner]*st->u[b][inner]);
				mm->I11_2[b1][vi1]=mm->I11_2[b2][vi2]=
						sign*(	mm->vcartvv[b][interface]*st->u[b][interface]-
								mm->vcartvv[b][inner]*st->u[b][inner]);
			}
		}
	}
}
void fillI11InteriorUpwindCC(Dimension*d,State*st,Momentum*mm){
	int b,i,j,vi;
	int vd1,vd0;
	for(b=0;b<d->tb;b++){
		vd0=d->v[b].i;
		vd1=d->v[b].j;
		for(j=0;j<vd1;j++){
			for(i=1;i<vd0-1;i++){
				vi=i+j*vd0;
				if(st->u[vi]>=0){
					mm->I11_1[b][vi]=
							mm->ucartvv[b][vi]*st->u[b][vi]-
							mm->ucartvv[b][vi-1]*st->u[b][vi-1];
				} else {
					mm->I11_1[b][vi]=
							mm->ucartvv[b][vi+1]*st->u[b][vi+1]-
							mm->ucartvv[b][vi]*st->u[b][vi];
				}
			}
		}
		for(j=0;j<vd1;j++){
			for(i=1;i<vd0-1;i++){
				vi=i+j*vd0;
				if(st->u[vi]>=0){
					mm->I11_2[b][vi]=
							mm->vcartvv[b][vi]*st->u[b][vi]-
							mm->vcartvv[b][vi-1]*st->u[b][vi-1];
				} else {
					mm->I11_2[b][vi]=
							mm->vcartvv[b][vi+1]*st->u[b][vi+1]-
							mm->vcartvv[b][vi]*st->u[b][vi];
				}
			}
		}
	}
}
void fillI22InterfaceUpwindCC(Dimension*d,State*st,Momentum*mm,
		BoundaryConditions*bc){
	//Sides are 0 or 2.
	int i,hi1,hi2;
	int sign;
	int ic,b1,b2,s1,s2;
	int b,interface,inner;
	Boundary*bo1,*bo2;
	for(i=0;i<bc->in;i++){
		s1=bc->i[i].s;
		if(s1%2==0){
			b1=bc->i[i].b;
			b2=bc->bid[b1].id[s1];
			s2=(s1+2)%4;
			bo1=&d->hside[b1].b[s1];
			bo2=&d->hside[b2].b[s2];
			if(s1==2){ sign=1; }
			else { sign=-1; }
			for(ic=0;ic<bo1->size;ic++){
				hi1=bo1->st+ic*bo1->iv;
				hi2=bo2->st+ic*bo2->iv;
				if(sign*st->v[b1][hi1]>=0){
					//inside b1
					b=b1;
					interface=hi1;
					inner=hi1+bo1->in;
				} else {
					//inside b2
					b=b2;
					interface=hi2;
					inner=hi2+bo2->in;
				}
				mm->I22_1[b1][hi1]=mm->I22_1[b2][hi2]=
						sign*(	mm->ucarthh[b][interface]*st->v[b][interface]-
								mm->ucarthh[b][inner]*st->v[b][inner]);
				mm->I22_2[b1][hi1]=mm->I22_2[b2][hi2]=
						sign*(	mm->vcarthh[b][interface]*st->v[b][interface]-
								mm->vcarthh[b][inner]*st->v[b][inner]);
			}
		}
	}
}
void fillI22InteriorUpwindCC(State*st,Dimension*d,Momentum*mm){
	int b,i,j,hi;
	int hd1,hd0;
	for(b=0;b<d->tb;b++){
		hd0=d->h[b].i;
		hd1=d->h[b].j;
		for(j=1;j<hd1-1;j++){
			for(i=0;i<hd0;i++){
				hi=i+j*hd0;
				if(st->v[b][hi]>=0){
					mm->I22_1[b][hi]=
							mm->ucarthh[b][hi]*st->v[b][hi]-
							mm->ucarthh[b][hi-hd0]*st->v[b][hi-hd0];
				} else {
					mm->I22_1[b][hi]=
							mm->ucarthh[b][hi+hd0]*st->v[b][hi+hd0]-
							mm->ucarthh[b][hi]*st->v[b][hi];
				}
			}
		}
		for(j=1;j<hd1-1;j++){
			for(i=0;i<hd0;i++){
				hi=i+j*hd0;
				if(st->v[b][hi]>=0){
					mm->I22_2[b][hi]=
							mm->vcarthh[b][hi]*st->v[b][hi]-
							mm->vcarthh[b][hi-hd0]*st->v[b][hi-hd0];
				} else {
					mm->I22_2[b][hi]=
							mm->vcarthh[b][hi+hd0]*st->v[b][hi+hd0]-
							mm->vcarthh[b][hi]*st->v[b][hi];
				}
			}
		}
	}
}
void fillI11UpwindCC(Dimension*d,Grid*g,State*s,Momentum*mm,BoundaryConditions*bc){
	fillI11InteriorUpwindCC(d,s,mm);
	fillI11ZeroGradient(g,mm,bc);
	fillI11InterfaceUpwindCC(d,s,mm,bc);
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
void fillI22UpwindCC(Dimension*d,Grid* g,State*s,Momentum*mm,BoundaryConditions*bc){
	fillI22InteriorUpwindCC(s,d,mm);
	fillI22ZeroGradient(g,mm,bc);
	fillI22InterfaceUpwindCC(d,s,mm,bc);
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
