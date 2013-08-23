void fillCInterface02(CSR*C,Dimension*d,BoundaryConditions*bc,int b1,int s1){
	//Assumes b1 is the primary block whose entity indices will be retained.
	int b2,s2,adj;
	int xi1,hi2,vi2;
	int hrow,vrow,newx;
	int i;

	b2=bc->bid[b1].id[s1];
	s2=(s1+2)%4;

	//Step through h edges.
	xi1=d->xside[b1].b[s1].st-1;
	hi2=d->hside[b2].b[s2].st-1;
	for(i=0;i<d->hside[b1].b[s1].size;i++){
		xi1++;hi2++;
		newx=d->x0[b1]+xi1;
		hrow=d->e0[b2]+d->v[b2].n+hi2;
		if((C->ri[2*hrow]!=hrow) | (C->ri[2*hrow+1]!=hrow)){ printf("hrow mismatch!\n"); };
		C->ci[2*hrow]=newx;
		C->ci[2*hrow+1]=newx+1;
	}
	//Step through v edges.
	xi1=d->xside[b1].b[s1].st-1;
	vi2=d->vside[b2].b[s2].st-1;
	if(s2==0){ adj=1; } else { adj=0; }
	for(i=0;i<d->vside[b1].b[s1].size;i++){
		xi1++;vi2++;
		newx=d->x0[b1]+xi1;
		vrow=d->e0[b2]+vi2;
		if((C->ri[2*vrow]!=vrow) | (C->ri[2*vrow+1]!=vrow)){ printf("vrow mismatch!\n"); };
		C->ci[2*vrow+adj]=newx;
	}
}
void fillCInterface13(CSR*C,Dimension*d,BoundaryConditions*bc,int b1,int s1){
	//Assumes b1 is the primary block whose entity indices will be retained.
	int b2,s2,adj;
	int xi1,hi2,vi2;
	int xiv1,hiv2,viv2;
	int hrow,vrow,newx;
	int i;

	b2=bc->bid[b1].id[s1];
	s2=(s1+2)%4;
	hiv2=d->hside[b2].b[s2].iv;
	viv2=d->vside[b2].b[s2].iv;
	xiv1=d->xside[b1].b[s1].iv;

	vi2=d->vside[b2].b[s2].st-viv2;
	xi1=d->xside[b1].b[s1].st-xiv1;
	//Step through v edges.
	for(i=0;i<d->vside[b1].b[s1].size;i++){
		xi1+=xiv1;vi2+=viv2;
		newx=d->x0[b1]+xi1;
		vrow=d->e0[b2]+vi2;
		if((C->ri[2*vrow]!=vrow) | (C->ri[2*vrow+1]!=vrow)){ printf("vrow mismatch!\n"); };
		C->ci[2*vrow]=newx+xiv1;
		C->ci[2*vrow+1]=newx;
	}
	//Step through h edges.
	if(s2==1){ adj=1; } else { adj=0; }
	hi2=d->hside[b2].b[s2].st-hiv2;
	xi1=d->xside[b1].b[s1].st-xiv1;
	for(i=0;i<d->hside[b1].b[s1].size;i++){
		xi1+=xiv1;hi2+=hiv2;
		newx=d->x0[b1]+xi1;
		hrow=d->e0[b2]+d->v[b2].n+hi2;
		if((C->ri[2*hrow]!=hrow) | (C->ri[2*hrow+1]!=hrow)){ printf("hrow mismatch!\n"); };
		C->ci[2*hrow+adj]=newx;
	}
}
void fillC(CSR*C,Dimension*d,BoundaryConditions*bc){
	//In C, nodes are in columns and edges are in rows.
	int b,vi,hi,xi,entry=0;
	int ii,b1,s1;
	C->n=2*d->te;
	mallocC(C,C->n);
	for(b=0;b<d->tb;b++){
		for(vi=0;vi<d->v[b].n;vi++){
			xi=vi;
			C->ri[entry]=d->e0[b]+vi;
			C->ci[entry]=d->x0[b]+xi+d->x[b].i;
			C->v[entry]=1;
			entry++;
			C->ri[entry]=d->e0[b]+vi;
			C->ci[entry]=d->x0[b]+xi;
			C->v[entry]=-1;
			entry++;
		}
		for(hi=0;hi<d->h[b].n;hi++){
			xi=hi+hi/d->h[b].i;
			C->ri[entry]=d->e0[b]+d->v[b].n+hi;
			C->ci[entry]=d->x0[b]+xi;
			C->v[entry]=1;
			entry++;
			C->ri[entry]=d->e0[b]+d->v[b].n+hi;
			C->ci[entry]=d->x0[b]+xi+1;
			C->v[entry]=-1;
			entry++;
		}
	}
	//Rearrange some entries to account for interfaces.
	for(ii=0;ii<bc->in;ii++){
		b1=bc->i[ii].b;
		s1=bc->i[ii].s;
		if(s1%2==0){ 	fillCInterface02(C,d,bc,b1,s1); }
		else { 			fillCInterface13(C,d,bc,b1,s1); }
	}
}
