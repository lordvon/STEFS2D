void upwindCcField(Dimension*d,Momentum*mm,BoundaryConditions*bc,
		double*fixed,double**cc,double**ccxi,double**cceta){
	int b,i,j,ci;
	//interior
	for(b=0;b<d->tb;b++){
		//xi
		for(j=0;j<d->c[b].j;j++){
			for(i=1;i<d->c[b].i-1;i++){
				ci=i+j*d->c[b].i;
				if(mm->ucc[b][ci]>=0){
					ccxi[b][ci]=cc[b][ci]-cc[b][ci-1];
				} else {
					ccxi[b][ci]=cc[b][ci+1]-cc[b][ci];
				}
			}
		}
		//eta
		for(j=1;j<d->c[b].j-1;j++){
			for(i=0;i<d->c[b].i;i++){
				ci=i+j*d->c[b].i;
				if(mm->vcc[b][ci]>=0){
					cceta[b][ci]=cc[b][ci]-cc[b][ci-d->c[b].i];
				} else {
					cceta[b][ci]=cc[b][ci+d->c[b].i]-cc[b][ci];
				}
			}
		}
	}
	//boundaries.
	int s,sign;
	Boundary*bo;
	double**uvcc,**ccdiff;
	//Zero Gradient
	int z;
	for(z=0;z<bc->zn;z++){
		b=bc->z[z].b;
		s=bc->z[z].s;
		bo=&d->cside[b].b[s];
		if(s%2==0){
			uvcc=mm->vcc;
			ccdiff=cceta;
			if(s==2){ sign=1; }
			else { sign=-1; }
		} else {
			uvcc=mm->ucc;
			ccdiff=ccxi;
			if(s==1){ sign=1; }
			else { sign=-1; }
		}
		for(ci=bo->st;ci<=bo->en;ci+=bo->iv){
			if(sign*uvcc[b][ci]>=0){
				ccdiff[b][ci]=sign*(cc[b][ci]-cc[b][ci+bo->in]);
			} else {
				ccdiff[b][ci]=0;
			}
		}
	}
	//Interfaces
	int ic,ci1,ci2,b1,b2,s1,s2,sign1,sign2;
	Boundary*bo1,*bo2;
	for(i=0;i<bc->in;i++){
		b1=bc->i[i].b;
		s1=bc->i[i].s;
		b2=bc->bid[b1].id[s1];
		s2=(s1+2)%4;
		bo1=&d->cside[b1].b[s1];
		bo2=&d->cside[b2].b[s2];
		if(s1%2==0){
			uvcc=mm->vcc;
			ccdiff=cceta;
			if(s1==2){ sign1=1; }
			else { sign1=-1; }
		} else {
			uvcc=mm->ucc;
			ccdiff=ccxi;
			if(s1==1){ sign1=1; }
			else { sign1=-1; }
		}
		for(ic=0;ic<bo1->size;ic++){
			ci1=bo1->st+ic*bo1->iv;
			if(sign1*uvcc[b1][ci1]>=0){
				ccdiff[b1][ci1]=sign1*(cc[b1][ci1]-cc[b1][ci1+bo1->in]);
			} else {
				ci2=bo2->st+ic*bo2->iv;
				ccdiff[b1][ci1]=sign1*(cc[b2][ci2]-cc[b1][ci1]);
			}
		}
		sign2=-sign1;
		for(ic=0;ic<bo2->size;ic++){
			ci2=bo2->st+ic*bo2->iv;
			if(sign2*uvcc[b2][ci2]>=0){
				ccdiff[b2][ci2]=sign2*(cc[b2][ci2]-cc[b2][ci2+bo2->in]);
			} else {
				ci1=bo1->st+ic*bo1->iv;
				ccdiff[b2][ci2]=sign2*(cc[b1][ci1]-cc[b2][ci2]);
			}
		}
	}
	//Fixed
	int f;
	for(f=0;f<bc->fn;f++){
		b=bc->f[f].b;
		s=bc->f[f].s;
		bo=&d->cside[b].b[s];
		if(s%2==0){
			uvcc=mm->vcc;
			ccdiff=cceta;
			if(s==2){ sign=1; }
			else { sign=-1; }
		} else {
			uvcc=mm->ucc;
			ccdiff=ccxi;
			if(s==1){ sign=1; }
			else { sign=-1; }
		}
		for(ci=bo->st;ci<=bo->en;ci+=bo->iv){
			if(sign*uvcc[b][ci]>=0){
				ccdiff[b][ci]=sign*(cc[b][ci]-cc[b][ci+bo->in]);
			} else {
				ccdiff[b][ci]=sign*2*(fixed[b]-cc[b][ci]);
			}
		}
	}
	//Walls
	int w;
	for(w=0;w<bc->wn;w++){
		b=bc->w[w].b;
		s=bc->w[w].s;
		bo=&d->cside[b].b[s];
		if(s%2==0){
			uvcc=mm->vcc;
			ccdiff=cceta;
			if(s==2){ sign=1; }
			else { sign=-1; }
		} else {
			uvcc=mm->ucc;
			ccdiff=ccxi;
			if(s==1){ sign=1; }
			else { sign=-1; }
		}
		for(ci=bo->st;ci<=bo->en;ci+=bo->iv){
			if(sign*uvcc[b][ci]>=0){
				ccdiff[b][ci]=sign*(cc[b][ci]-cc[b][ci+bo->in]);
			} else {
				ccdiff[b][ci]=sign*2*(-cc[b][ci]);
			}
		}
	}
}
