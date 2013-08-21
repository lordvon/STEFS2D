void initTimeDerivative(State*s,Grid* g,Momentum* mm,BoundaryConditions* bc){
	int b,i,j,os,vi1,vi2,hi1,hi2;
	int vd0,vd1,hd0,hd3;
	char bc0,bc1,bc2,bc3;
	for(b=0;b<g->totalblocks;b++){
		vd0=g->vdim[b][0];
		vd1=g->vdim[b][1];
		hd0=g->hdim[b][0];
		hd3=g->hdim[b][3];
		bc0=bc->id[b][0];
		bc1=bc->id[b][1];
		bc2=bc->id[b][2];
		bc3=bc->id[b][3];
		//vertical boundaries
		os=vd0-1;
		for(j=0;j<vd1;j++){
			vi1=j*vd0;
			vi2=j*vd0+os;
			if((bc3=='f') | (bc3=='w')){ s->ut[b][vi1]=0; }
			if((bc1=='f') | (bc1=='w')){ s->ut[b][vi2]=0; }
		}
		//horizontal boundaries
		os=hd3-hd0;
		for(i=0;i<hd0;i++){
			hi1=i;
			hi2=i+os;
			if((bc0=='f') | (bc0=='w')){ s->vt[b][hi1]=0; }
			if((bc2=='f') | (bc2=='w')){ s->vt[b][hi2]=0; }
		}
	}
}
void invariantZeroGradient(State*st,Grid*g,BoundaryConditions*bc){
	int b,s,i;
	Boundary*bo;
	double**field;
	for(b=0;b<g->totalblocks;b++){
		for(s=0;s<4;s++){
			if(bc->id[b][s]=='z'){
				if(s%2==0){
					bo=&g->sidehh[b].b[s];
					field=st->v;
				} else {
					bo=&g->sidevv[b].b[s];
					field=st->u;
				}
				for(i=bo->st;i<=bo->en;i+=bo->iv){
					field[b][i]=field[b][i+bo->in];
				}
			}
		}
	}
}
void invariantFixed(State*st,Grid*g,BoundaryConditions*bc){
	int b,s,i;
	Boundary*bo;
	double a[2],sign;
	double**field,**met1,**met2;
	for(b=0;b<g->totalblocks;b++){
		for(s=0;s<4;s++){
			if(bc->id[b][s]=='f'){
				if(s%2==0){
					bo=&g->sidehh[b].b[s];
					field=st->v;
					met1=g->a11;
					met2=g->a12;
					sign=-1;
				} else {
					bo=&g->sidevv[b].b[s];
					field=st->u;
					met1=g->a21;
					met2=g->a22;
					sign=1;
				}
				for(i=bo->st;i<=bo->en;i+=bo->iv){
					a[0]=met1[b][i];
					a[1]=met2[b][i];
					field[b][i]=sign*cross(bc->fv[b],a);
				}
			}
		}
	}
}
void invariantWall(State*st,Grid*g,BoundaryConditions*bc){
	int b,s,i;
	Boundary*bo;
	double**field;
	for(b=0;b<g->totalblocks;b++){
		for(s=0;s<4;s++){
			if(bc->id[b][s]=='w'){
				if(s%2==0){
					bo=&g->sidehh[b].b[s];
					field=st->v;
				} else {
					bo=&g->sidevv[b].b[s];
					field=st->u;
				}
				for(i=bo->st;i<=bo->en;i+=bo->iv){
					field[b][i]=0;
				}
			}
		}
	}
}
void enforceInvariantBoundary(State*s,Grid*g,BoundaryConditions*bc){
	//Interfaces are not handled, since interface edges are part of the solution.
	//Zero Gradient
	invariantZeroGradient(s,g,bc);
	//Fixed
	invariantFixed(s,g,bc);
	//Wall
	invariantWall(s,g,bc);
}
void initState(Grid*g,Momentum*mm,LinearSystem*ls,State*s,Initial*in){
	int b,i,vi,hi,ci;
	int vd3,hd3,cd3;
	double a1[2],a2[2];
	for(b=0;b<g->totalblocks;b++){
		vd3=g->vdim[b][3];
		hd3=g->hdim[b][3];
		cd3=g->cdim[b][3];
		//velocity field initialization.
		for(vi=0;vi<vd3;vi++){
			a2[0]=g->a21[b][vi];
			a2[1]=g->a22[b][vi];
			s->u[b][vi]=cross(in->initcart,a2);
		}
		for(hi=0;hi<hd3;hi++){
			a1[0]=g->a11[b][hi];
			a1[1]=g->a12[b][hi];
			s->v[b][hi]=-cross(in->initcart,a1);
		}
		for(ci=0;ci<cd3;ci++){
			s->sanu[b][ci]=in->sanu0;
		}
		/*
		//init turbulence quantities
		for(ci=0;ci<cd3;ci++){
			mm->nut[b][ci]=in->nut0;
		}
		for(ci=0;ci<cd3;ci++){
			s->tke[b][ci]=in->tke0;
			s->tdr[b][ci]=in->tdr0;
		}
		*/
		//discrete streamfunction vector
		for(i=0;i<ls->ccols;i++){
			ls->s[i]=0;
		}
	}
}
void initFluxes(Grid* g,Momentum* mm){
	int b,vi,hi;
	int vd3,hd3;
	for(b=0;b<g->totalblocks;b++){
		vd3=g->vdim[b][3];
		hd3=g->hdim[b][3];
		for(vi=0;vi<vd3;vi++){
			mm->I11_1[b][vi]=0;
			mm->I11_2[b][vi]=0;
			mm->I12_1[b][vi]=0;
			mm->I12_2[b][vi]=0;
			mm->I14_1[b][vi]=0;
			mm->I14_2[b][vi]=0;
			mm->I15_1[b][vi]=0;
			mm->I15_2[b][vi]=0;
		}
		for(hi=0;hi<hd3;hi++){
			mm->I21_1[b][hi]=0;
			mm->I21_2[b][hi]=0;
			mm->I22_1[b][hi]=0;
			mm->I22_2[b][hi]=0;
			mm->I24_1[b][hi]=0;
			mm->I24_2[b][hi]=0;
			mm->I25_1[b][hi]=0;
			mm->I25_2[b][hi]=0;
		}
	}
}
void cctoxxz(double**cc,double**xx,Grid*g,char**bcid,int***sideinfo){
	int b,st,iv,en,side,adj,ci0,civ,cif;
	int ci,xi;
	int xd0,cd0;
	for(b=0;b<g->totalblocks;b++){
		xd0=g->xdim[b][0];
		cd0=g->cdim[b][0];
		for(side=0;side<4;side++){
			if(bcid[b][side]=='z'){
				switch(side){
				case 0: adj=0;		break;
				case 1: adj=-1; 	break;
				case 2: adj=-cd0; 	break;
				case 3: adj=0;		break;
				}
				switch(side%2){
				case 0: civ=1;		break;
				case 1: civ=cd0; 	break;
				}
				st=sideinfo[b][side][0];
				iv=sideinfo[b][side][1];
				en=sideinfo[b][side][2];
				ci0=st-st/xd0+adj;
				cif=en-en/xd0+adj-civ;
				xx[b][st]=cc[b][ci0];
				xx[b][en]=cc[b][cif];
				for(xi=st+iv;xi<=en-iv;xi+=iv){
					ci=xi-xi/xd0+adj;
					xx[b][xi]=0.5*(cc[b][ci]+cc[b][ci-civ]);
				}
			}
		}
	}
}
void cctoxxwfsymlap(double**cc,double**xx,Grid*g,char**bcid,int***sideinfo){
	//only fixed value of zero currently available.
	int b,st,iv,en,side;
	int xi;
	for(b=0;b<g->totalblocks;b++){
		for(side=0;side<4;side++){
			if((bcid[b][side]=='w') | (bcid[b][side]=='f')){
				st=sideinfo[b][side][0];
				iv=sideinfo[b][side][1];
				en=sideinfo[b][side][2];
				for(xi=st;xi<=en;xi+=iv){
					xx[b][xi]=0;
				}
			}
		}
	}
}
double scalarSanuToNut(SpalartAllmaras*sa,Property*p,double sanu){
	double chi3=pow(sanu*p->nuinv,3);
	return sanu*chi3/(chi3+sa->cv1_3);
}
void cctoxxwfSA(double**cc,double**xx,
		SpalartAllmaras*sa,Property*p,
		Grid*g,char**bcid,int***sideinfo){
	int b,st,iv,en,side;
	int xi;
	for(b=0;b<g->totalblocks;b++){
		for(side=0;side<4;side++){
			if(bcid[b][side]=='w'){
				st=sideinfo[b][side][0];
				iv=sideinfo[b][side][1];
				en=sideinfo[b][side][2];
				for(xi=st;xi<=en;xi+=iv){
					xx[b][xi]=0;
				}
			}
			if(bcid[b][side]=='f'){
				st=sideinfo[b][side][0];
				iv=sideinfo[b][side][1];
				en=sideinfo[b][side][2];
				for(xi=st;xi<=en;xi+=iv){
					xx[b][xi]=scalarSanuToNut(sa,p,sa->fixed[b]);
				}
			}
		}
	}
}
int checkForW(int*icorner,char**bcids){
	//Checks to see if wall is present at a given corner.
	int wall=0;
	int c,s1,s2,b;
	for(c=0;c<4;c++){
		b=icorner[c];
		if(b>=0){
			s1=c;
			s2=(c+3)%4;
			if((bcids[b][s1]=='w') | (bcids[b][s2]=='w')){
				wall=1;
				break;
			}
		}
	}
	return wall;
}
void enforcecctoxxicorners(double**cc,double**xx,char**bcids,Grid*g,Interfaces*is){
	int ic,c,b,ci,xi;
	double interp,sum;int n;//interpolation=sum/n
	int tic=is->totalicorners;
	for(ic=0;ic<tic;ic++){
		if(checkForW(is->icorners[ic],bcids)){
			for(c=0;c<4;c++){
				b=is->icorners[ic][c];
				if(b>=0){
					xi=getcornerxi(g,b,c);
					xx[b][xi]=0;
				}
			}
		} else {
			sum=0;n=0;
			for(c=0;c<4;c++){
				b=is->icorners[ic][c];
				if(b>=0){
					ci=getcornercellci(g,b,c);
					sum+=cc[b][ci];n++;
				}
			}
			interp=sum/n;
			for(c=0;c<4;c++){
				b=is->icorners[ic][c];
				if(b>=0){
					xi=getcornerxi(g,b,c);
					xx[b][xi]=interp;
				}
			}
		}
	}
}
void cctoxxi(double**cc,double**xx,char**bcids,Grid*g,Interfaces*is){
	int i,side;
	int ip,b1,b2,xi1,xi2;
	int adj1,adj2,ci1,ci2;
	int cd01,cd02,civ1,civ2;
	for(i=0;i<is->totalinterfaces;i++){
		b1=is->block[i];
		b2=is->adjacentBlock[i];
		side=is->side[i];
		cd01=g->cdim[b1][0];
		cd02=g->cdim[b2][0];
		switch(side){
		case 0: adj1=0;adj2=-cd02; 	break;
		case 1: adj1=-1;adj2=0; 	break;
		case 2: adj1=-cd01;adj2=0;	break;
		case 3: adj1=0;adj2=-1; 	break;
		}
		switch(side%2){
		case 0:	civ1=civ2=1; break;
		case 1: civ1=cd01;civ2=cd02; break;
		}
		//start and end points.
		xi1=is->blockx[i][0]-g->xstart[b1];ci1=xi1-xi1/g->xdim[b1][0]+adj1;
		xi2=is->adjacentBlockx[i][0]-g->xstart[b2];ci2=xi2-xi2/g->xdim[b2][0]+adj2;
		xx[b1][xi1]=xx[b2][xi2]=0.5*(
				cc[b1][ci1]+cc[b2][ci2]);
		xi1=is->blockx[i][is->xsize[i]-1]-g->xstart[b1];ci1=xi1-xi1/g->xdim[b1][0]+adj1;
		xi2=is->adjacentBlockx[i][is->xsize[i]-1]-g->xstart[b2];ci2=xi2-xi2/g->xdim[b2][0]+adj2;
		xx[b1][xi1]=xx[b2][xi2]=0.5*(
				cc[b1][ci1-civ1]+cc[b2][ci2-civ2]);
		//the rest
		for(ip=1;ip<is->xsize[i]-1;ip++){
			xi1=is->blockx[i][ip]-g->xstart[b1];
			xi2=is->adjacentBlockx[i][ip]-g->xstart[b2];
			ci1=xi1-xi1/g->xdim[b1][0]+adj1;
			ci2=xi2-xi2/g->xdim[b2][0]+adj2;
			xx[b1][xi1]=xx[b2][xi2]=0.25*(
					cc[b1][ci1]+cc[b1][ci1-civ1]+
					cc[b2][ci2]+cc[b2][ci2-civ2]);
		}
	}
	enforcecctoxxicorners(cc,xx,bcids,g,is);
}
void cctoxxInterior(double**cc,double**xx,Grid*g){
	int b,i,j,xi,ci;
	int xd0,xd1,cd0;//,cd1;
	for(b=0;b<g->totalblocks;b++){
		xd0=g->xdim[b][0];
		xd1=g->xdim[b][1];
		cd0=g->cdim[b][0];
		//cd1=g->cdim[b][1];
		//interior
		for(j=1;j<xd1-1;j++){
			for(i=1;i<xd0-1;i++){
				xi=i+j*xd0;
				ci=xi-j;
				xx[b][xi]=0.25*(cc[b][ci]+
						cc[b][ci-1]+
						cc[b][ci-cd0]+
						cc[b][ci-cd0-1]);
			}
		}
	}
}
void cctoxxAllZG(double**cc,double**xx,Grid*g,Interfaces*is,
		int***sideinfo){
	int b,s;
	char **bcid=malloc(sizeof(char*)*g->totalblocks);
	for(b=0;b<g->totalblocks;b++){
		bcid[b]=malloc(sizeof(char)*4);
		for(s=0;s<4;s++){
			bcid[b][s]='z';
		}
	}
	cctoxxInterior(cc,xx,g);
	cctoxxz(cc,xx,g,bcid,sideinfo);
	cctoxxi(cc,xx,bcid,g,is);
}
void cctoxxCustom(double**cc,double**xx,Grid*g,Interfaces*is,
		char**bcid,int***sideinfo){
	cctoxxInterior(cc,xx,g);
	cctoxxi(cc,xx,bcid,g,is);
	cctoxxz(cc,xx,g,bcid,sideinfo);
	cctoxxwfsymlap(cc,xx,g,bcid,sideinfo);
}
void cctoxxSL(double**cc,double**xx,Grid*g,BoundaryConditions*bc,Interfaces*is){
	cctoxxInterior(cc,xx,g);
	cctoxxi(cc,xx,bc->id,g,is);
	cctoxxz(cc,xx,g,bc->id,bc->side);
	cctoxxwfsymlap(cc,xx,g,bc->id,bc->side);
}
void cctoxxSA(double**cc,double**xx,
		SpalartAllmaras*sa,Property*p,
		Grid*g,BoundaryConditions*bc,Interfaces*is){
	cctoxxInterior(cc,xx,g);
	cctoxxi(cc,xx,bc->id,g,is);
	cctoxxz(cc,xx,g,bc->id,bc->side);
	cctoxxwfSA(cc,xx,sa,p,g,bc->id,bc->side);
}
