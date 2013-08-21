
void symlapInteriorNaive(SymLap*s,Grid*g,Interfaces*is){
	//currently, doesnt do fixed boundary values.
	int b,tc;
	int i,j,ci,hi,vi,xi;
	int cd0,cd1,cd3,hd0,xd0;
	//naming:
	//p: phi
	//underscore: negative
	//number: displacement direction relative to cell center.
	double p,p1,p_1,p2,p_2;
	double p12,p_1_2,p_12,p1_2;
	double *B11,*B22,*B12;
	int md,cs;
	int bottom,top,right,left;
	tc=g->totalcells;
	for(b=0;b<g->totalblocks;b++){
		B11=g->B11[b];
		B22=g->B22[b];
		B12=g->B12[b];
		cd0=g->cdim[b][0];
		cd1=g->cdim[b][1];
		cd3=g->cdim[b][3];
		hd0=g->hdim[b][0];
		xd0=g->xdim[b][0];
		//interior
		for(j=0;j<cd1;j++){
			for(i=0;i<cd0;i++){
				hi=ci=i+j*cd0;
				xi=vi=ci+j;

				p1=B11[vi+1];
				p2=B22[hi+hd0];
				p_1=B11[vi];
				p_2=B22[hi];

				p12=0.5*B12[xi+1+xd0];
				p1_2=-0.5*B12[xi+1];
				p_12=-0.5*B12[xi+xd0];
				p_1_2=0.5*B12[xi];

				p=-p1-p_1-p_2-p2
						-p_12-p1_2-p_1_2-p12;

				cs=g->cstart[b]+ci;
				md=(tc+1)*cs;
				s->A[md]=s->d0[b][ci]+=p;

				bottom=j==0;
				top=j==cd1-1;
				left=i==0;
				right=i==cd0-1;

				if(!top){ 				s->A[md+cd0]+=p2;if(ci<cd3-cd0){ s->d3[b][ci]+=p2; } }
				if(!bottom){ 			s->A[md-cd0]+=p_2; }
				if(!right){ 			s->A[md+1]+=p1;if(ci<cd3-1){ s->d1[b][ci]+=p1; } }
				if(!left){ 				s->A[md-1]+=p_1; }
				if(!left & !bottom){ 	s->A[md-cd0-1]+=p_1_2; }
				if(!right & !bottom){ 	s->A[md-cd0+1]+=p1_2; }
				if(!top & !right){ 		s->A[md+cd0+1]+=p12;if(ci<cd3-cd0-1){ s->d4[b][ci]+=p12; } }
				if(!top & !left){ 		s->A[md+cd0-1]+=p_12;if(ci<cd3-cd0+1){ s->d2[b][ci]+=p_12; } }
			}
		}
	}
}
void symlapBoundaryNaive(SymLap*s,Grid*g,BoundaryConditions*bc){
	int b,side,adj,sign;
	int xi,hi,vi,ci,xd0;
	int st,iv,en,coli,tc;
	double bbb,met1,met2;
	tc=g->totalcells;
	//fixed/wall boundaries, to be called last.
	for(b=0;b<g->totalblocks;b++){
		xd0=g->xdim[b][0];
		for(side=0;side<4;side++){
			if((bc->id[b][side]=='f') & (bc->id[b][side]=='w')){
				st=bc->side[b][side][0];
				iv=bc->side[b][side][1];
				en=bc->side[b][side][2];
				switch(side){
				case 0: adj=0;				sign=1;break;
				case 1: adj=-1;				sign=-1;break;
				case 2: adj=-g->cdim[b][0];	sign=-1;break;
				case 3: adj=0;				sign=1;break;
				}
				for(xi=st;xi<=en;xi++){
					ci=hi=xi-xi/xd0;
					vi=xi;
					switch(side%2){
					case 0: bbb=g->B22[b][hi]; break;
					case 1: bbb=g->B11[b][vi]; break;
					}
					coli=g->cstart[b]+ci+adj;
					met1=g->B12[b][xi];
					met2=g->B12[b][xi+iv];
					//additions to center cell
					s->A[coli*(tc+1)]+=0.5*sign*(met2-met1);
					s->A[coli*(tc+1)]+=-bbb;
					s->d0[b][ci]+=0.5*sign*(met2-met1);
					s->d0[b][ci]+=-bbb;
				}
			}
		}
	}
}
void symlapInterfaceNaive(SymLap*s,Grid*g,Interfaces* is){
	int i,b1,b2,side,ip,tc;
	int xi1,xi2,vi1,vi2,ci1,ci2,hi1,hi2;
	int civ1,civ2,xiv,adj1,adj2,rowi,coli,slin;
	double bbb;
	//int i1,i2;
	//int boundary;
	tc=g->totalcells;
	//interfaces, currently ignores corners of multiple interfaces.
	for(i=0;i<is->totalinterfaces;i++){
		b1=is->block[i];
		side=is->side[i];
		b2=is->adjacentBlock[i];
		for(ip=0;ip<is->xsize[i]-1;ip++){
			xi1=vi1=is->blockx[i][ip]-g->xstart[b1];
			xi2=vi2=is->adjacentBlockx[i][ip]-g->xstart[b2];
			ci1=hi1=xi1-xi1/g->xdim[b1][0];
			ci2=hi2=xi2-xi2/g->xdim[b2][0];
			if(side%2==0){
				civ1=civ2=xiv=1;
				bbb=g->B22[b1][hi1];
			} else {
				xiv=g->xdim[b1][0];
				civ1=g->cdim[b1][0];
				civ2=g->cdim[b2][0];
				bbb=g->B11[b1][vi1];
			}
			if(side==0){
				adj1=0;
				adj2=-g->cdim[b2][0];
			} else if(side==2){
				adj1=-g->cdim[b1][0];
				adj2=0;
			} else if(side==1){
				adj1=-1;
				adj2=0;
			} else if(side==3){
				adj1=0;
				adj2=-1;
			}
			switch(0){
			case 0:
				//boundary=(ip>0) & (ip<is->xsize[i]-2);
				//b1
				rowi=(g->cstart[b1]+ci1+adj1)*tc;
				coli=g->cstart[b2]+ci2+adj2;
				slin=rowi+coli;
				if(ip>0){		s->A[slin-civ2]=	g->B12[b1][xi1]; }
				s->A[slin]=			bbb;
				if(ip<is->xsize[i]-2){	 	s->A[slin+civ2]=	g->B12[b1][xi1+xiv]; }
				//b2
				rowi=(g->cstart[b2]+ci2+adj2)*tc;
				coli=g->cstart[b1]+ci1+adj1;
				slin=rowi+coli;
				if(ip>0){		s->A[slin-civ1]=	g->B12[b1][xi1]; }
				s->A[slin]=			bbb;
				if(ip<is->xsize[i]-2){	 	s->A[slin+civ1]=	g->B12[b1][xi1+xiv]; }
				break;
			case 1:
				//column method.
				/*
				i1=g->cstart[b1]+ci1+adj1;
				i2=g->cstart[b2]+ci2+adj2;
				//printf("%d,%d\n",i1,i2);
				//s->iloc[i][ip][0]=i1;
				//s->iloc[i][ip][1]=i2;
				if((ip>0) & (ip<is->xsize[i]-2)){
					s->iloc[i][ip][2]=civ2;
					s->valcol[i][ip][0]=g->B12[b1][xi1];
					s->valcol[i][ip][2]=g->B12[b1][xi1+xiv];
					//if(ip>0){
						s->A[i1+(i2-civ2)*tc]=g->B12[b1][xi1];
						s->A[i1*tc+i2-civ2]=g->B12[b1][xi1];
					//}
					//if(ip<is->xsize[i]-2){
						s->A[i1+(i2+civ2)*tc]=g->B12[b1][xi1+xiv];
						s->A[i1*tc+i2+civ2]=g->B12[b1][xi1+xiv];
					//}
				} else {
					s->iloc[i][ip][2]=0;
				}
				s->valcol[i][ip][1]=bbb;
				s->A[i1+i2*tc]=bbb;
				s->A[i1*tc+i2]=bbb;
				*/
				break;
			}
		}
	}
}

void slconstructNaive(SymLap*s,Grid*g,BoundaryConditions*bc,Interfaces*i){
	prepSymlap(s,g,i);
	int tc=g->totalcells;
	s->A=malloc(sizeof(double)*tc*tc);
	initv(s->A,tc*tc,0);
	symlapInteriorNaive(s,g,i);
	symlapInterfaceNaive(s,g,i);
	symlapBoundaryNaive(s,g,bc);
}
