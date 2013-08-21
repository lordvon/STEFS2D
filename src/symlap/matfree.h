#ifndef MATFREE_H_
#define MATFREE_H_

void symlapInterior(SymLap*s,Grid*g,Interfaces*is){
	//currently, doesnt do fixed boundary values.
	int b;
	int i,j,ci,hi,vi,xi;
	int cd0,cd1,cd3,hd0,xd0;
	//naming:
	//p: phi
	//underscore: negative
	//number: displacement direction relative to cell center.
	double p,p1,p_1,p2,p_2;
	double p12,p_1_2,p_12,p1_2;
	double *B11,*B22,*B12;
	int top,right,left;
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

				s->d0[b][ci]+=p;

				top=j==cd1-1;
				left=i==0;
				right=i==cd0-1;

				if(!top){ 				if(ci<cd3-cd0){ s->d3[b][ci]+=p2; } }
				if(!right){ 			if(ci<cd3-1){ s->d1[b][ci]+=p1; } }
				if(!top & !right){ 		if(ci<cd3-cd0-1){ s->d4[b][ci]+=p12; } }
				if(!top & !left){ 		if(ci<cd3-cd0+1){ s->d2[b][ci]+=p_12; } }
			}
		}
	}
}
int countSymlapBoundaryCorners(Interfaces*is){
	int count;
	int ic,c,b,oppositec,oppositeb;
	int tic=is->totalicorners;

	count=0;
	//cycle through interface intersections.
	for(ic=0;ic<tic;ic++){
		for(c=0;c<2;c++){
			b=is->icorners[ic][c];
			if(b>=0){
				oppositec=(c+2)%4;
				oppositeb=is->icorners[ic][oppositec];
				if (oppositeb>=0){ count++; }
			}
		}
	}
	return count;
}
int**identifyDiagonalCorners(Interfaces*is,int diagcorners){
	//Identifies the blocks involved at the interface corners where two cells are diagonally adjacent.
	int count;
	int ic,c,b,oppositec,oppositeb;
	int tic=is->totalicorners;

	int **diagonalCorners=malloc(sizeof(int *)*diagcorners);
	for(c=0;c<diagcorners;c++){ diagonalCorners[c]=malloc(sizeof(int)*4); }
	count=0;
	//cycle through interface intersections.
	for(ic=0;ic<tic;ic++){
		for(c=0;c<2;c++){
			b=is->icorners[ic][c];
			if(b>=0){
				oppositec=(c+2)%4;
				oppositeb=is->icorners[ic][oppositec];
				if (oppositeb>=0){
					diagonalCorners[count][0]=b;
					diagonalCorners[count][1]=oppositeb;
					diagonalCorners[count][2]=c;
					diagonalCorners[count][3]=oppositec;
					count++;
				}
			}
		}
	}
	return diagonalCorners;
}
int checkCouch(int icorner,int c,Interfaces*is){
	//Returns -1 for no couch, otherwise the block id.
	//icorner is the id for the previously identified interface corner, i.e. index to icorners.
	//c is the corner id, i.e. 0-3.
	int couched=-1;
	int b1,b2,b3;
	b1=is->icorners[icorner][(c+3)%4];
	b2=is->icorners[icorner][c];
	b3=is->icorners[icorner][(c+1)%4];
	if((b1>=0) & (b2>=0) & (b3>=0)){
		couched=b2;
	}
	return couched;
}
void symlapBoundaryCorners(SymLap*s,Grid*g,Interfaces*is){
	//fixed/wall boundary corner enforcement in blocks that have walls at corners (but no wall sides, rather interfaces)
	//diagcorners: interface corners where there are diagonally adjacent cells.
	int ic,c,b,couch,oppositeb,oppositec;
	int mc,mr;
	int cornerci,cornerxi,oppositecornerci;
	double sign;
	int diagcorners,dc;
	int tic=is->totalicorners;

	//Process wall corners.
	for(ic=0;ic<tic;ic++){
		for(c=0;c<4;c++){
			b=is->icorners[ic][c];
			if(b>=0){
				oppositec=(c+2)%4;
				oppositeb=is->icorners[ic][oppositec];
				couch=checkCouch(ic,c,is);
				cornerci=getcornercellci(g,b,c);
				cornerxi=getcornerxi(g,b,c);
				if((couch>=0) & (oppositeb<0)){//couched and no diagonally adjacent cell.
					fprintf(is->out,"Corner procedure reached.\n");
					if(c%2==0){ sign=-1; } else { sign=1; }
					s->d0[b][cornerci]-=0.5*sign*g->B12[b][cornerxi];
				}
			}
		}
	}
	//Process interface corners with diagonal adjacency.
	diagcorners=countSymlapBoundaryCorners(is);
	fprintf(is->out,"symlap diagonal adjacency corners: %d\n",diagcorners);
	mallocCSR2(s,diagcorners);
	int**diagonalCorners=identifyDiagonalCorners(is,diagcorners);
	for(dc=0;dc<diagcorners;dc++){
		b=diagonalCorners[dc][0];
		oppositeb=diagonalCorners[dc][1];
		c=diagonalCorners[dc][2];
		oppositec=diagonalCorners[dc][3];
		cornerci=getcornercellci(g,b,c);
		oppositecornerci=getcornercellci(g,oppositeb,oppositec);
		mc=g->cstart[b]+cornerci;
		mr=g->cstart[oppositeb]+oppositecornerci;
		s->csr2.ri[dc]=mr;
		s->csr2.ci[dc]=mc;
		if(c%2==0){ sign=1; } else { sign=-1; }
		s->csr2.v[dc]=0.5*sign*g->B12[b][cornerxi];
	}
	//free memory.
	for(dc=0;dc<diagcorners;dc++){
		free(diagonalCorners[dc]);
	}
	free(diagonalCorners);
}
void symlapBoundary(SymLap*s,Grid*g,BoundaryConditions*bc,Interfaces*is){
	//fixed/wall boundaries, to be called last.
	//CURRENTLY HAS A FLAWED BOUNDARY TRAVERSAL (need adj for ci)
	//ALSO THE CONDITIONAL IS BOTCHED: should be 'or' instead of 'and'.
	int b,side,sign;
	int xi,hi,vi,ci,xd0;
	int st,iv,en;
	double bbb,met1,met2;
	for(b=0;b<g->totalblocks;b++){
		xd0=g->xdim[b][0];
		for(side=0;side<4;side++){
			if((bc->id[b][side]=='f') & (bc->id[b][side]=='w')){
				printf("Boundary procedure reached.\n");
				st=bc->side[b][side][0];
				iv=bc->side[b][side][1];
				en=bc->side[b][side][2];
				switch(side){
				case 0: //adj=0;
					sign=1;break;
				case 1: //adj=-1;
					sign=-1;break;
				case 2: //adj=-g->cdim[b][0];
					sign=-1;break;
				case 3: //adj=0;
					sign=1;break;
				}
				for(xi=st;xi<=en;xi++){
					ci=hi=xi-xi/xd0;
					vi=xi;
					switch(side%2){
					case 0: bbb=g->B22[b][hi]; break;
					case 1: bbb=g->B11[b][vi]; break;
					}
					met1=g->B12[b][xi];
					met2=g->B12[b][xi+iv];
					//additions to center cell
					s->d0[b][ci]+=0.5*sign*(met2-met1);
					s->d0[b][ci]+=-bbb;
				}
			}
		}
	}
}
void symlapInterfaceCSR(SymLap*s,Grid*g,Interfaces* is){
	int i,b1,b2,side,ip;
	int adj1,adj2;
	double bbb;
	int mc,mr;//matrix column, matrix row.
	int xi10,xi20,ci10,ci20;
	//Since Bxx will be identical on the interface for both blocks, only the first block indices are needed.
	int xi1,vi1,hi1,xiv1;
	int ci1,ci2,civ1,civ2;
	double sign;
	int entryi=0;
	fprintf(is->out,"Cell indices for each interface:\n");
	for(i=0;i<is->totalinterfaces;i++){
		b1=is->block[i];
		side=is->side[i];
		b2=is->adjacentBlock[i];
		xi10=is->blockx[i][0]-g->xstart[b1];
		xi20=is->adjacentBlockx[i][0]-g->xstart[b2];
		switch(side%2){
		case 0: civ1=civ2=xiv1=1; break;
		case 1: civ1=g->cdim[b1][0];civ2=g->cdim[b2][0];
				xiv1=g->xdim[b1][0]; break;
		}
		switch(side){
		case 0:
			adj1=0;
			adj2=-g->cdim[b2][0];
			sign=1;
			break;
		case 1:
			adj1=-1;
			adj2=0;
			sign=-1;
			break;
		case 2:
			adj1=-g->cdim[b1][0];
			adj2=0;
			sign=-1;
			break;
		case 3:
			adj1=0;
			adj2=-1;
			sign=1;
			break;
		}
		ci10=xi10-xi10/g->xdim[b1][0]+adj1;
		ci20=xi20-xi20/g->xdim[b2][0]+adj2;
		fprintf(is->out,"Interface %d:\n",i);
		for(ip=0;ip<is->xsize[i]-1;ip++){
			ci1=ci10+ip*civ1;
			ci2=ci20+ip*civ2;
			vi1=xi1=xi10+ip*xiv1;
			hi1=xi1-xi1/g->xdim[b1][0];
			fprintf(is->out,"%d\t%d\n",ci1,ci2);

			switch(side%2){
			case 0: bbb=g->B22[b1][hi1]; break;
			case 1: bbb=g->B11[b1][vi1]; break;
			}

			//column is associated with b1, row with b2.
			mc=g->cstart[b1]+ci1;
			mr=g->cstart[b2]+ci2;

			if(ip>0){
				s->csr.v[entryi]=0.5*sign*g->B12[b1][xi1];
				s->csr.ci[entryi]=mc;
				s->csr.ri[entryi]=mr-civ2;
				entryi++;
			}
			s->csr.v[entryi]=bbb;
			s->csr.ci[entryi]=mc;
			s->csr.ri[entryi]=mr;
			entryi++;
			if(ip<is->xsize[i]-2){
				s->csr.v[entryi]=-0.5*sign*g->B12[b1][xi1+xiv1];
				s->csr.ci[entryi]=mc;
				s->csr.ri[entryi]=mr+civ2;
				entryi++;
			}
		}
	}
}
void slrhs(SymLap*s,Grid*g){
	int ci,b;
	for(b=0;b<g->totalblocks;b++){
		for(ci=0;ci<g->cdim[b][3];ci++){
			s->b[g->cstart[b]+ci]=-g->Jcc[b][ci];
		}
	}
}
void slconstruct(SymLap*s,Grid*g,BoundaryConditions*bc,Interfaces*i){
	prepSymlap(s,g,i);
	symlapInterior(s,g,i);
	symlapInterfaceCSR(s,g,i);
	symlapBoundary(s,g,bc,i);
	symlapBoundaryCorners(s,g,i);
	slrhs(s,g);
	initv(s->x,s->dim,0);
}

#endif /* MATFREE_H_ */
