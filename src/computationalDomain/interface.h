/*
 * interface.h
 *
 *  Created on: Jul 15, 2013
 *      Author: lordvon
 */

#ifndef INTERFACE_H_
#define INTERFACE_H_

int blockMatch(int xi1,int xi2,int b1,int b2,
		Grid* grid, Interfaces* interfaces){
	double diff;
	int match=0;
	diff=sqrt(pow(grid->x[b1][xi1]-grid->x[b2][xi2],2)+pow(grid->y[b1][xi1]-grid->y[b2][xi2],2));
	if(diff<(interfaces->matchtolerance)){
		match=1;
	}
	return match;
}
int getSideXdim(int side, int *xdim){
	int sidedim;
	if((side==0) | (side==2)){
		sidedim=xdim[0];
	} else {
		sidedim=xdim[1];
	}
	return sidedim;
}
int getInterfaceBlock(int block,int side,
		Grid* grid,Interfaces* interfaces){
	//from block and side, deduce adjacent block.
	int otherblock,matchedblock;
	int xi1,xi2;
	if(side==0){ xi1=0; }
	else if(side==1){ xi1=grid->xdim[block][0]-1; }
	else if(side==2){ xi1=grid->xdim[block][3]-grid->xdim[block][0]; }
	else if(side==3){ xi1=0; }
	matchedblock=-1;
	//search through all other blocks to find interfaced block.
	for(otherblock=0;otherblock<grid->totalblocks;otherblock++){
		if(side==0){ xi2=grid->xdim[otherblock][3]-grid->xdim[otherblock][0]; }
		else if(side==1){ xi2=0; }
		else if(side==2){ xi2=0; }
		else if(side==3){ xi2=grid->xdim[otherblock][0]-1; }
		if(blockMatch(xi1,xi2,block,otherblock,grid,interfaces)>0){
			matchedblock=otherblock;
			break;
		}
	}
	if(matchedblock==-1){
		printf("Error: Interfacing block not found at for block %d (side %d)!\n",block,side);
	} else if(0){
		printf("Matched block for block %d (side %d): %d\n",block,side,matchedblock);
	}
	return matchedblock;
}
void establishInterfaces(Grid* grid,Interfaces* interfaces,BoundaryConditions* bc){
	int b,otherb,s,i,n,othern,row;
	i=0;row=grid->totaledges;
	interfaces->totalinterfacepoints=0;
	//iterate through all blocks.
	for(b=0;b<grid->totalblocks;b++){
		//iterate through all sides to find declared interfaces.
		for(s=0;s<4;s++){
			if(bc->id[b][s]=='i'){
				otherb=getInterfaceBlock(b,s,grid,interfaces);
				if(otherb>=b){
					//add to interfaces
					interfaces->block[i]=b;
					interfaces->side[i]=s;
					interfaces->adjacentBlock[i]=otherb;
					n=getSideXdim(s,grid->xdim[b]);
					othern=getSideXdim(s,grid->xdim[b]);
					if(n!=othern){ printf("Error: Interface dimension does not match (interface %d)!\n",i); }
					interfaces->xsize[i]=n;
					interfaces->totalinterfacepoints+=n;
					//printf("interfacepoints: %d\n",n);
					interfaces->rowstart[i]=row;
					i++;
					row+=n;
				}
			}
		}
	}
}
void fillInterfaceCoords(Grid* g,Interfaces* is){
	int i,ip,iv1,iv2,s1,s2,b1,b2,side;
	for(i=0;i<is->totalinterfaces;i++){
		b1=is->block[i];
		b2=is->adjacentBlock[i];
		side=is->side[i];
		is->blockx[i]=malloc(is->xsize[i]*sizeof(int));
		is->adjacentBlockx[i]=malloc(is->xsize[i]*sizeof(int));
		if(side==0){
			iv1=iv2=1;
			s1=0;
			s2=g->xdim[b2][3]-g->xdim[b2][0];
		} else if(side==1){
			iv1=g->xdim[b1][0];
			iv2=g->xdim[b2][0];
			s1=g->xdim[b1][0]-1;
			s2=0;
		} else if(side==2){
			iv1=iv2=1;
			s1=g->xdim[b1][3]-g->xdim[b1][0];
			s2=0;
		} else if(side==3){
			iv1=g->xdim[b1][0];
			iv2=g->xdim[b2][0];
			s1=0;
			s2=g->xdim[b2][0]-1;
		}
		for(ip=0;ip<is->xsize[i];ip++){
			//interfaces->blockx[i][ip]=s1+ip*iv1;
			//interfaces->adjacentBlockx[i][ip]=s2+ip*iv2;
			is->blockx[i][ip]=g->xstart[b1]+s1+ip*iv1;
			is->adjacentBlockx[i][ip]=g->xstart[b2]+s2+ip*iv2;
		}
	}
}
void fillInterfaceEdges(Grid* g,Interfaces* is){
	int i,ie,iv1,iv2,s1,s2,b1,b2,side;
	int es;
	for(i=0;i<is->totalinterfaces;i++){
		b1=is->block[i];
		b2=is->adjacentBlock[i];
		side=is->side[i];
		es=is->xsize[i]-1;
		is->e[i]=malloc(es*sizeof(int));
		is->ae[i]=malloc(es*sizeof(int));
		if(side==0){
			iv1=iv2=1;
			s1=g->vdim[b1][3];
			s2=g->vdim[b2][3]+g->hdim[b2][3]-g->hdim[b2][0];
		} else if(side==1){
			iv1=g->vdim[b1][0];
			iv2=g->vdim[b2][0];
			s1=g->vdim[b1][0]-1;
			s2=0;
		} else if(side==2){
			iv1=iv2=1;
			s1=g->vdim[b1][3]+g->hdim[b1][3]-g->hdim[b1][0];
			s2=g->vdim[b2][3];
		} else if(side==3){
			iv1=g->vdim[b1][0];
			iv2=g->vdim[b2][0];
			s1=0;
			s2=g->vdim[b2][0]-1;
		}
		for(ie=0;ie<es;ie++){
			//interfaces->blockx[i][ip]=s1+ip*iv1;
			//interfaces->adjacentBlockx[i][ip]=s2+ip*iv2;
			is->e[i][ie]=g->estart[b1]+s1+ie*iv1;
			is->ae[i][ie]=g->estart[b2]+s2+ie*iv2;
		}
	}
}
void fillConnectivity(Grid* g,Interfaces* is){
	int b,i,b1,b2,s1,s2;
	is->conn=malloc(g->totalblocks*sizeof(int*));
	for(b=0;b<g->totalblocks;b++){
		is->conn[b]=malloc(4*sizeof(int));
	}
	for(i=0;i<is->totalinterfaces;i++){
		b1=is->block[i];
		b2=is->adjacentBlock[i];
		s1=is->side[i];
		if((s1-2)<0){ s2=s1+2; } else { s2=s1-2; }
		is->conn[b1][s1]=b2;
		is->conn[b2][s2]=b1;
	}
}
void addcornercell(Interfaces*is,int icornerid,
		int block,int side1,int side2){
	int cornerid=getcornerid(side1,side2);
	is->icorners[icornerid][cornerid]=block;
}
int processInterfaceEnd(Interfaces*is,int interface,
		int end,int icornerid,
		int***ifromb,int *ichecklist){
	//end: 0 or 1. Corresponds to each end of the interface line.
	int b1,b2,s1,s2,adjacentSide,oppositeAdjacentSide;
	int i1,i2,icheck1,icheck2,is1,is2;
	b1=is->block[interface];
	b2=is->adjacentBlock[interface];
	s1=is->side[interface];
	s2=(s1+2)%4;
	adjacentSide=(s1+1+2*end)%4;//adjacent side to interface (same number for both blocks)
	oppositeAdjacentSide=(adjacentSide+2)%4;
	i1=ifromb[b1][adjacentSide][0];
	i2=ifromb[b2][adjacentSide][0];
	is1=i1 >= 0;
	is2=i2 >= 0;
	icheck1=icheck2=1;
	if(is1){ icheck1=ichecklist[i1]<1; }
	if(is2){ icheck2=ichecklist[i2]<1; }
	//icheck1=(i1 >= 0) && (ichecklist[i1]<1);
	//icheck2=(i2 >= 0) && (ichecklist[i2]<1);
	if((is1 | is2) & icheck1 & icheck2){
		//printf("An interface corner is processed.\n");
		//if(is1){ ichecklist[i1]=1; }
		//if(is2){ ichecklist[i2]=1; }
		//corners from both sides of primary interface.
		addcornercell(is,icornerid,b1,s1,adjacentSide);
		addcornercell(is,icornerid,b2,s2,adjacentSide);
		//additional corners from additional interfaces.
		if(is1 & icheck1){ addcornercell(is,icornerid,
				is->conn[b1][adjacentSide],s1,oppositeAdjacentSide); }
		if(is2 & icheck2){ addcornercell(is,icornerid,
				is->conn[b2][adjacentSide],s2,oppositeAdjacentSide); }
		icornerid++;
	}
	return icornerid;
}
int processInterface(Interfaces*is,int interface,
		int icornerid,int***ifromb,int *ichecklist){
	icornerid=processInterfaceEnd(is,interface,0,
			icornerid,ifromb,ichecklist);
	icornerid=processInterfaceEnd(is,interface,1,
			icornerid,ifromb,ichecklist);
	ichecklist[interface]=1;
	return icornerid;
}
int***makeifromb(Interfaces*is,Grid*g){
	int i,j,b,b1,b2,s1,s2;
	int ti=is->totalinterfaces;
	int tb=g->totalblocks;
	int ***ifromb=malloc(sizeof(int **)*tb);//[b]{per side}{interface, other b}
	for(b=0;b<tb;b++){
		ifromb[b]=malloc(sizeof(int *)*4);
		for(j=0;j<4;j++){
			ifromb[b][j]=malloc(sizeof(int)*2);
			ifromb[b][j][0]=ifromb[b][j][1]=-1;//will be initialized to id value if interface present.
		}
	}
	for(i=0;i<ti;i++){
		b1=is->block[i];
		b2=is->adjacentBlock[i];
		s1=is->side[i];
		s2=(s1+2)%4;
		ifromb[b1][s1][0]=i;
		ifromb[b2][s2][0]=i;
		ifromb[b1][s1][1]=b2;
		ifromb[b2][s2][1]=b1;
	}
	/*
	for(b=0;b<tb;b++){
		printf("Block %d :::\n",b);
		for(j=0;j<4;j++){
			printf("\ti: %d, b: %d\n",ifromb[b][j][0],ifromb[b][j][1]);
		}
	}
	*/
	return ifromb;
}
int counticornerssub(Interfaces* is,int i,int cornercount,int end,
		int***ifromb,int*ichecklist){
	//end: 0 or 1. Corresponds to each end of the interface line.
	int b1,b2,adjacentSide,s1,i1,i2,icheck1,icheck2,is1,is2;
	b1=is->block[i];
	b2=is->adjacentBlock[i];
	s1=is->side[i];
	adjacentSide=(s1+1+2*end)%4;//adjacent side to interface (same number for both blocks)
	i1=ifromb[b1][adjacentSide][0];
	i2=ifromb[b2][adjacentSide][0];
	//printf("i: %d, i1: %d, i2: %d\n",i,i1,i2);
	//printf("ichecklist: %d, %d, %d, %d\n",ichecklist[0],ichecklist[1],ichecklist[2],ichecklist[3]);
	//check if the adjacent sides are interfaces.
	is1=i1>=0;
	is2=i2>=0;
	//check to see if they have been checked off.
	icheck1=icheck2=1;
	if(is1){ icheck1=ichecklist[i1]<1; }
	if(is2){ icheck2=ichecklist[i2]<1; }
	if((is1 | is2) & icheck1 & icheck2){//at least one is an interface, and both are not checked off.
		cornercount++;
	}
	return cornercount;
}
int counticorners(Interfaces*is,int***ifromb){
	int i;
	int cornercount=0;
	int ti=is->totalinterfaces;
	int*ichecklist=malloc(sizeof(int)*ti);//boolean
	for(i=0;i<ti;i++){ ichecklist[i]=0; }//initialize check list.
	for(i=0;i<ti;i++){
		//end 1
		cornercount=counticornerssub(is,i,cornercount,1,ifromb,ichecklist);
		//end 0
		cornercount=counticornerssub(is,i,cornercount,0,ifromb,ichecklist);
		//finish this interface
		ichecklist[i]=1;
	}
	free(ichecklist);
	return cornercount;
}
void fillicorners(Interfaces*is,Grid*g){
	int i,icornerid,ic,c;
	int ti=is->totalinterfaces;
	int tic=is->totalicorners;

	int ***ifromb=makeifromb(is,g);
	tic=is->totalicorners=counticorners(is,ifromb);
	//fprintf(is->out,"Interface corners counted: %d\n",tic);
	if(tic>0){
		int*ichecklist=malloc(sizeof(int)*ti);//boolean
		for(i=0;i<ti;i++){ ichecklist[i]=0; }//zero the checklist.
		//malloc
		is->icorners=malloc(sizeof(int *)*tic);
		for(ic=0;ic<tic;ic++){
			is->icorners[ic]=malloc(sizeof(int)*4);
			for(c=0;c<4;c++){
				is->icorners[ic][c]=-1;//these values will change to real block ids, if available interpolation is present.
			}
		}
		//fill
		icornerid=0;
		for(i=0;i<ti;i++){
			//printf("tic\n");
			icornerid=processInterface(is,i,icornerid,ifromb,ichecklist);
		}
		free(ichecklist);
	}
	int b,j;
	for(b=0;b<g->totalblocks;b++){
		for(j=0;j<4;j++){
			free(ifromb[b][j]);
		}
		free(ifromb[b]);
	}
	free(ifromb);
}

#endif /* INTERFACE_H_ */
