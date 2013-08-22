void countAndIdentifyBC(BoundaryConditions*bc){
	int b,success,maxchars=200,i;
	char bcid[4],line[maxchars],filename[maxchars];
	bc->wn=0;
	bc->fn=0;
	bc->in=0;
	bc->zn=0;
	bc->sn=0;
	//Read blockDict.
	FILE * file = fopen ("in/blockDict", "rt");
	while(fgets(line, maxchars, file) != NULL){
		success = sscanf (line, "%d %s %c %c %c %c",&b,filename,
				&bcid[0],&bcid[1],&bcid[2],&bcid[3]);
		if(line[0]=='#'){ break; }
		if(success>0){
			for(i=0;i<4;i++){
				bc->bcid[b].id[i]=bcid[i];
				switch(bcid[i]){
				case 'w': bc->wn++; break;
				case 'f': bc->fn++; break;
				case 'i': bc->in++; break;
				case 'z': bc->zn++; break;
				case 's': bc->sn++; break;
				}
			}
		}
	}
	if((bc->in%2)!=0){ printf("Error: Individual interface declarations not even!\n"); }
	else { bc->in/=2; }
	fprintf(bc->out,"Boundary condition tally:\n");
	fprintf(bc->out,"\tWalls: %d\n",bc->wn);
	fprintf(bc->out,"\tFixed Velocity: %d\n",bc->fn);
	fprintf(bc->out,"\tInterfaces: %d\n",bc->in);
	fprintf(bc->out,"\tZero Gradient: %d\n",bc->zn);
	fprintf(bc->out,"\tSliding Interfaces: %d\n",bc->sn);
	fclose(file);
}
int adjacentBlockCheck(Grid* grid,int xi1,int xi2,int b1,int b2){
	double diff=sqrt(
			pow(grid->x[b1][xi1]-grid->x[b2][xi2],2)+
			pow(grid->y[b1][xi1]-grid->y[b2][xi2],2));
	int match=0; if(diff<(1.0e-5)){ match=1; }
	return match;
}
int getAdjacentBlock(int b1,int s1,Dimension*d,Grid*g){
	//from block and side, deduce adjacent block. Returns -1 if no block.
	int b2,s2,matchedblock;
	int xi1,xi2;
	xi1=d->xside[b1].b[s1].st;
	//search through all other blocks to find interfaced block.
	matchedblock=-1;
	s2=(s1+2)%4;
	for(b2=0;b2<d->tb;b2++){
		xi2=d->xside[b2].b[s2].st;
		if(adjacentBlockCheck(g,xi1,xi2,b1,b2)){
			matchedblock=b2;
			break;
		}
	}
	return matchedblock;
}
void blockConnectivity(BoundaryConditions*bc,Dimension*d,Grid*g){
	int b,s;
	fprintf(bc->out,"Block connectivities:\n");
	for(b=0;b<d->tb;b++){
		for(s=0;s<4;s++){
			if(bc->bcid[b].id[s]=='i'){
				bc->bid[b].id[s]=getAdjacentBlock(b,s,d,g);
			} else {
				bc->bid[b].id[s]=-1;
			}
		}
		fprintf(bc->out,"\t%d\t%d\t%d\t%d\n",bc->bid[b].id[0],bc->bid[b].id[1],bc->bid[b].id[2],bc->bid[b].id[3]);
	}
}
int checkCommonCornerRepeat(int**list,int*corner,int nlist){
	int i,c,match,repeated=0;
	for(i=0;i<nlist;i++){
		match=1;
		for(c=0;c<4;c++){
			if(list[i][c]!=corner[c]){ match=0; break; }
		}
		if(match){ repeated=1; break; }
	}
	return repeated;
}
void fillCorners(int*corner1,int*corner2,int b1,int s1,
		int*bid1,int*bid2){
	int b2,s2;
	int adj11,adj12,adj21,adj22;
	int c11,c12;
	b2=bid1[s1];
	s2=(s1+2)%4;
	//end 1 is clockwise from end 2.
	//end 1
	c11=s1;
	adj11=(s1+3)%4;
	adj22=(s2+1)%4;
	corner1[(c11+3)%4]=b2;
	corner1[c11]=b1;
	corner1[(c11+1)%4]=bid1[adj11];
	corner1[(c11+2)%4]=bid2[adj22];
	//end 2
	c12=(s1+1)%4;
	adj12=(s1+1)%4;
	adj21=(s2+3)%4;
	corner2[(c12+3)%4]=bid1[adj12];
	corner2[c12]=b1;
	corner2[(c12+1)%4]=b2;
	corner2[(c12+2)%4]=bid2[adj21];
}
int checkCommonCorner(int*corner){
	int iscc,c,b=0;
	for(c=0;c<4;c++){ if(corner[c]>=0){ b++; } }
	if(b>2){ iscc=1; } else { iscc=0; }
	return iscc;
}
int checkCornerType(char id,int*corner,BoundaryConditions*bc){
	int istype=0,c,s1,s2;
	for(c=0;c<4;c++){
		if(corner[c]>=0){
			s1=c;
			s2=(s1+3)%4;
			if((bc->bcid[corner[c]].id[s1]==id) | (bc->bcid[corner[c]].id[s2]==id)){
				istype=1;
				break;
			}
		}
	}
	return istype;
}
int**addCornerToList(int*corner,int**list,int oldamount){
	int**newlist=malloc(sizeof(int*)*(oldamount+1));
	int i;
	for(i=0;i<oldamount;i++){
		newlist[i]=list[i];
	}
	free(list);
	newlist[oldamount]=malloc(sizeof(int)*4);
	for(i=0;i<4;i++){
		newlist[oldamount][i]=corner[i];
	}
	return newlist;
}
void countCommonCorners(BoundaryConditions*bc,int totalblocks){
	int wcn=0,fcn=0,icn=0;
	int **wc=NULL,**fc=NULL,**ic=NULL;
	int corner1[4],corner2[4];//represents blocks adjacent to corner; indexed by corner id.
	int b,s,b2;
	for(b=0;b<totalblocks;b++){
		for(s=0;s<4;s++){
			if(bc->bcid[b].id[s]=='i'){
				b2=bc->bid[b].id[s];
				fillCorners(corner1,corner2,b,s,bc->bid[b].id,bc->bid[b2].id);
				if(checkCommonCorner(corner1)){
					if(checkCornerType('w',corner1,bc)){
						if(!checkCommonCornerRepeat(wc,corner1,wcn)){
							wc=addCornerToList(corner1,wc,wcn);wcn++;
						}
					} else if(checkCornerType('f',corner1,bc)){
						if(!checkCommonCornerRepeat(fc,corner1,fcn)){
							fc=addCornerToList(corner1,fc,fcn);fcn++;
						}
					} else {
						if(!checkCommonCornerRepeat(fc,corner1,fcn)){
							ic=addCornerToList(corner1,ic,icn);icn++;
						}
					}
				}
				if(checkCommonCorner(corner2)){
					if(checkCornerType('w',corner2,bc)){
						if(!checkCommonCornerRepeat(wc,corner2,wcn)){
							wc=addCornerToList(corner2,wc,wcn);wcn++;
						}
					} else if(checkCornerType('f',corner2,bc)){
						if(!checkCommonCornerRepeat(fc,corner2,fcn)){
							fc=addCornerToList(corner2,fc,fcn);fcn++;
						}
					} else {
						if(!checkCommonCornerRepeat(ic,corner2,icn)){
							ic=addCornerToList(corner2,ic,icn);icn++;
						}
					}
				}
			}
		}
	}
	int c;
	fprintf(bc->out,"Interface Corners:\n");
	for(c=0;c<icn;c++){
		fprintf(bc->out,"\t%d\t%d\t%d\t%d\n",ic[c][0],ic[c][1],ic[c][2],ic[c][3]);
	}
	fprintf(bc->out,"Fixed Corners:\n");
	for(c=0;c<fcn;c++){
		fprintf(bc->out,"\t%d\t%d\t%d\t%d\n",fc[c][0],fc[c][1],fc[c][2],fc[c][3]);
	}
	fprintf(bc->out,"Wall Corners:\n");
	for(c=0;c<wcn;c++){
		fprintf(bc->out,"\t%d\t%d\t%d\t%d\n",wc[c][0],wc[c][1],wc[c][2],wc[c][3]);
	}
}
void fillBCs(BoundaryConditions*bc,Dimension*d,Grid*g){
	mallocBC1(bc,d->tb);
	countAndIdentifyBC(bc);
	mallocBC2(bc,bc->wn,bc->fn,bc->in,bc->zn,bc->sn);
	blockConnectivity(bc,d,g);
	countCommonCorners(bc,d->tb);
	mallocBC3(bc,bc->wcn,bc->fcn,bc->icn);
}

/*
void fillBCSub(char bcid,int b,int s,BoundaryConditions*bc,int*counter){
	switch(bcid){
	case 'w':
		bc->w[counter[0]].b=b;
		bc->w[counter[0]].s=s;
		counter[0]++;
		break;
	case 'f':
		bc->f[counter[1]].b=b;
		bc->f[counter[1]].s=s;
		counter[1]++;
		break;
	case 'i':
		bc->i[counter[2]].b=b;
		bc->i[counter[2]].s=s;
		counter[2]++;
		break;
	case 'z':
		bc->z[counter[3]].b=b;
		bc->z[counter[3]].s=s;
		counter[3]++;
		break;
	case 's':
		bc->s[counter[4]].b=b;
		bc->s[counter[4]].s=s;
		counter[4]++;
		break;
	}
}
void fillBC(BoundaryConditions*bc,char*inputfile,int totalblocks){
	int b,success,maxchars=200,counter[5];
	char bc0,bc1,bc2,bc3,line[maxchars],filename[maxchars];
	countBC(bc);
	mallocBC1(bc,totalblocks);
	mallocBC2(bc,bc->wn,bc->fn,bc->in,bc->zn,bc->sn);
	mallocBC3(bc,bc->wcn,bc->fcn,bc->icn);
	//Read blockDict.
	counter[0]=counter[1]=counter[2]=counter[3]=counter[4]=0;
	FILE * file = fopen (inputfile, "rt");
	while(fgets(line, maxchars, file) != NULL){
		success = sscanf (line, "%d %s %c %c %c %c",&b,filename,&bc0,&bc1,&bc2,&bc3);
		if(line[0]=='#'){ break; }
		if(success>0){

		}
	}
	if((bc->in%2)!=0){ printf("Error: Individual interface declarations not even!\n"); }
	else { bc->in/=2; }
	fprintf(bc->out,"Boundary condition tally:\n");
	fprintf(bc->out,"Walls: %d\n",bc->wn);
	fprintf(bc->out,"Fixed Velocity: %d\n",bc->fn);
	fprintf(bc->out,"Interfaces: %d\n",bc->in);
	fprintf(bc->out,"Zero Gradient: %d\n",bc->zn);
	fprintf(bc->out,"Sliding Interfaces: %d\n",bc->sn);
	fclose(file);
}
void checkCornerSides(int c,int b,BoundaryConditions*bc){
	//int s1,s2;
}
int getPrecedence(char bc){
	//The higher the more imposing.
	int p=0;
	switch(bc){
	case 'w': p=5; break;
	case 'f': p=4; break;
	case 'i': p=3; break;
	case 'z': p=2; break;
	case 's': p=1; break;
	default: printf("BC not recognized!\n"); break;
	}
	return p;
}
void countIC(Interfaces*is,BoundaryConditions*bc){
	//counts interface corners.
	int i,c;
	for(i=0;i<is->totalicorners;i++){
		for(c=0;c<4;c++){
			if(is->icorners[i][c]>-1){
			}
		}
	}
}
*/
