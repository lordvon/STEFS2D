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
	fprintf(bc->out,"Walls: %d\n",bc->wn);
	fprintf(bc->out,"Fixed Velocity: %d\n",bc->fn);
	fprintf(bc->out,"Interfaces: %d\n",bc->in);
	fprintf(bc->out,"Zero Gradient: %d\n",bc->zn);
	fprintf(bc->out,"Sliding Interfaces: %d\n",bc->sn);
	fclose(file);
}
int adjacentBlockCheck(Grid* grid,int xi1,int xi2,int b1,int b2){
	double diff=sqrt(
			pow(grid->x[b1][xi1]-grid->x[b2][xi2],2)+
			pow(grid->y[b1][xi1]-grid->y[b2][xi2],2));
	int match=0; if(diff<1e-6){ match=1; }
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
		if(adjacentBlockCheck(g,xi1,xi2,b1,b2)){ matchedblock=b2; break; }
	}
	return matchedblock;
}
void blockConnectivity(BoundaryConditions*bc,Dimension*d,Grid*g){
	int b,s;
	for(b=0;b<d->tb;b++){
		for(s=0;s<4;s++){
			bc->bid[b].id[s]=getAdjacentBlock(b,s,d,g);
		}
	}
}
void countCommonCorners(BoundaryConditions*bc,int totalblocks){
	int wcn=0,fcn=0,icn=0;
	int cornerid[4];
	int b,s1;
	for(b=0;b<totalblocks;b++){
		for(s1=0;s1<4;s1++){
			if(bc->bcid[b].id[s1]=='i'){

			}
		}
	}
}
void fillBCs(BoundaryConditions*bc,Dimension*d,Grid*g){
	mallocBC1(bc,d->tb);
	countAndIdentifyBC(bc);
	blockConnectivity(bc,d,g);
	mallocBC2(bc,bc->wn,bc->fn,bc->in,bc->zn,bc->sn);
	countCommonCorners(bc,d->tb);
	mallocBC3(bc,bc->wcn,bc->fcn,bc->icn);
}

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
