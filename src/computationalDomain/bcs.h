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
void segregatedBC(BoundaryConditions*bc,int totalblocks){
	int b,s,wi=0,fi=0,ii=0,si=0,zi=0;
	for(b=0;b<totalblocks;b++){
		for(s=0;s<4;s++){
			switch(bc->bcid[b].id[s]){
			case 'w':
				bc->w[wi].b=b;
				bc->w[wi].s=s;
				wi++;
				break;
			case 'f':
				bc->f[fi].b=b;
				bc->f[fi].s=s;
				fi++;
				break;
			case 'i':
				if(bc->bid[b].id[s]>b){//If the connecting block is greater than the current block.
					bc->i[ii].b=b;
					bc->i[ii].s=s;
					ii++;
				}
				break;
			case 'z':
				bc->z[zi].b=b;
				bc->z[zi].s=s;
				zi++;
				break;
			case 's':
				bc->s[si].b=b;
				bc->s[si].s=s;
				si++;
				break;
			}
		}
	}
	if(bc->wn>0){
		fprintf(bc->out,"Walls (block side):\n");
		for(wi=0;wi<bc->wn;wi++){
			fprintf(bc->out,"\t%d\t%d\n",bc->w[wi].b,bc->w[wi].s);
		}
	}
	if(bc->fn>0){
		fprintf(bc->out,"Fixed (block side):\n");
		for(fi=0;fi<bc->fn;fi++){
			fprintf(bc->out,"\t%d\t%d\n",bc->f[fi].b,bc->f[fi].s);
		}
	}
	if(bc->in>0){
		fprintf(bc->out,"Interfaces (block side):\n");
		for(ii=0;ii<bc->in;ii++){
			fprintf(bc->out,"\t%d\t%d\n",bc->i[ii].b,bc->i[ii].s);
		}
	}
	if(bc->zn>0){
		fprintf(bc->out,"Zero Gradients (block side):\n");
		for(zi=0;zi<bc->zn;zi++){
			fprintf(bc->out,"\t%d\t%d\n",bc->z[zi].b,bc->z[zi].s);
		}
	}
	if(bc->sn>0){
		fprintf(bc->out,"Sliding Interfaces (block side):\n");
		for(si=0;si<bc->sn;si++){
			fprintf(bc->out,"\t%d\t%d",bc->s[si].b,bc->s[si].s);
		}
	}
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
void countAndIdentifyCommonCorners(BoundaryConditions*bc,int totalblocks){
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
	bc->wcn=wcn;
	bc->fcn=fcn;
	bc->icn=icn;
	mallocBC3(bc,bc->wcn,bc->fcn,bc->icn);
	int c;
	if(bc->icn>0){
		fprintf(bc->out,"Interface Corners:\n");
		for(c=0;c<icn;c++){
			fprintf(bc->out,"\t%d\t%d\t%d\t%d\n",ic[c][0],ic[c][1],ic[c][2],ic[c][3]);
			for(b=0;b<4;b++){ bc->ic[c].block[b]=ic[c][b]; }
		}
	}
	if(bc->fcn>0){
		fprintf(bc->out,"Fixed Corners:\n");
		for(c=0;c<fcn;c++){
			fprintf(bc->out,"\t%d\t%d\t%d\t%d\n",fc[c][0],fc[c][1],fc[c][2],fc[c][3]);
			for(b=0;b<4;b++){ bc->fc[c].block[b]=fc[c][b]; }
		}
	}
	if(bc->wcn>0){
		fprintf(bc->out,"Wall Corners:\n");
		for(c=0;c<wcn;c++){
			fprintf(bc->out,"\t%d\t%d\t%d\t%d\n",wc[c][0],wc[c][1],wc[c][2],wc[c][3]);
			for(b=0;b<4;b++){ bc->wc[c].block[b]=wc[c][b]; }
		}
	}
}
void fillBCs(BoundaryConditions*bc,Dimension*d,Grid*g){
	mallocBC1(bc,d->tb);
	countAndIdentifyBC(bc);
	mallocBC2(bc,bc->wn,bc->fn,bc->in,bc->zn,bc->sn);
	blockConnectivity(bc,d,g);//need this to sort out interfaces in segregatedBC and identifying common corners.
	segregatedBC(bc,d->tb);
	countAndIdentifyCommonCorners(bc,d->tb);
}
