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
void incrementBC(char bcchar,BoundaryConditions*bc){
	switch(bcchar){
	case 'w': bc->w.n++; break;
	case 'f': bc->f.n++; break;
	case 'i': bc->i.n++; break;
	case 'z': bc->z.n++; break;
	case 's': bc->s.n++; break;
	}
}
void countBC(BoundaryConditions*bc,char*inputfile){
	int b,success,maxchars=200;
	char bc0,bc1,bc2,bc3,line[maxchars],filename[maxchars];
	bc->w.n=0;
	bc->f.n=0;
	bc->i.n=0;
	bc->z.n=0;
	bc->s.n=0;
	//Read blockDict.
	FILE * file = fopen (inputfile, "rt");
	while(fgets(line, maxchars, file) != NULL){
		success = sscanf (line, "%d %s %c %c %c %c",&b,filename,&bc0,&bc1,&bc2,&bc3);
		if(line[0]=='#'){ break; }
		if(success>0){
			incrementBC(bc0,bc);
			incrementBC(bc1,bc);
			incrementBC(bc2,bc);
			incrementBC(bc3,bc);
		}
	}
	//Interface requires further processing.
	//if((bc->i.n%2)!=0){ printf("Error: Individual interface declarations not even!\n"); }
	//else { bc->i.n/=2; }
	fprintf(bc->out,"Boundary condition tally:\n");
	fprintf(bc->out,"Walls: %d\n",bc->w.n);
	fprintf(bc->out,"Fixed Velocity: %d\n",bc->f.n);
	fprintf(bc->out,"Interfaces: %d\n",bc->i.n);
	fprintf(bc->out,"Zero Gradient: %d\n",bc->z.n);
	fprintf(bc->out,"Sliding Interfaces: %d\n",bc->s.n);
	fclose(file);
}
void fillBCSub(char bcid,int b,int s,BoundaryConditions*bc,int*counter){
	switch(bcid){
	case 'w':
		bc->w.id[counter[0]].b=b;
		bc->w.id[counter[0]].s=s;
		counter[0]++;
		break;
	case 'f':
		bc->f.id[counter[1]].b=b;
		bc->f.id[counter[1]].s=s;
		counter[1]++;
		break;
	case 'i':
		bc->i.id[counter[2]].b=b;
		bc->i.id[counter[2]].s=s;
		counter[2]++;
		break;
	case 'z':
		bc->z.id[counter[3]].b=b;
		bc->z.id[counter[3]].s=s;
		counter[3]++;
		break;
	case 's':
		bc->s.id[counter[4]].b=b;
		bc->s.id[counter[4]].s=s;
		counter[4]++;
		break;
	}
}
void fillBC(BoundaryConditions*bc,char*inputfile){
	int b,success,maxchars=200,counter[5];
	char bc0,bc1,bc2,bc3,line[maxchars],filename[maxchars];
	countBC(bc,inputfile);
	bc->w.id=malloc(bc->w.n*sizeof(BoundaryID));
	bc->f.id=malloc(bc->f.n*sizeof(BoundaryID));
	bc->i.id=malloc(bc->i.n*sizeof(BoundaryID));
	bc->z.id=malloc(bc->z.n*sizeof(BoundaryID));
	bc->s.id=malloc(bc->s.n*sizeof(BoundaryID));
	//Read blockDict.
	counter[0]=counter[1]=counter[2]=counter[3]=counter[4]=0;
	FILE * file = fopen (inputfile, "rt");
	while(fgets(line, maxchars, file) != NULL){
		success = sscanf (line, "%d %s %c %c %c %c",&b,filename,&bc0,&bc1,&bc2,&bc3);
		if(line[0]=='#'){ break; }
		if(success>0){

		}
	}
	if((bc->i.n%2)!=0){ printf("Error: Individual interface declarations not even!\n"); }
	else { bc->i.n/=2; }
	fprintf(bc->out,"Boundary condition tally:\n");
	fprintf(bc->out,"Walls: %d\n",bc->w.n);
	fprintf(bc->out,"Fixed Velocity: %d\n",bc->f.n);
	fprintf(bc->out,"Interfaces: %d\n",bc->i.n);
	fprintf(bc->out,"Zero Gradient: %d\n",bc->z.n);
	fprintf(bc->out,"Sliding Interfaces: %d\n",bc->s.n);
	fclose(file);
}
void checkCornerSides(int c,int b,BoundaryConditions*bc){
	//int s1,s2;
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
