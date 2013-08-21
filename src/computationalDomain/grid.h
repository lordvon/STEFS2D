void determineTotalBlocks(Grid* g,char *inputfile){
	int maxchars=200;
	int blocknumber,totalblocks;
	char bc0,bc1,bc2,bc3;
	char line[maxchars];
	char filename[maxchars];
	FILE * file = fopen (inputfile, "rt");
	int success;
	//Read total number of blocks.
	totalblocks=0;
	while(fgets(line, maxchars, file) != NULL){
		success = sscanf (line, "%d %s %c %c %c %c",&blocknumber,filename,
				&bc0,&bc1,&bc2,&bc3);
		//printf("line: %c %c %c %c\n",line[0],line[1],line[2],line[3]);
		if(line[0]=='#'){ break; }
		if(success>0){ totalblocks++; }
	}
	g->totalblocks=totalblocks;
	fclose(file);
}
void readBlockDict(Grid* grid,Interfaces* interfaces,BoundaryConditions* bc,
		char *inputfile){
	int maxchars=200;
	int blocknumber,interfacecount,i;
	char bc0,bc1,bc2,bc3;
	char line[maxchars];
	char filename[maxchars];
	double fixeducart,fixedvcart;
	FILE * file = fopen (inputfile, "rt");
	int success;
	//Read total number of blocks.
	determineTotalBlocks(grid,inputfile);
	grid->blockNames=charmalloc(grid->totalblocks,200);
	bc->id=malloc(grid->totalblocks*sizeof(char*));
	bc->fv=malloc(grid->totalblocks*sizeof(double *));
	bc->side=malloc(grid->totalblocks*sizeof(int **));
	int b,s;
	for(b=0;b<grid->totalblocks;b++){
		bc->id[b]=malloc(4*sizeof(char));
		bc->fv[b]=malloc(2*sizeof(double));
		bc->side[b]=malloc(4*sizeof(int *));
		for(s=0;s<4;s++){
			bc->side[b][s]=malloc(4*sizeof(int));
		}
	}
	//Read file names for each block.
	interfacecount=0;
	while(fgets(line, maxchars, file) != NULL){
		success = sscanf (line, "%d %s %c %c %c %c %lf %lf",&blocknumber,filename,
				&bc0,&bc1,&bc2,&bc3,&fixeducart,&fixedvcart);
		if(line[0]=='#'){ break; }
		if(success>0){
			//printf("%s\n",filename);
			strcpy(grid->blockNames[blocknumber],filename);
					//strcat(blockdirectory,strcat(filename,blockextension)));
			bc->id[blocknumber][0]=bc0;
			bc->id[blocknumber][1]=bc1;
			bc->id[blocknumber][2]=bc2;
			bc->id[blocknumber][3]=bc3;
			for(i=0;i<4;i++){
				if(bc->id[blocknumber][i]=='i'){ interfacecount++; }
			}
			if((bc0=='f') | (bc1=='f') | (bc2=='f') | (bc3=='f')){
				bc->fv[blocknumber][0]=fixeducart;
				bc->fv[blocknumber][1]=fixedvcart;
			}
		}
	}
	if((interfacecount%2)!=0){
		printf("Error: Individual interface declarations not even!\n");
	}
	mallocInterfaces(interfaces,interfacecount/2);
	fclose(file);
}
void readGrids(Grid* grid,BoundaryConditions* bc){
	grid->xdim=mbinfo(grid->totalblocks);
	grid->xstart=(int *)malloc((grid->totalblocks)*sizeof(int));
	grid->cstart=(int *)malloc((grid->totalblocks)*sizeof(int));
	grid->estart=(int *)malloc((grid->totalblocks)*sizeof(int));
	grid->x=(double **)malloc((grid->totalblocks)*sizeof(double *));
	grid->y=(double **)malloc((grid->totalblocks)*sizeof(double *));
	int b,success,dim[3],maxchars=200,currentxstart;
	char line[maxchars];
	double xcoord,ycoord;
	currentxstart=0;
	for(b=0;b<grid->totalblocks;b++){
		//readSingleBlockGrid(grid->blockNames[b],grid->xdim[b],grid->x[b],grid->y[b]);
		FILE * file = fopen (grid->blockNames[b], "rt");
		//Get grid dimensions.
		while(fgets(line, maxchars, file) != NULL)
		{
			success = sscanf (line, "%d %d %d", &dim[0], &dim[1], &dim[2]);
			if(success>0){
				//printf("%s dimensions: %d %d %d.\n",grid->blockNames[b],dim[0],dim[1],dim[2]);
				break;
			}
		}
		grid->xdim[b][0]=dim[0];
		grid->xdim[b][1]=dim[1];
		grid->xdim[b][2]=dim[2];
		grid->xdim[b][3]=grid->xdim[b][0]*grid->xdim[b][1]*grid->xdim[b][2];
		grid->xstart[b]=currentxstart;
		//printf("xstart: %d\n",grid->xstarts[b]);
		currentxstart+=grid->xdim[b][3];
		//fill boundary traverse info.
		int s;
		s=0;
		bc->side[b][s][0]=0;
		bc->side[b][s][1]=1;
		bc->side[b][s][2]=grid->xdim[b][0]-1;
		bc->side[b][s][3]=grid->xdim[b][0];
		s=1;
		bc->side[b][s][0]=grid->xdim[b][0]-1;
		bc->side[b][s][1]=grid->xdim[b][0];
		bc->side[b][s][2]=grid->xdim[b][3]-1;
		bc->side[b][s][3]=-1;
		s=2;
		bc->side[b][s][0]=grid->xdim[b][3]-grid->xdim[b][0];
		bc->side[b][s][1]=1;
		bc->side[b][s][2]=grid->xdim[b][3]-1;
		bc->side[b][s][3]=-grid->xdim[b][0];
		s=3;
		bc->side[b][s][0]=0;
		bc->side[b][s][1]=grid->xdim[b][0];
		bc->side[b][s][2]=grid->xdim[b][3]-grid->xdim[b][0];
		bc->side[b][s][3]=1;

		//Read grid.
		grid->x[b]=malloc(sizeof(double)*grid->xdim[b][3]);
		grid->y[b]=malloc(sizeof(double)*grid->xdim[b][3]);
		int xi=0;
		while(fgets(line, maxchars, file) != NULL) {
			success = sscanf (line, "%lf %lf", &xcoord, &ycoord);
			if(success>0){
				if(xi>=grid->xdim[b][3]){ printf("Warning: More grid points were available than specified.\n"); }
				grid->x[b][xi]=xcoord;
				grid->y[b][xi]=ycoord;
				xi++;
			}
		}
		if(xi<grid->xdim[b][3]){ printf("Warning: more grid points were expected.\n"); }
		fclose(file);
	}
}
void getEntityDim(Grid* grid){
	grid->vdim=mbinfo(grid->totalblocks);
	grid->hdim=mbinfo(grid->totalblocks);
	grid->cdim=mbinfo(grid->totalblocks);
	int b;
	for(b=0;b<grid->totalblocks;b++){
		grid->cdim[b][0]=grid->xdim[b][0]-1;
		grid->cdim[b][1]=grid->xdim[b][1]-1;
		grid->cdim[b][2]=1;
		grid->cdim[b][3]=grid->cdim[b][0]*grid->cdim[b][1]*grid->cdim[b][2];
		grid->hdim[b][0]=grid->cdim[b][0];
		grid->hdim[b][1]=grid->xdim[b][1];
		grid->hdim[b][2]=1;
		grid->hdim[b][3]=grid->hdim[b][0]*grid->hdim[b][1]*grid->hdim[b][2];
		grid->vdim[b][0]=grid->xdim[b][0];
		grid->vdim[b][1]=grid->cdim[b][1];
		grid->vdim[b][2]=1;
		grid->vdim[b][3]=grid->vdim[b][0]*grid->vdim[b][1]*grid->vdim[b][2];
	}
	//calculate total edges and nodes
	grid->totaledges=0;
	grid->totalnodes=0;
	grid->totalcells=0;
	int currentestart=0,currentcstart=0;
	for(b=0;b<grid->totalblocks;b++){
		grid->totaledges+=grid->vdim[b][3];
		grid->totaledges+=grid->hdim[b][3];
		grid->totalnodes+=grid->xdim[b][3];
		grid->totalcells+=grid->cdim[b][3];
		grid->cstart[b]=currentcstart;
		grid->estart[b]=currentestart;
		currentcstart+=grid->cdim[b][3];
		currentestart+=grid->vdim[b][3]+grid->hdim[b][3];
	}
}
void fillGridBoundaryDataSub(BlockBoundaries*bbo,int d0,int d3){
	bbo->b[0].st=0;
	bbo->b[0].en=d0-1;
	bbo->b[1].st=d0-1;
	bbo->b[1].en=d3-1;
	bbo->b[2].st=d3-d0;
	bbo->b[2].en=d3-1;
	bbo->b[3].st=0;
	bbo->b[3].en=d3-d0;
	bbo->b[0].iv=bbo->b[2].iv=1;
	bbo->b[1].iv=bbo->b[3].iv=d0;
	bbo->b[0].in=d0;
	bbo->b[1].in=-1;
	bbo->b[2].in=-d0;
	bbo->b[3].in=1;
	bbo->b[0].size=bbo->b[2].size=d0;
	bbo->b[1].size=bbo->b[3].size=d3/d0;
}
void fillGridBoundaryData(Grid*g){
	int b;
	int xd0,xd3,vd0,vd3,hd0,hd3,cd0,cd3;
	g->sidexx=malloc(sizeof(BlockBoundaries)*g->totalblocks);
	g->sidevv=malloc(sizeof(BlockBoundaries)*g->totalblocks);
	g->sidehh=malloc(sizeof(BlockBoundaries)*g->totalblocks);
	g->sidecc=malloc(sizeof(BlockBoundaries)*g->totalblocks);
	for(b=0;b<g->totalblocks;b++){
		xd0=g->xdim[b][0];
		xd3=g->xdim[b][3];
		vd0=g->vdim[b][0];
		vd3=g->vdim[b][3];
		hd0=g->hdim[b][0];
		hd3=g->hdim[b][3];
		cd0=g->cdim[b][0];
		cd3=g->cdim[b][3];
		fillGridBoundaryDataSub(&g->sidexx[b],xd0,xd3);
		fillGridBoundaryDataSub(&g->sidevv[b],vd0,vd3);
		fillGridBoundaryDataSub(&g->sidehh[b],hd0,hd3);
		fillGridBoundaryDataSub(&g->sidecc[b],cd0,cd3);
	}
}
void writeMultiBlockGrid(Grid* grid,char * name){
	//Writes out formatted (ASCII) grid in whole, multi-block PLOT3D format.
	FILE * file = fopen(name, "w");
	//Write total block count.
	fprintf(file,"%d\n",grid->totalblocks);//Number of blocks.
	//Write grid dimensions.
	int b;
	for(b=0;b<grid->totalblocks;b++){
		fprintf(file,"%d\t%d\t%d\n",grid->xdim[b][0],grid->xdim[b][1],grid->xdim[b][2]);
	}
	//Write grids.
	int xi;
	int i,j,k;
	for(b=0;b<grid->totalblocks;b++){
		for (k=0;k<grid->xdim[b][2];k++) {
			for (j=0;j<grid->xdim[b][1];j++){
				for (i=0;i<grid->xdim[b][0];i++){
					xi=i+j*grid->xdim[b][0];
					fprintf(file,"\t%f",grid->x[b][xi]);
				}
				fprintf(file,"\n");
			}
			fprintf(file,"\n");
		}
		for (k=0;k<grid->xdim[b][2];k++) {
			for (j=0;j<grid->xdim[b][1];j++){
				for (i=0;i<grid->xdim[b][0];i++){
					xi=i+j*grid->xdim[b][0];
					fprintf(file,"\t%f",grid->y[b][xi]);
				}
				fprintf(file,"\n");
			}
			fprintf(file,"\n");
		}
		for (k=0;k<grid->xdim[b][2];k++) {
			for (j=0;j<grid->xdim[b][1];j++){
				for (i=0;i<grid->xdim[b][0];i++){
					fprintf(file,"\t%f",0.0);
				}
				fprintf(file,"\n");
			}
			fprintf(file,"\n");
		}
		fprintf(file,"\n");
	}
	fclose(file);
}
