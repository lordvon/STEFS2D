int countBlocks(){
	int blocknumber,totalblocks,maxchars=200,success;
	char bc0,bc1,bc2,bc3,line[maxchars],filename[maxchars];
	//Read total number of blocks.
	totalblocks=0;
	FILE * file = fopen ("in/blockDict", "rt");
	while(fgets(line, maxchars, file) != NULL){
		success = sscanf (line, "%d %s %c %c %c %c",&blocknumber,filename,
				&bc0,&bc1,&bc2,&bc3);
		//printf("line: %c %c %c %c\n",line[0],line[1],line[2],line[3]);
		if(line[0]=='#'){ break; }
		if(success>0){ totalblocks++; }
	}
	fclose(file);
	return totalblocks;
}
void fillBlockFileNames(char**blockFileNames){
	int b,maxchars=200;
	char line[maxchars],filename[maxchars];
	FILE * file = fopen ("in/blockDict", "rt");
	int success;
	//Read file names for each block.
	while(fgets(line, maxchars, file) != NULL){
		success = sscanf (line, "%d %s",&b,filename);
		if(line[0]=='#'){ break; }
		if(success>0){ strcpy(blockFileNames[b],filename); }
	}
	fclose(file);
}
void countNodes(char**blockFileNames,Dimension*d){
	int b,success,dim[4],maxchars=200,currentxstart;
	char line[maxchars];
	currentxstart=0;
	for(b=0;b<d->tb;b++){
		//readSingleBlockGrid(grid->blockNames[b],grid->xdim[b],grid->x[b],grid->y[b]);
		FILE * file = fopen (blockFileNames[b], "rt");
		//Get grid dimensions.
		while(fgets(line, maxchars, file) != NULL)
		{
			success = sscanf (line, "%d %d %d", &dim[0], &dim[1], &dim[2]);
			if(success>0){ break; }
		}
		dim[3]=dim[0]*dim[1]*dim[2];
		d->x[b].i=dim[0];
		d->x[b].j=dim[1];
		d->x[b].k=dim[2];
		d->x[b].n=dim[3];
		d->x0[b]=currentxstart;
		currentxstart+=dim[3];
		fclose(file);
	}
	d->tn=currentxstart;
}
void inferEntityCounts(Dimension*d){
	int b,te,tc;
	te=tc=0;
	for(b=0;b<d->tb;b++){
		d->c[b].i=d->x[b].i-1;
		d->c[b].j=d->x[b].j-1;
		d->v[b].i=d->x[b].i;
		d->v[b].j=d->x[b].j-1;
		d->h[b].i=d->x[b].i-1;
		d->h[b].j=d->x[b].j;

		d->c[b].k=1;
		d->v[b].k=1;
		d->h[b].k=1;
		d->c[b].n=d->c[b].i*d->c[b].j*d->c[b].k;
		d->v[b].n=d->v[b].i*d->v[b].j*d->v[b].k;
		d->h[b].n=d->h[b].i*d->h[b].j*d->h[b].k;

		d->e0[b]=te;te+=d->v[b].n+d->h[b].n;
		d->c0[b]=tc;tc+=d->c[b].n;
	}
	d->te=te;
	d->tc=tc;
}
void fillSidesSub(BlockBoundaries*bo,int d0,int d3){
	bo->b[0].st=0;
	bo->b[0].en=d0-1;
	bo->b[1].st=d0-1;
	bo->b[1].en=d3-1;
	bo->b[2].st=d3-d0;
	bo->b[2].en=d3-1;
	bo->b[3].st=0;
	bo->b[3].en=d3-d0;
	bo->b[0].iv=bo->b[2].iv=1;
	bo->b[1].iv=bo->b[3].iv=d0;
	bo->b[0].in=d0;
	bo->b[1].in=-1;
	bo->b[2].in=-d0;
	bo->b[3].in=1;
	bo->b[0].size=bo->b[2].size=d0;
	bo->b[1].size=bo->b[3].size=d3/d0;
}
void fillSides(Dimension*d){
	int b;
	for(b=0;b<d->tb;b++){
		fillSidesSub(&d->xside[b],d->x->i,d->x->n);
		fillSidesSub(&d->vside[b],d->v->i,d->v->n);
		fillSidesSub(&d->hside[b],d->h->i,d->h->n);
		fillSidesSub(&d->cside[b],d->c->i,d->c->n);
	}
}
void fillDimension(Dimension*d){
	d->tb=countBlocks();
	char**blockFileNames=charmalloc(d->tb,200);
	mallocDimension(d);
	fillBlockFileNames(blockFileNames);
	countNodes(blockFileNames,d);
	inferEntityCounts(d);
	fillSides(d);
	//free
	int b;
	for(b=0;b<d->tb;b++){
		free(blockFileNames[b]);
	}
	free(blockFileNames);
}
