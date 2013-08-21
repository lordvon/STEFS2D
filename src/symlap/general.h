
void prepSymlap(SymLap*s,Grid*g,Interfaces*is){
	mallocSymlap(s,g,is);
	int b,cd3;
	int tb=g->totalblocks;
	int tc=g->totalcells;
	//Initialization.
	initv(s->csr.v,s->csr.nn,0);
	initvi(s->csr.ci,s->csr.nn,0);
	initvi(s->csr.ri,s->csr.nn,0);
	initv(s->b,tc,0);
	for(b=0;b<tb;b++){
		cd3=g->cdim[b][3];
		initv(s->d0[b],cd3,0);
		initv(s->d1[b],cd3,0);
		initv(s->d2[b],cd3,0);
		initv(s->d3[b],cd3,0);
		initv(s->d4[b],cd3,0);
	}
}
