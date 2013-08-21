void inline diagmult(double*d,double*x,double*b,int dlen,int offset,int r){
	int i,row,col;
	for(i=0;i<dlen;i++){
		row=r+i;col=row+offset;
		b[row]+=x[col]*d[i];
		//if(row==0){ printf("%g x %g = %g\n",x[col],d[i],x[col]*d[i]); }
	}
}
void inline symmdiagmult(double*d,double*x,double*b,int dlen,int offset,int r){
	int i,row,col;
	for(i=0;i<dlen;i++){
		row=r+i;col=row+offset;
		b[row]+=x[col]*d[i];
		//if(row==0){ printf("%g x %g = %g\n",x[col],d[i],x[col]*d[i]); }
		b[col]+=x[row]*d[i];
	}
}
void slmultCSR(double*x,double*b,SymLap*s,Grid*g,Interfaces*is){
	int block,cd0,cd3,cs;
	initv(b,s->dim,0);
	for(block=0;block<g->totalblocks;block++){
		cs=g->cstart[block];
		cd0=g->cdim[block][0];
		cd3=g->cdim[block][3];
		diagmult	(s->d0[block],x,b,cd3,		0,cs);
		symmdiagmult(s->d1[block],x,b,cd3-1,	1,cs);
		symmdiagmult(s->d2[block],x,b,cd3-cd0+1,cd0-1,cs);
		symmdiagmult(s->d3[block],x,b,cd3-cd0,	cd0,cs);
		symmdiagmult(s->d4[block],x,b,cd3-cd0-1,cd0+1,cs);
	}
	csrmult1symmetric(&s->csr,x,b);
	csrmult1symmetric(&s->csr2,x,b);
}
void slCG(SymLap* sl,LinearSystem*ls,Numerics*n,Grid*g,Interfaces*is){
	//ls->s evolves. Its current value is the starting guess.
	int dim=sl->dim;
	int maxiter=n->maxiter;
	double tol=n->tol;
	double *x=sl->x;
	double *b=sl->b;
	double *r=ls->r;
	double *d=ls->d;
	double *q=ls->q;
	double *bp=ls->bp;
	int refresh=10;
	int i;
	slmultCSR(x,bp,sl,g,is);//checkposdef(x,bp,ls);//simpleMult(x,bp);
	a_plus_b_eq_c(b,bp,r,1,-1,dim);
	copyv(r,d,dim);
	double dnew=dot(r,r,dim);
	double d0=dnew;
	double dold;
	double thresh=pow(tol,2)*d0;
	double alpha,beta;
	for(i=0;i<maxiter;i++){
		slmultCSR(d,q,sl,g,is);//ctcmult(d,q,g,is,ls);//checkposdef(x,bp,ls);//simpleMult(d,q);
		alpha=dnew/dot(d,q,dim);
		addon(x,alpha,d,dim);
		if(i%refresh==0){
			slmultCSR(x,bp,sl,g,is);//ctcmult(x,bp,g,is,ls);//checkposdef(x,bp,ls);//simpleMult(x,bp);
			a_plus_b_eq_c(b,bp,r,1,-1,dim);
		} else {
			addon(r,-alpha,q,dim);
		}
		dold=dnew;
		dnew=dot(r,r,dim);
		beta=dnew/dold;
		a_plus_b_eq_c(r,d,d,1,beta,dim);
		i++;
		if(dnew<thresh){break;}
		if(isnan(dnew)){break;}
	}
	ls->residual=sqrt(dnew/d0);
	printf("Symmetric Laplacian CG converged to a residual of %g in %d iterations.\n",ls->residual,i);
}
