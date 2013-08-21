void checkposdef(double *s,double *ctrhs,LinearSystem* ls){
	if(dot(s,ctrhs,ls->ccols)<0){
		printf("Not positive definite!\n");
	}
}

void efsCG(LinearSystem* ls,Numerics* n,Grid* g,Interfaces* is){
	//ls->s evolves. Its current value is the starting guess.
	int dim=ls->ccols;
	int maxiter=n->maxiter;
	double tol=n->tol;
	double *x=ls->s;
	double *b=ls->ctrhs;
	double *r=ls->r;
	double *d=ls->d;
	double *q=ls->q;
	double *bp=ls->bp;
	//initv(x,dim,0);
	int refresh=10;
	int i;
	ctcmult(x,bp,g,is,ls);//checkposdef(x,bp,ls);//simpleMult(x,bp);
	a_plus_b_eq_c(b,bp,r,1,-1,dim);
	copyv(r,d,dim);
	double dnew=dot(r,r,dim);
	double d0=dnew;
	double dold;
	double thresh=pow(tol,2)*d0;
	double alpha,beta;
	for(i=0;i<maxiter;i++){
		ctcmult(d,q,g,is,ls);//checkposdef(x,bp,ls);//simpleMult(d,q);
		alpha=dnew/dot(d,q,dim);
		addon(x,alpha,d,dim);
		if(i%refresh==0){
			ctcmult(x,bp,g,is,ls);//checkposdef(x,bp,ls);//simpleMult(x,bp);
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
	printf("CG converged to a residual of %g in %d iterations.\n",ls->residual,i);
}
