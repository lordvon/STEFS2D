void cg(EFS*efs,double*u,double*rhs){
	double d0,dold,dnew,alpha,dnewThreshold;
	int i=0;
	subtract(rhs,u,efs->diff,efs->totalEdges);
	csrTranspose(efs->C,efs->diff,efs->r,efs->totalNodes);
	copyVector(efs->r,efs->d,efs->totalNodes);
	dnew=d0=dotSelf(efs->r,efs->totalNodes);
	dnewThreshold=pow(efs->residualTolerance,2)*d0;
	while((dnew>=dnewThreshold) & (i<=efs->maxIterations)){
		i++;
		csr(efs->C,efs->d,efs->uu,efs->totalEdges);
		csrTranspose(efs->C,efs->uu,efs->q,efs->totalNodes);
		alpha=dnew/dot(efs->d,efs->q,efs->totalNodes);
		addAssign(u,alpha,efs->uu,efs->totalEdges);
		if(i%efs->refresh==0){
			subtract(rhs,u,efs->diff,efs->totalEdges);
			csrTranspose(efs->C,efs->diff,efs->r,efs->totalNodes);
		} else {
			addAssign(efs->r,-alpha,efs->q,efs->totalNodes);
		}
		dold=dnew;
		dnew=dotSelf(efs->r,efs->totalNodes);
		addMultiply(efs->r,dnew/dold,efs->d,efs->d,efs->totalNodes);
	}
	efs->residual=sqrt(dnew/d0);
	efs->iterations=i;
}
