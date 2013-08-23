void cg(EFS*efs,double*u,double*rhs){
	double d0,dold,dnew,alpha,dnewThreshold;
	int i=1;
	subtract(rhs,u,efs->diff,efs->totalEdges);
	csrTranspose(efs->C,efs->diff,efs->r);
	copyVector(efs->r,efs->d,efs->totalNodes);
	dnew=d0=dotSelf(efs->r,efs->totalNodes);
	dnewThreshold=pow(efs->residualTolerance,2)*d0;
	while((dnew>=dnewThreshold) & (i<=efs->maxIterations)){
		csr(efs->C,efs->d,efs->uu);
		csrTranspose(efs->C,efs->uu,efs->q);
		alpha=dnew/dot(efs->d,efs->q,efs->totalNodes);
		addAssign(u,alpha,efs->uu,efs->totalEdges);
		if(i%efs->refresh==0){
			subtract(rhs,u,efs->diff,efs->totalEdges);
			csrTranspose(efs->C,efs->diff,efs->r);
		} else {
			addAssign(efs->r,-alpha,efs->q,efs->totalNodes);
		}
		dold=dnew;
		dnew=dotSelf(efs->r,efs->totalNodes);
		addMultiply(efs->r,dnew/dold,efs->d,efs->d,efs->totalNodes);
	}
}
