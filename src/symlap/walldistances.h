/*
 * walldistances.h
 *
 *  Created on: Jul 28, 2013
 *      Author: lordvon
 */

#ifndef WALLDISTANCES_H_
#define WALLDISTANCES_H_

void fillWDBC(WallDistances*wd,Grid*g,BoundaryConditions*bc){
	int b,s;
	for(b=0;b<g->totalblocks;b++){
		for(s=0;s<4;s++){
			wd->bcid[b][s]=bc->id[b][s];
			if(bc->id[b][s]=='z'){
				wd->bcid[b][s]='f';
			}
		}
	}
}
void prepWd(WallDistances*wd,Grid*g,BoundaryConditions*bc){
	mallocWd(wd,g);
	fillWDBC(wd,g,bc);
}
void fillMinDist(WallDistances*wd,Grid*g){
	int b,ci,cd3;
	double wterm1,wterm2,wterm3,A;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			wterm1=pow(g->xixcc[b][ci],2)+pow(g->xiycc[b][ci],2);
			wterm2=pow(g->etaxcc[b][ci],2)+pow(g->etaycc[b][ci],2);
			wterm3=2*(	g->xixcc[b][ci]*g->etaxcc[b][ci]+
						g->xiycc[b][ci]*g->etaycc[b][ci]	);
			A=		wterm1*pow(wd->phixi[b][ci],2)+
					wterm2*pow(wd->phieta[b][ci],2)+
					wterm3*wd->phixi[b][ci]*wd->phieta[b][ci];
			wd->cc[b][ci]=sqrt(A+2*wd->phi[b][ci])-sqrt(A);
		}
	}
}
void initWd(SymLap*s,Grid*g,BoundaryConditions*bc,Interfaces*is,
		LinearSystem*ls,Numerics*n,WallDistances*wd){
	prepWd(wd,g,bc);
	slconstruct(s,g,bc,is);
	slCG(s,ls,n,g,is);
	blockcc(s->x,wd->phi,g);
	CentralDiffCcField(is,g,
			bc->id,
			wd->fixed,
			wd->phi,
			wd->phixi,
			wd->phieta
			);
	fillMinDist(wd,g);
	cctoxxSL(wd->cc,wd->xx,g,bc,is);
	writeMultiBlockCustomSolution("out/walldistances.q",g,wd->xx);
}

#endif /* WALLDISTANCES_H_ */
