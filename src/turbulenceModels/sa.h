/*
 * sa.h
 *
 *  Created on: Aug 4, 2013
 *      Author: lordvon
 */

#ifndef SA_H_
#define SA_H_

void fillconstantsSA(SpalartAllmaras*sa){
	sa->cb1=0.1355;
	sa->cb2=0.622;
	sa->sigma=2.0/3.0,
	sa->kappa=0.41,
	sa->cw1=sa->cb1/pow(sa->kappa,2)+(1+sa->cb2)/sa->sigma;
	sa->cw2=0.3;
	sa->cw3=2;
	sa->cv1=7.1;sa->cv1_3=pow(sa->cv1,3);
	sa->cv2=0.7;
	sa->cv3=0.9;
	sa->ct1=1;
	sa->ct2=2;
	sa->ct3=1.2;
	sa->ct4=0.5;
	sa->rlim=10;
	sa->cn1=16;
}
void fillSABC(SpalartAllmaras*sa,Grid*g,BoundaryConditions*bc){
	int b,s;
	for(b=0;b<g->totalblocks;b++){
		for(s=0;s<4;s++){
			sa->bcid[b][s]=bc->id[b][s];
			if(bc->id[b][s]=='z'){
				sa->bcid[b][s]='f';
			}
		}
	}
}
void prepSA(SpalartAllmaras*sa,Grid*g,BoundaryConditions*bc,Initial*i){
	fillconstantsSA(sa);
	mallocSA(sa,g);
	fillSABC(sa,g,bc);
}

void fillChi(State*s,SpalartAllmaras*sa,Property*p,Grid*g){
	//needs sanu
	int b,cd3,ci;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			sa->chi[b][ci]=s->sanu[b][ci]*p->nuinv;
		}
	}
}
void fillft2(SpalartAllmaras*sa,Grid*g){
	//needs chi
	int b,cd3,ci;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			sa->ft2[b][ci]=sa->ct3*exp(-sa->ct4*pow(sa->chi[b][ci],2));
		}
	}
}
void fillfv1(SpalartAllmaras*sa,Grid*g){
	//needs chi
	int b,cd3,ci;
	double chi3;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			chi3=pow(sa->chi[b][ci],3);
			sa->fv1[b][ci]=chi3/(chi3+sa->cv1_3);
		}
	}
}
void fillfv2(SpalartAllmaras*sa,Grid*g){
	//needs fv1,chi
	int b,cd3,ci;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			sa->fv2[b][ci]=1-(sa->chi[b][ci]/(1+sa->chi[b][ci]*sa->fv1[b][ci]));
		}
	}
}
void fillS(SpalartAllmaras*sa,Grid*g,Momentum*mm){
	int b,cd3,ci;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			sa->S[b][ci]=abs(mm->uycc[b][ci]-mm->vxcc[b][ci]);
		}
	}
}
void fillSM(State*s,SpalartAllmaras*sa,Grid*g,double **wd){
	//assumes wall distance filled.
	//needs S,fv2,sanu
	int b,cd3,ci;
	double Sbar,k2d2,S;
	double cv22=pow(sa->cv2,2);
	double cvdiff=sa->cv3-2*sa->cv2;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			S=sa->S[b][ci];
			k2d2=pow(sa->kappa,2)*pow(wd[b][ci],2);
			Sbar=s->sanu[b][ci]/k2d2*sa->fv2[b][ci];
			if(Sbar>=-sa->cv2*S){
				sa->SM[b][ci]=S+Sbar;
			} else {
				sa->SM[b][ci]=S+S*(cv22*S+sa->cv3*Sbar)/(cvdiff*S-Sbar);
			}
		}
	}
}
void fillFn(SpalartAllmaras*sa,Grid*g){
	//needs chi
	int b,cd3,ci;
	double chi3;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			chi3=pow(sa->chi[b][ci],3);
			sa->fn[b][ci]=(sa->cn1+chi3)/(sa->cn1-chi3);
		}
	}
}
void fillProduction(State*s,SpalartAllmaras*sa,Grid*g){
	//needs ft2,SM,sanu
	int b,cd3,ci;
	double coeff;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			if(s->sanu[b][ci]>=0){
				coeff=1-sa->ft2[b][ci];
			} else {
				coeff=1-sa->ct3;
			}
			sa->production[b][ci]=sa->cb1*
					coeff*
					sa->SM[b][ci]*
					s->sanu[b][ci];
		}
	}
}
void fillR(State*s,SpalartAllmaras*sa,Grid*g,double**wd){
	//needs wd,SM,sanu
	int b,cd3,ci;
	double k2d2,min;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			k2d2=pow(sa->kappa,2)*pow(wd[b][ci],2);
			min=s->sanu[b][ci]/sa->SM[b][ci]/k2d2;
			if((sa->rlim<min) | (isnan(min))){ min=sa->rlim; }
			sa->r[b][ci]=min;
		}
	}
}
void fillG(SpalartAllmaras*sa,Grid*g){
	//needs r
	int b,cd3,ci;
	double r;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			r=sa->r[b][ci];
			sa->g[b][ci]=r+sa->cw2*(pow(r,6)-r);
		}
	}
}
void fillFw(SpalartAllmaras*sa,Grid*g){
	//needs g
	int b,cd3,ci;
	double cw3_6=pow(sa->cw3,6);
	double q;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			q=(1+cw3_6)/(pow(sa->g[b][ci],6)+cw3_6);
			sa->fw[b][ci]=sa->g[b][ci]*pow(q,1.0/6.0);
		}
	}
}
void fillWallDestruction(State*s,SpalartAllmaras*sa,Grid*g,double**wd){
	//needs wd,ft2,fw,sanu
	int b,cd3,ci;
	double k2inv=1/pow(sa->kappa,2);
	double coeff;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			if(s->sanu[b][ci]>=0){
				coeff=sa->cw1*sa->fw[b][ci]-sa->cb1*k2inv*sa->ft2[b][ci];
			} else {
				coeff=-sa->cw1;
			}
			sa->wallDestruction[b][ci]=
					coeff*pow(s->sanu[b][ci]/wd[b][ci],2);
		}
	}
}
void fillDissTerms(State*s,SpalartAllmaras*sa,Grid*g,Property*p,Interfaces*is){
	//needs
	int b,cd3,ci;
	double nuterm;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			if(s->sanu[b][ci]>=0){
				nuterm=p->nu+s->sanu[b][ci];
			} else {
				nuterm=p->nu+s->sanu[b][ci]*sa->fn[b][ci];
			}
			sa->dissTermx[b][ci]=nuterm*sa->sanux[b][ci];
			sa->dissTermy[b][ci]=nuterm*sa->sanuy[b][ci];
		}
	}
	CentralDiffCcFieldFOB(is,g,
			sa->dissTermx,sa->dissTermxxi,sa->dissTermxeta);
	CentralDiffCcFieldFOB(is,g,
			sa->dissTermy,sa->dissTermyxi,sa->dissTermyeta);
	//now get cartesian derivatives
	fillCartDerivativeXCc(g,sa->dissTermxx,
			sa->dissTermxxi,sa->dissTermxeta);
	fillCartDerivativeYCc(g,sa->dissTermyy,
			sa->dissTermyxi,sa->dissTermyeta);
}
void fillDissipation(SpalartAllmaras*sa,Grid*g){
	int b,cd3,ci;
	double sigmaInv=1/sa->sigma;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			sa->dissipation[b][ci]=sigmaInv*(
					sa->dissTermxx[b][ci]+
					sa->dissTermyy[b][ci]+
					sa->cb2*(
							pow(sa->sanux[b][ci],2)+
							pow(sa->sanuy[b][ci],2)
							)
					);
		}
	}
}
void fillAdvectionCD(SpalartAllmaras*sa,Grid*g,Momentum*mm){
	//CD: central difference.
	int b,cd3,ci;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			sa->advection[b][ci]=
					mm->ucartcc[b][ci]*sa->sanux[b][ci]+
					mm->vcartcc[b][ci]*sa->sanuy[b][ci];
		}
	}
}
void computeSanut(State*s,SpalartAllmaras*sa,Property*p,Grid*g,Interfaces*is,Momentum*mm,
		double**wd){
	//assumes cart derivatives are updated.
	//spatial derivatives of sanu.
	CentralDiffCcField(is,g,sa->bcid,sa->fixed,s->sanu,sa->sanuxi,sa->sanueta);
	fillCartDerivativeCc(g,sa->sanux,sa->sanuy,sa->sanuxi,sa->sanueta);
	//fill sa terms.
	fillChi(s,sa,p,g);
	fillft2(sa,g);
	fillfv1(sa,g);
	fillfv2(sa,g);
	fillS(sa,g,mm);
	fillSM(s,sa,g,wd);
	fillFn(sa,g);
	fillProduction(s,sa,g);
	fillR(s,sa,g,wd);
	fillG(sa,g);
	fillFw(sa,g);
	fillWallDestruction(s,sa,g,wd);
	fillDissTerms(s,sa,g,p,is);
	fillDissipation(sa,g);
	fillAdvectionCD(sa,g,mm);
	//compute time derivative.
	int b,cd3,ci;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			s->sanut[b][ci]=sa->production[b][ci]-sa->wallDestruction[b][ci]+
					sa->dissipation[b][ci]-sa->advection[b][ci];
		}
	}
}
void advanceSanu(State*s,Numerics*n,Grid*g){
	int b,cd3,ci;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			s->sanu[b][ci]+=s->sanut[b][ci]*n->dt;
		}
	}
}
void updateNut(State*s,SpalartAllmaras*sa,Grid*g,Property*p,Momentum*mm,
		BoundaryConditions*bc,Interfaces*is){
	int b,cd3,ci;
	for(b=0;b<g->totalblocks;b++){
		cd3=g->cdim[b][3];
		for(ci=0;ci<cd3;ci++){
			if(s->sanu[b][ci]>=0){
				mm->nut[b][ci]=s->sanu[b][ci]*sa->fv1[b][ci];
			} else {
				mm->nut[b][ci]=0;
			}
		}
	}
	cctoxxSA(mm->nut,mm->nutxx,sa,p,g,bc,is);
}

#endif /* SA_H_ */
