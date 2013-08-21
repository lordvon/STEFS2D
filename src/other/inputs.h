/*
 * inputs.h
 *
 *  Created on: Jul 23, 2013
 *      Author: lordvon
 */

#ifndef INPUTS_H_
#define INPUTS_H_

void initialDict(Initial* i,char *var, double* val){
	int matchlimit=4;
	if(strncasecmp(var, "uv0", matchlimit)==0){
		i->initcart[0]=val[0];i->initcart[1]=val[1]; }
	if(strncasecmp(var, "nut0", matchlimit)==0){ i->nut0=val[0]; }
	if(strncasecmp(var, "sanu0", matchlimit)==0){ i->sanu0=val[0]; }
}
void numericsDict(Numerics* n,SpalartAllmaras*s,char *var, double* val){
	int matchlimit=4;
	if(strncasecmp(var, "dt", matchlimit)==0){ n->dt=val[0]; }
	if(strncasecmp(var, "end", matchlimit)==0){ n->end=val[0]; }
	if(strncasecmp(var, "tol", matchlimit)==0){ n->tol=val[0]; }
	if(strncasecmp(var, "maxiter", matchlimit)==0){ n->maxiter=val[0]; }
	if(strncasecmp(var, "ptol", matchlimit)==0){ n->ptol=val[0]; }
	if(strncasecmp(var, "pmaxiter", matchlimit)==0){ n->pmaxiter=val[0]; }
	if(strncasecmp(var, "fixedSanu", matchlimit)==0){ s->fixedSanu=val[0]; }
}
void propertyDict(Property*p,char *var, double* val){
	int matchlimit=4;
	if(strncasecmp(var, "nu", matchlimit)==0){ p->nu=val[0];p->nuinv=1/val[0]; }
	if(strncasecmp(var, "rho", matchlimit)==0){ p->rho=val[0]; }
}
void switchDict(Switches* s,char *var, double* val){
	int matchlimit=4;
	if(strncasecmp(var, "turbmod", matchlimit)==0){ s->turbmod=val[0]; }
	if(strncasecmp(var, "wallfun", matchlimit)==0){ s->wallfun=val[0]; }
}
void readInitial(Initial* i,char* fn){
	int maxchars=200;
	char line[maxchars];
	char var[maxchars];
	double val[2];
	FILE * file = fopen (fn, "rt");
	int success;
	//Read total number of blocks.
	while(fgets(line, maxchars, file) != NULL){
		success = sscanf(line, "%s %lf %lf",var,&val[0],&val[1]);
		if(success>0){ initialDict(i,var,val); }
	}
	fclose(file);
}
void readNumerics(Numerics* n,SpalartAllmaras*s,char* fn){
	int maxchars=200;
	char line[maxchars];
	char var[maxchars];
	double val[2];
	FILE * file = fopen (fn, "rt");
	int success;
	//Read total number of blocks.
	while(fgets(line, maxchars, file) != NULL){
		success = sscanf(line, "%s %lf %lf",var,&val[0],&val[1]);
		if(success>0){ numericsDict(n,s,var,val); }
	}
	fclose(file);
}
void readProperty(Property* p,char* fn){
	int maxchars=200;
	char line[maxchars];
	char var[maxchars];
	double val[2];
	FILE * file = fopen (fn, "rt");
	int success;
	//Read total number of blocks.
	while(fgets(line, maxchars, file) != NULL){
		success = sscanf(line, "%s %lf %lf",var,&val[0],&val[1]);
		if(success>0){ propertyDict(p,var,val); }
	}
	fclose(file);
}
void readSwitch(Switches* s,char* fn){
	int maxchars=200;
	char line[maxchars];
	char var[maxchars];
	double val[2];
	FILE * file = fopen (fn, "rt");
	int success;
	//Read total number of blocks.
	while(fgets(line, maxchars, file) != NULL){
		success = sscanf(line, "%s %lf %lf",var,&val[0],&val[1]);
		if(success>0){ switchDict(s,var,val); }
	}
	fclose(file);
}

#endif /* INPUTS_H_ */
