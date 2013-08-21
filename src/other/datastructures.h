/*
 * datastructures.h
 *
 *  Created on: Jul 13, 2013
 *      Author: lordvon
 */

#ifndef DATASTRUCTURES_H_
#define DATASTRUCTURES_H_

char **charmalloc(int strings,int length){
	char **chararray=malloc(sizeof(char *)*strings);
	int i;
	for(i=0;i<strings;i++){
		chararray[i]=malloc(sizeof(char)*length);
	}
	return chararray;
}
int **mbinfo(int blocks){
	int b;
	int**mbi=malloc(sizeof(*mbi)*blocks);
	for(b=0;b<blocks;b++){
		//printf("%d\n",sizeof());
		mbi[b]=malloc(4*sizeof(**mbi));
	}
	return mbi;
}
void blockcc(double*ccvec,double**ccblock,Grid*g){
	int b,c;
	for(b=0;b<g->totalblocks;b++){
		for(c=0;c<g->cdim[b][3];c++){
			ccblock[b][c]=ccvec[g->cstart[b]+c];
		}
	}
}
int permuteCompare(int s1,int s2,int v1,int v2){
	return (((s1==v1) & (s2==v2)) | ((s1==v2) & (s2==v1)));
}
int getcornerid(int s1, int s2){
	int corner;
	if(permuteCompare(s1,s2,0,3)){
		corner=0;
	} else if(permuteCompare(s1,s2,0,1)){
		corner=1;
	} else if(permuteCompare(s1,s2,1,2)){
		corner=2;
	} else if(permuteCompare(s1,s2,2,3)){
		corner=3;
	} else {
		printf("Error: corner not identified! (getcornerid)\n");
	}
	return corner;
}
int getcornerxi(Grid*g,int block,int cornerid){
	int returnval;
	int xd0=g->xdim[block][0];
	int xd3=g->xdim[block][3];
	switch(cornerid){
	case 0:	returnval=0;		break;
	case 1:	returnval=xd0-1;	break;
	case 2:	returnval=xd3-1;	break;
	case 3:	returnval=xd3-xd0;	break;
	}
	return returnval;
}
int getcornercellci(Grid*g,int block,int cornerid){
	int returnval;
	int cd0=g->cdim[block][0];
	int cd3=g->cdim[block][3];
	switch(cornerid){
	case 0: returnval=0;		break;
	case 1:	returnval=cd0-1;	break;
	case 2:	returnval=cd3-1;	break;
	case 3:	returnval=cd3-cd0;	break;
	}
	return returnval;
}

#endif /* DATASTRUCTURES_H_ */
