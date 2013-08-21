/*
 * vectorops.h
 *
 *  Created on: Jul 16, 2013
 *      Author: lordvon
 */

#ifndef VECTOROPS_H_
#define VECTOROPS_H_

void initvi(int * vec, int len, int val){
	int i;
	for (i=0;i<len;i++){
		vec[i]=val;
	}
}
void initv(double*vec,int len,double val){
	int i;
	for (i=0;i<len;i++){
		vec[i]=val;
	}
}
double dot(double * a, double * b, int l){
	int i;
	double tot=0;
	for(i=0;i<l;i++){
		tot+=a[i]*b[i];
	}
	return tot;
}
void a_plus_b_eq_c(double * a, double * b, double * c, double asign, double bsign, int length){
	int i;
	for(i=0;i<length;i++){
		c[i]=asign*a[i]+bsign*b[i];
	}
}
void addon(double * base, double coeff, double * toadd, int length){
	int i;
	for(i=0;i<length;i++){
		base[i]+=coeff*toadd[i];
	}
}
void copyv(double * a, double * b, int length){
	//copy a to b.
	int i;
	for(i=0;i<length;i++){
		b[i]=a[i];
	}
}

#endif /* VECTOROPS_H_ */
