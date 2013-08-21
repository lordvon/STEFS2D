/*
 * auxmath.h
 *
 *  Created on: Jul 18, 2013
 *      Author: lordvon
 */

#ifndef AUXMATH_H_
#define AUXMATH_H_

double cross(double * a, double * b){
	//for 2x1 vectors.
	return a[0]*b[1]-a[1]*b[0];
}

#endif /* AUXMATH_H_ */
