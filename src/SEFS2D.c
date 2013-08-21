#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef struct {//Boundary
	int st,en,iv,in;//start,end,traverse interval,inner interval; indexed by side numbers.
	int size;
} Boundary;
typedef struct {//BlockBoundaries
	Boundary b[4];
} BlockBoundaries;
typedef struct {//BlockGridDimension
	int i,j,k,n;//(i,j,k)-dimensions,total number
} BlockGridDimension;
typedef struct {//Dimension
	int tb;//total blocks
	int tn;//total nodes
	int te;//total edges
	int tc;//total cells
	BlockGridDimension*x;
	BlockGridDimension*v;
	BlockGridDimension*h;
	BlockGridDimension*c;
	BlockBoundaries*sidexx;
	BlockBoundaries*sidevv;
	BlockBoundaries*sidehh;
	BlockBoundaries*sidecc;
	int *x0;//node starting indices for each block.
	int *e0;//edge starting indices for each block.
	int *c0;//cell starting indices for each block.
} Dimension;
typedef struct {//InvariantVelocity
	double u,v;
} InvariantVelocity;
typedef struct {//CartesianVelocity
	double u,v;
} CartesianVelocity;
typedef struct {//Property
	double nu,rho,nuinv;
} Property;
typedef struct {//State
	double **u,**v,**sanu;//the invariant velocity components u and v are actually pointers to appropriate locations in one contiguous vector that is used in the linear system solve.
	double *uvec;//to be used in the linear system. contains room for all velocity components and interface stitching relations. **u and **v point to locations within this vector.
	//Time derivatives
	double **ut,**vt,**sanut;
} State;
typedef struct {//Initial
	double initcart[2];//Initial Cartesian Velocity
	double nut0;//Eddy Viscosity
	double sanu0;//Spalart Allmaras
} Initial;
typedef struct {//Node
	double x,y;
} Node;
typedef struct {//Covariant
	double dxi[2],deta[2];
} Covariant;
typedef struct {//Contravariant
	double xi[2],eta[2];
} Contravariant;
typedef struct {//Switches
	int turbmod,wallfun;
} Switches;
typedef struct {//Numerics
	double dt,end,tol,maxiter,ptol,pmaxiter;
} Numerics;
typedef struct {//Grid
	char **blockNames;
	double **x;
	double **y;
	//Dimensions
	int totalblocks;
	int totalnodes;
	int totaledges;
	int totalcells;
	int **xdim;
	int **vdim;
	int **hdim;
	int **cdim;
	int *estart;
	int *xstart;
	int *cstart;
	BlockBoundaries*sidexx;
	BlockBoundaries*sidevv;
	BlockBoundaries*sidehh;
	BlockBoundaries*sidecc;
	//Transformation
	double **a11,**a12,**a21,**a22,
			**a11cc,**a12cc,**a21cc,**a22cc,
			**a11xx,**a12xx,**a21xx,**a22xx,
			**a11vv,**a12vv,**a21hh,**a22hh;
	double **Jcc,**Jxx,**Jvv,**Jhh,
			**vareas,**hareas;
	double **covc11,**covc12,**covc21,**covc22,
			**contrac11,**contrac12,**contrac21,**contrac22,
			**xixxx,**xiyxx,**etaxxx,**etayxx,
			**xixcc,**xiycc,**etaxcc,**etaycc;
	double **B11,**B22,**B12;
} Grid;
typedef struct {//CommonCorner
	int block[4];//block numbers involved; indexed by corner id of the respective block. -1 if no block present at the corner.
} CommonCorner;
typedef struct {//InterfaceCorners
	//Where three or more blocks interface at a corner, and is not a wall corner.
	int n;
	CommonCorner*corners;//of size n.
} InterfaceCorners;
typedef struct {//Interfaces
	FILE*out;
	int totalinterfaces;
	int totalinterfacepoints;
	double matchtolerance;//For matching interface points.
	int *block;
	int *side;
	int *adjacentBlock;
	int *xsize;
	int *rowstart;
	int **blockx;
	int **adjacentBlockx;
	int **e;//edges at interface for block.
	int **ae;//edges at interface for adjacent block;
	int **conn;//conn[b][0,1,2,3]: ab, blocks at sides
	int totalicorners;//number of corners that require special interpolation treatment (corners that are at intersection of multiple blocks).
	int **icorners;//[corner id][block (indexed by block corner id)] ::: information for interpolating cell-centered values at vertices that are at the intersection of multiple interfaces.
	//InterfaceCorners wcs;
	//InterfaceCorners ics;
} Interfaces;
typedef struct {//LinearSystem
	int crows;
	int ccols;
	double *ctcmain,*s,*ctrhs;
	double *cs,*rhs;
	double *bp,*r,*d,*q;
	double residual;
} LinearSystem;
typedef struct {//Momentum
	FILE*out;
	double**ucc,**vcc,**uxx,**vxx,**vvv,**uhh;
	double**ucartxi,**ucarteta,**vcartxi,**vcarteta;//taken at vertices.
	double**uxcc,**vxcc,**uycc,**vycc,**uxxx,**vxxx,**uyxx,**vyxx;
	double**ucartcc,**ucartxx,**ucartvv,**ucarthh,
		**vcartcc,**vcartxx,**vcartvv,**vcarthh;
	double**term14_1,**term14_2,**term15_1,**term15_2,
		**term24_1,**term24_2,**term25_1,**term25_2;
	double**I11_1,**I11_2,**I12_1,**I12_2,**I14_1,**I14_2,**I15_1,**I15_2,
		**I21_1,**I21_2,**I22_1,**I22_2,**I24_1,**I24_2,**I25_1,**I25_2;
	//Viscosity
	double**nut,**nutxx;
} Momentum;
typedef struct {//BoundaryID
	int b,s;//block,side
} BoundaryID;
typedef struct {//BC
	int n;//total number of boundaries for this BC.
	BoundaryID*id;
} BC;
typedef struct {//BoundaryConditions
	FILE*out;
	char **id;//bc id based on block and side. ( i.e. id[b][s] )
	double **fv;//fixed values for each block. Currently, one fixed vector value per block. (i.e. bc->f[b][0,1] for x- and y- components)
	int ***side;//side[block][side][start,interval,end]
	BC w;//wall
	BC f;//fixed
	BC i;//interface
	BC z;//zerogradient
	BC s;//sliding
	//The following contain corners at which 3 or 4 blocks are joined, and at which the corresponding boundary condition is of the highest precedence.
	InterfaceCorners wc;
	InterfaceCorners fc;
	InterfaceCorners ic;
} BoundaryConditions;
typedef struct {//CSR
	//int nr,*p,*r;//number of rows, pointers, row indices for the pointers (length r)
	int nn;//number of non-zeros
	int *ci,*ri;//column indices, row indices (length n)
	double *v;//values (length n)
	//combi 1 uses nn,v,ci,ri.
} CSR;
typedef struct {//SymLap
	double *A;
	double *x;
	double *b;
	int dim;
	double **d0,**d1,**d2,**d3,**d4;
	CSR csr;//for (adjacent) interface values.
	CSR csr2;//for corner interface values.
	double **cc;//cc block-divided field.
	double **xx;//field from *x (cell-centered) to vertices.
} SymLap;
typedef struct {//SpalartAllmaras
	char **bcid;//[b][s] = enum{'w','f','i'}
	double *fixed;//by block
	double fixedSanu;
	double cb1,cb2,
			sigma,
			kappa,
			cw1,cw2,cw3,
			cv1,cv2,cv3,cv1_3,
			ct1,ct2,ct3,ct4,
			rlim,
			cn1;
	double **sanuxi,**sanueta;//xi and eta derivatives.
	double **sanux,**sanuy;//x and y derivatives.
	double **dissTermx,**dissTermy;//sub-parts of dissipation term.
	double **dissTermxxi,**dissTermyxi,**dissTermxeta,**dissTermyeta;
	double **dissTermxx,**dissTermyy;//sub-parts of dissipation term.
	double **advection,**production,**wallDestruction,**dissipation;//all equation terms (except for time derivative)
	double **chi,**fv1,**fv2;
	double **ft2,**S2M,**SM,**S2,**S;
	double **r,**g,**fw;
	double **fn;
} SpalartAllmaras;
typedef struct {//WallDistances
	char **bcid;
	double *fixed;//currently just set all to zero.
	double **phi,**phixi,**phieta;
	double **cc;//wall distance field
	double **xx;//wall distance field
} WallDistances;
typedef struct {//RungeKuttaStage
	double **u,**v,**sanu;
} RungeKuttaStage;
typedef struct {//RungeKutta
	int stages;
	State*f;
} RungeKutta;
/*
typedef struct {//Pressure
	double **p;
} Pressure;
 */

#include "memory.h"

#include "other/datastructures.h"
#include "other/auxmath.h"
#include "other/state.h"
#include "other/inputs.h"

#include "linearsystem/vectorops.h"
#include "linearsystem/matrixops.h"
#include "linearsystem/efs.h"
#include "linearsystem/cmbops.h"
#include "linearsystem/cg.h"
#include "linearsystem/csr.h"

#include "invariantInterpolation/xx.h"
#include "invariantInterpolation/uhh.h"
#include "invariantInterpolation/vvv.h"
#include "invariantInterpolation/main.h"

#include "invariantToCartesian/main.h"

#include "velocityDerivatives/diff.h"
#include "velocityDerivatives/main.h"

#include "finiteVolume/convective.h"
#include "finiteVolume/viscous.h"

#include "finiteDifference/turb.h"

#include "computationalDomain/dimension.h"
#include "computationalDomain/grid.h"
#include "computationalDomain/interface.h"
#include "computationalDomain/transformation.h"
#include "computationalDomain/solution.h"
#include "computationalDomain/bcs.h"

#include "turbulenceModels/sa.h"

#include "symlap/general.h"
#include "symlap/naive.h"
#include "symlap/matfree.h"
#include "symlap/cg.h"
#include "symlap/walldistances.h"

#include "steps.h"

#include "timeAdvancement/general.h"
#include "timeAdvancement/RungeKutta.h"

#include "debug/checks.h"
#include "debug/matrix.h"

int main(void){
	State st;
	Property p;
	Numerics n;
	Switches sw;
	Grid g;
	Initial in;
	Interfaces is;
	LinearSystem ls;
	Momentum mm;
	BoundaryConditions bc;
	SymLap sl;
	SpalartAllmaras sa;
	WallDistances wd;
	RungeKutta rk;

	preprocessing(&st,&g,&bc,&is,&ls,&in,&p,&n,&sw,&sa,&rk);
	initialization(&g,&mm,&ls,&bc,&n,&st,&in,&is,&sl,&wd);
	int mode=1;

	int i;
	for(i=1;i<=n.end/n.dt;i++){
		printf("Iteration %d: ",i);
		enforceInvariantBoundary(&st,&g,&bc);

		if(mode){
			RK3(&st,&rk,&g,&bc,&is,&wd,&mm,&sa,&p,&sw,&n);
		} else {
			fillTimeDerivative(&st,&g,&bc,&is,&wd,&mm,&sa,&p,&sw,&n);
		}
		fillRhs(&g,&is,&st,&n,&ls);

		ctmult(ls.rhs,ls.ctrhs,&g,&is,&ls);
		efsCG(&ls,&n,&g,&is,&st,&bc);
		if(isnan(ls.residual)){ printf("Simulation diverged.\n"); break; }
		cmult(ls.s,st.uvec,&g,&is,&ls);
		if(sw.turbmod>0){ advanceSanu(&st,&n,&g); }
		//update(&st,&g,&ls,&bc);
	}

	writeMultiBlockGrid(&g,"out/grid.xyz");
	writeMultiBlockStateSolution("out/solution.q",&g,&st,&mm,&is,&bc,&ls);
	writeMultiBlockCustomSolution("out/nutxx.q",&g,mm.nutxx);

	writeMultiBlockCustomSolution("out/uxx.q",&g,mm.uxx);
	//writeMultiBlockCustomSolution("out/uxx.q",&g,mm.uxx);
	//writeMultiBlockCustomSolution("out/vxx.q",&g,mm.vxx);




	/*
	blockcc(sl.x,sl.cc,&g);cctoxxSL(sl.cc,wd.xx,&g,&bc,&is);
	writeMultiBlockCustomSolution("out/phi.q",&g,wd.xx);
	blockcc(sl.b,sl.cc,&g);cctoxxSL(sl.cc,wd.xx,&g,&bc,&is);
	writeMultiBlockCustomSolution("out/slrhs.q",&g,wd.xx);
	printf("Total interface points: %d\n",is.totalinterfacepoints);
	int ic=3;
	printf("icorner %d block: %d %d %d %d\n",ic,is.icorners[ic][0],is.icorners[ic][1],is.icorners[ic][2],is.icorners[ic][3]);
	FILE *f = fopen("file.txt", "w");
	if (f == NULL)
	{
	    printf("Error opening file!\n");
	    exit(1);
	}
	SymLap sl2;
	slconstructNaive(&sl2,&g,&bc,&is);
	printcsrfile(sl.csr.nn,sl.csr.v,sl.csr.ri,sl.csr.ci);//sl.csr);
	printmfile(sl2.A,sl2.dim,sl2.dim,"SL");
	printf("Symmetry: %d\n",checkSymmetry(sl2.A,sl2.dim));
	compareCsrSlA(&sl,&sl2);

	writeMultiBlockCustomSolution("out/xixxx.q",&g,g.xixxx);


	double**testcc=sa.advection;
	double**testxx=sl.xx;
	printcc(testcc,&g);
	cctoxxAllZG(testcc,testxx,&g,&is,bc.side);
	writeMultiBlockCustomSolution("out/testDebug.q",&g,testxx);



	printf("Solve Loop...\n");
	int i;
	double residual;
	for(i=1;i<=n.end/n.dt;i++){
		printf("Iteration %d: ",i);
	}
	double *Cmb=constructCmb(&g,&is);
	printmi(Cmb,g.totaledges+is.totalinterfacepoints,g.totalnodes,"Cmb");
	double *CmbT=makeAT(Cmb,ls.crows,ls.ccols);
	printmi(CmbT,ls.ccols,ls.crows,"CmbT");
	double *CTC=makeATA(Cmb,g.totaledges+is.totalinterfacepoints,g.totalnodes);
	printmi(CTC,g.totalnodes,g.totalnodes,"CTC");

	printf("Symmetry: %d\n",checkSymmetry(CTC,g.totalnodes));
	printf("Main diag match: %d\n",comparemaindiag(CTC,&ls));
	printf("ctcmult match: %d\n",comparectcmult(CTC,&g,&interfaces,&ls));
	printf("ctmult match: %d\n",comparectmult(CmbT,&g,&interfaces,&ls));
	printf("cmult match: %d\n",comparecmult(Cmb,&g,&interfaces,&ls));

	printinterfacepoints(&interfaces);
	printinterfaceedges(&interfaces);
	printJxx(&g);
	printhareas(&g);
	printvareas(&g);
	printu(&st,&g);
	printv(&st,&g);
	printcontrac(&g);
	printcartdiffvv(&mm,&g);
	printuxxx(&mm,&g);
	printucc(&mm,&g);
	printvcc(&mm,&g);
	printucartxx(&mm,&g);
	printvcartxx(&mm,&g);
	printucartcc(&mm,&g);
	printvcartcc(&mm,&g);
	printut(&mm,&g);
	printuxx(&mm,&g);
	printvxx(&mm,&g);
	printuhh(&mm,&g);
	printvvv(&mm,&g);
	printf("%g\n",n.dt);
	printcarthh(&mm,&g);
	printcartvv(&mm,&g);
	printcartdiffxx(&mm,&g);
	printf("%g\n",st.nu);
	printI2(&mm,&g);
	printf("info %d\n",is.side[0]);
	printrhs(&ls,&g);
	printctrhs(&ls,&g);
	printf("%d %d",g.cstart[0],g.cstart[1]);
	printf("\n\n\n%d\n\n\n",sl.csr.nr);
	printm(sl.A,sl.dim,sl.dim,"SL");
	printf("Symmetry: %d\n",checkSymmetry(sl.A,sl.dim));
	printf("slmultCSR match: %d\n",compareslmult(&g,&is,&sl));
	printsl(&sl,&g);
	checkcsr(&sl);
	slconstruct(&sl,&g,&bc,&is);
	//slconstructNaive(&sl,&g,&bc,&is);
	int *** ifromb=malloc(sizeof(int **)*2);
	ifromb[0]=malloc(sizeof(int *)*4);
	ifromb[0][2]=malloc(sizeof(int)*2);
	ifromb[0][2][0]=9;
	ifromb[0][2][1]=13;
	printf("It works: %d\n",ifromb[0][2][0]);
	 */

	freeAll(&st,&n,&sw,&g,&in,&is,&ls,&mm,&bc,&sl,&sa,&wd,&rk);

	printf("Finished!");
	return EXIT_SUCCESS;
}
