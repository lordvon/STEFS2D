void preprocessing(State*s,Grid*g,BoundaryConditions*bc,Interfaces*is,LinearSystem*ls,
		Initial*in,Property*p,Numerics*n,Switches*sw,SpalartAllmaras*sa,RungeKutta*rk){
	//Read inputs.
	readBlockDict(g,is,bc,"in/blockDict");
	readGrids(g,bc);
	getEntityDim(g);
	fillGridBoundaryData(g);
	readInitial(in,"in/initial");
	readProperty(p,"in/property");
	readNumerics(n,sa,"in/numerics");
	readSwitch(sw,"in/switch");
	//Interfacing.
	is->matchtolerance=1e-6;
	establishInterfaces(g,is,bc);
	fillConnectivity(g,is);
	fillInterfaceCoords(g,is);
	fillInterfaceEdges(g,is);
	fillicorners(is,g);

	//Ax=b system
	prepls(g,is,ls);
	fillmaindiag(g,is,ls);

	//Grid quantities
	transformation(g,bc,is);

	//Turbulence Models
	prepSA(sa,g,bc,in);

	//Time stepping schemes.
	rk->stages=3;
	mallocRK(g,rk);
}
void initialization(Grid*g,Momentum*mm,LinearSystem*ls,BoundaryConditions*bc,
		Numerics*n, State*s,Initial*in,Interfaces*is,SymLap*sl,
		WallDistances*wd){
	mallocMomentum(s,g,mm,is);
	initState(g,mm,ls,s,in);
	//Wall distances
	initWd(sl,g,bc,is,ls,n,wd);
}
void interpolateState(State* s,Momentum* mm,Grid* g,
		Interfaces* is,BoundaryConditions* bc){
	fillcc(s,mm,g);
	fillxx(s,mm,g,bc,is);
	filluhh(s,g,mm,bc,is);
	fillvvv(s,g,mm,bc,is);
}
void cartesianConvert(State* s,Grid* g,Momentum* mm){
	//Assumes invariant field has been interpolated to all locations.
	fillCartcc(g,mm);
	fillCartxx(g,mm);
	fillCarthh(s,g,mm);
	fillCartvv(s,g,mm);
}
void cartesianDerivatives(Grid* g,Momentum* mm,Interfaces* is){
	fillDerivativecc(g,mm);
	fillDerivativexx(g,mm,is);
}
void convection(Grid* g,Momentum* mm,BoundaryConditions* bc,Interfaces* is){
	fillI11(g,mm,bc,is);
	fillI12(g,mm);
	fillI21(g,mm);
	fillI22(g,mm,bc,is);
}
void viscosity(Property*p,Grid* g,Momentum* mm,BoundaryConditions* bc,
		Interfaces* is){
	fillterms14and25(p,g,mm);
	fillterms15and24(p,g,mm);
	fillI14(g,mm,bc,is);
	fillI15(g,mm);
	fillI24(g,mm);
	fillI25(g,mm,bc,is);
}
