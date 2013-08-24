double dot(double*a,double*b,int l){
	int i;
	double tot=0;
	for(i=0;i<l;i++){
		tot+=a[i]*b[i];
	}
	return tot;
}
double dotSelf(double*a,int l){
	int i;
	double tot=0;
	for(i=0;i<l;i++){
		tot+=a[i]*a[i];
	}
	return tot;
}
void add(double*a,double*b,double*c,int length){
	//performs a+b=c (c is overwritten)
	int i;
	for(i=0;i<length;i++){
		c[i]=a[i]+b[i];
	}
}
void addAssign(double*a,double coeff,double*b,int length){
	//performs a+=coeff*b (a is overwritten)
	int i;
	for(i=0;i<length;i++){
		a[i]+=coeff*b[i];
	}
}
void addMultiply(double*a,double coeff,double*b,double*c,int length){
	//performs a+coeff*b=c (c is overwritten)
	int i;
	for(i=0;i<length;i++){
		c[i]=a[i]+coeff*b[i];
	}
}
void subtract(double*a,double*b,double*c,int length){
	//performs a-b=c (c is overwritten)
	int i;
	for(i=0;i<length;i++){
		c[i]=a[i]-b[i];
	}
}
void copyVector(double*a,double*b,int length){
	//copy a to b.
	int i;
	for(i=0;i<length;i++){
		b[i]=a[i];
	}
}
void initializeVector(double*a,int length,double value){
	int i;
	for(i=0;i<length;i++){
		a[i]=value;
	}
}
