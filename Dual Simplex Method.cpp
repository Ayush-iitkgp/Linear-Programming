#include<iostream>
#include<stdlib.h>
using namespace std;

void ShowMat(double **A,int m,int n) {
	for(int i=0;i<m;i++) {
		for(int j=0;j<n;j++) cout<<A[i][j]<<"\t";
		cout<<"\n";
	}
}

void Tablaeu(double **A,double *Z,int m,int n) {
	cout<<"Enter Constraints: ";
	char sign,g='>';
	for(int i=0;i<m;i++) {
		for(int j=0;j<n;j++) cin>>A[i][j];
		cin>>sign;
		cin>>A[i][n];
		if(sign==g) {
			for(int j=0;j<=n;j++) A[i][j]*=-1;
		}
	}
	for(int i=0;i<n;i++) A[m][i]=-1*Z[i];
	A[m][n]=0;
}

void Operate(double** A,int m,int n,int a,int b) {
	double p;
	double **B;
	B=(double **)malloc((m+1)*sizeof(double *));
	for(int i=0;i<=m;i++) 
		B[i]=(double *)malloc((n+1)*sizeof(double));
	for(int i=0;i<=m;i++)
		for(int j=0;j<=n;j++) B[i][j]=A[i][j];

	p=B[a][b];
	A[a][b]=1/p;
	for(int i=0;i<=m;i++) 
		if(i!=a) A[i][b]=-1*B[i][b]/p;
	for(int j=0;j<=n;j++) 
		if(j!=b) A[a][j]=B[a][j]/p;
	for(int i=0;i<=m;i++) {
		for(int j=0;j<=n;j++) {
			if((i!=a)&&(j!=b)) 
				A[i][j]= ((B[a][b]*B[i][j])-(B[a][j]*B[i][b]))/p;
		}
	}
}

double Primal(double **A,double *M,int m,int n,int im) {
	int J=-1,I;
	double min=0,PI;
	for(int j=0;j<n;j++) if(A[m][j]<min) {min=A[m][j]; J=j;}
	double sum=-1,ratio;
	for(int i=0;i<m;i++) {
		if(A[i][J]>0) {
			ratio=A[i][n]/A[i][J];
			if((sum==-1)||(sum>ratio)) {sum=ratio; I=i;}
		}
	}
	if(im==1) {
		PI=A[I][n]*A[m][J]/A[I][J];
		if(PI<0) PI*=-1;	cout<<"PI="<<PI<<"\n";
		return PI;
	}
	Operate(A,m,n,I,J);
	int temp;
	temp=M[I+n];
	M[I+n]=M[J];
	M[J]=temp;
	return 0;
}

double Dual(double **A,double *M,int m,int n,int im) {
	int I=-1,J;
	double min=0,DI;
	for(int i=0;i<m;i++) if(A[i][n]<min) {min=A[i][n]; I=i;}
	double sum=-1,ratio;
	for(int j=0;j<m;j++) {
		if(A[I][j]<0) {
			ratio=A[m][j]/A[I][j];
			if((sum==-1)||(sum>ratio)) {sum=ratio; J=j;}
		}
	}
	if(im==1) {
		DI=A[I][n]*A[m][J]/A[I][J];
		if(DI<0) DI*=-1; cout<<"DI="<<DI<<"\n";
		return DI;
	}
	Operate(A,m,n,I,J);
	int temp;
	temp=M[I+n];
	M[I+n]=M[J];
	M[J]=temp;
	return 0;
}

void simplex(double** A,double *M,int m,int n) {
	int c1=0,c2=0;
	for(int i=0;i<m;i++) if(A[i][n]<0) c1=1;	
	for(int j=0;j<n;j++) if(A[m][j]<0) c2=1;	
	if((c1==0)&&(c2==0)) return;
	else if(c1==0) Primal(A,M,m,n,0);
	else if(c2==0) Dual(A,M,m,n,0);
	else {						
	 	double PI,DI;
		PI=Primal(A,M,m,n,1);
		DI=Dual(A,M,m,n,1);
		if(PI>=DI) Primal(A,M,m,n,0);
		else Dual(A,M,m,n,0);
	} 
	simplex(A,M,m,n);
}

void Result(double**A,double*M,int m,int n) {
	cout<<"Solution:\nZ = "<<A[m][n];
	cout<<"\nAt:\n";
	for(int i=0;i<m;i++) cout<<"X"<<M[i+n]<<" = "<<A[i][n]<<"\n";
	for(int i=0;i<n;i++) cout<<"X"<<M[i]<<"=0\n";
}

int main() {
	int i,n,m;
	cout<<"Enter 0 for max 1 for min: ";
	int check=0;
	cin>>check;
	cout<<"Enter no. of variables:";
	cin>>n;
	cout<<"Enter no. of constraints:";
	cin>>m;
	double *Z,*M,**A;
	Z=(double *)malloc(n*sizeof(double));

	A=(double **)malloc((m+1)*sizeof(double *));
	for(i=0;i<=m;i++) 
		A[i]=(double *)malloc((n+1)*sizeof(double));

	M=(double *)malloc((n+m)*sizeof(double));
	for(i=0;i<m+n;i++) M[i]=i;
	
	cout<<"Enter coefficients of objective function: ";
	for(i=0;i<n;i++) cin>>Z[i];
	if(check==1)
		for(i=0;i<n;i++) Z[i]=-1*Z[i];

	Tablaeu(A,Z,m,n);

	simplex(A,M,m,n);
	if(check==1) A[m][n]*=-1;
	cout<<"\n";
	ShowMat(A,m+1,n+1);
	cout<<"\n";
	Result(A,M,m,n);
	cout<<"\n\n";
	return 0;
}
