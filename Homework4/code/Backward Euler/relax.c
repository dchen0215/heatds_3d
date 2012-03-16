#include "nrutil.h"
#include "mg.h"

void relax(float ***u, float ***rhs, int n){
	/*
	Red-black Gauss-Seidel relaxation for model problem. Updates the current value of the solution
	u[1..n][1..n][1..n], using the right-hand side function rhs[1..n][1..n][1..n].
	*/

	float dx, C, dt = 0.001, alpha = 1.0;
	dx = 1.0/(n - 1);
	C = alpha*dt/(dx*dx);
	jacobi(u, rhs, n, C);	
}

void gaussSeidel_RB(float ***u, float ***rhs, int n, float C){
	int i, ipass, isw, j, jsw=1, k , ksw=1;
	/* Red and black sweeps.*/
	/* jsw and isw toggle between 1 and 2 and
	determine starting row in each column
	for given pass 
	*/

	/*Gauss-Seidel with black and red ordering formula.*/
	for (ipass=1;ipass<=2;ipass++,ksw=3-ksw) {
		jsw = ksw;	
		for (k=2;k<n;k++,jsw=3-jsw)
			isw = jsw;
			for (j=2;j<n;j++,isw=3-isw)
				for (i=isw+1;i<n;i+=2) 
					u[i][j][k] = C/(6*C + 1)*(u[i+1][j][k] + u[i-1][j][k] + u[i][j+1][k]
				+ u[i][j-1][k] + u[i][j][k-1] + u[i][j][k+1]) - rhs[i][j][k]/(6*C + 1);
	}
}

void gaussSeidel(float ***u, float ***rhs, int n, float C){
	int i, j, k;
	/* Gauss-Seidel formula. */
	for (i=2;i<n;i++)
		for (j=2;j<n;j++)
			for (k=2;k<n;k++)
				u[i][j][k] = C/(6*C + 1)*(u[i+1][j][k] + u[i-1][j][k] + u[i][j+1][k]
				+ u[i][j-1][k] + u[i][j][k-1] + u[i][j][k+1]) - rhs[i][j][k]/(6*C + 1);

}

void jacobi(float ***u, float ***rhs, int n, float C){
	int i, j, k;
	float ***unew = f3tensor(1, n, 1, n, 1, n);
	/* Jacobi formula
	   it requires declaring an extra matrix to store the result. */
	fill0(unew, n);
	for (i=2;i<n;i++)
		for (j=2;j<n;j++)
			for (k=2;k<n;k++)
				unew[i][j][k] = C/(6*C + 1)*(u[i+1][j][k] + u[i-1][j][k] + u[i][j+1][k]
				+ u[i][j-1][k] + u[i][j][k-1] + u[i][j][k+1]) - rhs[i][j][k]/(6*C + 1);
	copy(u, unew, n);
}
