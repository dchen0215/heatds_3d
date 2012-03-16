#include "mg.h"
#include "nrutil.h"
#include "math.h"
#include <stdio.h>
#include <time.h>
#include <windows.h>

#define e 2.71828183

int main(int argc, char **argv) {
	int i, j, k;
	int start, end;
	float ***f;
	int n = 257;
	int ncycle = 2;
	f = f3tensor(1, n, 1, n, 1, n);
	/* initialize the matrix with gaussian and periodic boundary conditions */
	initialization(f, n);

	start = (int)clock();
	mglin(f, n, ncycle);
	end = (int)clock();
	printf("\n**Gauss-Seidel(R&B) + BE**\nProblem size: %d\nRunning time: %d ms\n\n", n, end - start);
	/* write the result into text file */
	printToFile(f, n);
	system("pause");
}

void printToFile(float ***u, int n){
	int i, j;
	FILE *fp;
    if ((fp=fopen("soln.txt", "w"))==NULL){
         printf("Can't open file.\n");
	}
	for (i=1;i<=n;i++){
		for (j=1;j<=n;j++){
			fprintf(fp, "%f ", u[i][j][2]);
		}
		fprintf(fp, "%\n");
	}
    fclose(fp);
}

void initialization(float ***u, int n){
	int i, j, k;
	float xx, yy, zz;
	float dx = 1.0/(n-1);
	for (i=1;i<=n;i++){
		for (j=1;j<=n;j++){
			for (k=1;k<=n;k++){
				xx = pow(5.0*dx*(i-1)-2.5, 2);
				yy = pow(5.0*dx*(j-1)-2.5, 2);
				zz = pow(5.0*dx*(k-1)-2.5, 2);
				u[i][j][k] = pow((float)e, (float)(0-xx-yy-zz));
			}
		}
	}
	for (j=1;j<=n;j++){
		for (k=1;k<=n;k++){
			u[1][j][k] = 0.0;
			u[n][j][k] = 0.0;
		}
	}
	for (i=1;i<=n;i++){
		for (k=1;k<=n;k++){
			u[i][1][k] = 0.0;
			u[i][n][k] = 0.0;
		}
	}
	for (i=1;i<=n;i++){
		for (j=1;j<=n;j++){
			u[i][j][1] = 0.0;
			u[i][j][n] = 0.0;
		}
	}
}