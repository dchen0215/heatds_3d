/*	3D HEAT DIFFUSION SOLVER
Use FTCS and Crank-Nicolson to establish the 3D heat equation, then
implement Jacobi, Gauss-Seidel, and SOR iterative methods to solve
the linear system created by CN.
There is a ouput of the middle slice of the cubic to the text file
in oder to draw the diffusion plot using Matlab.
Change to the domain size and partitioning can be made to the front
of the program.
*/

#include "nrutil.h"
#include "Windows.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* physical domain */
#define X 1
#define Y 1
#define Z 1

/* partitioning */
#define nx 100
#define ny 100
#define nz 100

float dx = 0.01;
float dy = 0.01;
float dz = 0.01;
float threshhold = 0.000001;

#define e 2.71828183
#define nsteps 100
#define dt 0.05	/* time step */
#define alpha 0.01	/* diffusivity */
#define w 1.65

#define Cx alpha*dt/(dx*dx)
#define Cy alpha*dt/(dy*dy)
#define Cz alpha*dt/(dz*dz)


int FTCS(float ***T);		/* forward-time central-space */
int CN(float ***T);			/* crank-nicolson */
void jacobi(float ***T);	/* jacobi method */
void gaussSeidel(float ***T);	/* gauss-seidel method */
void SOR(float ***T);
float error(float ***T1, float ***T2);	/* calculate the error during iteration */
void setvalue(float ***T1, float ***T2);
void gaussnoise(float ***T);	/* Gaussian+noise */
void dirichletBC(float ***T, float t);	/* dirichlet (constant) conditions */
void periodicBC(float ***T, float x, float y, float z);	/* periodic boundary conditions */
void sourceterm(float ***T);	/* source term */
void printToFile(float ***T, char *filename);

int main(int argc, char** argv){
	int i, j;
	float ***T = f3tensor(1, nx+1, 1, ny+1, 1, nz+1);
	gaussnoise(T);		/* initial condition */
	dirichletBC(T, 0);	/* boundary conditions */
	printToFile(T, "output.txt");
	printf("\n**Crank-Nicolson**\nProblem size: %d %d %d\nRunning time: %d ms\n\n", nx, ny, nz, CN(T));
	system("pause");
}

int FTCS(float ***T){
	int start, end;		/* timer */
	int n, i, j, k;
	float ***Tnew = f3tensor(1, nx+1, 1, ny+1, 1, nz+1);
	for (n=1;n<=nsteps;n++){
		start = (int)clock();	/* time starts */
		for (i=2;i<=nx;i++){
			for (j=2;j<=ny;j++){
				for (k=2;k<=nz;k++){
					Tnew[i][j][k] = T[i][j][k] + Cx*(T[i+1][j][k] + T[i-1][j][k] - 2*T[i][j][k])
						+ Cy*(T[i][j+1][k] + T[i][j-1][k] - 2*T[i][j][k])
						+ Cz*(T[i][j][k+1] + T[i][j][k-1] - 2*T[i][j][k]);
				}
			}
		}
		setvalue(T, Tnew);	/* assign the Tnew's value to T */
		sourceterm(T);		/* source term added */
		dirichletBC(T, 0);	/* boundary conditions */
		if (n%10==0)	/* write into file every 1000 time steps */
			printToFile(T, "output.txt");
		end = (int)clock();	/* time ends */
	}
	free_f3tensor(Tnew, 1, nx+1, 1, ny+1, 1, nz+1);
	return (end-start);		/* return the running time */
}

int CN(float ***T){
	int start, end;		/* timer */
	int i, j, k, n;
	float ***Told = f3tensor(1, nx+1, 1, ny+1, 1, nz+1);
	dirichletBC(Told, 0);
	for (n=1;n<=nsteps;n++){
		/* apply the change to the right hand side matrix -- the old value
		   which used to a (nx+1)*(ny+1)*(nz+1) vector */
		for (i=2;i<=nx;i++){
			for (j=2;j<=ny;j++){
				for (k=2;k<=nz;k++){
					Told[i][j][k] = T[i][j][k] + Cx/2*(T[i+1][j][k] + T[i-1][j][k] - 2*T[i][j][k])
						+ Cy/2*(T[i][j+1][k] + T[i][j-1][k] - 2*T[i][j][k])
						+ Cz/2*(T[i][j][k+1] + T[i][j][k-1] - 2*T[i][j][k]);
				}
			}
		}
		setvalue(T, Told);
		start = (int)clock();	/* time starts */
		SOR(T);
		end = (int)clock();	/* time ends */
		if (n%10==0)	/* write into file every 1000 time steps */
			printToFile(T, "output.txt");
		
	}
	free_f3tensor(Told, 1, nx+1, 1, ny+1, 1, nz+1);
	return (end-start);		/* return the running time */
}

void jacobi(float ***T){
	float C = Cx;
	int i, j, k, n;
	int MAX_ITER = 1000;
	float ***Told = f3tensor(1, nx+1, 1, ny+1, 1, nz+1);
	float ***Tnew = f3tensor(1, nx+1, 1, ny+1, 1, nz+1);
	setvalue(Told, T);
	dirichletBC(Tnew, 0);
	for (n=1;n<=MAX_ITER;n++)
	{
		for (i=2;i<=nx;i++){
			for (j=2;j<=ny;j++){
				for (k=2;k<=nz;k++){
					Tnew[i][j][k] = C/2/(3*C + 1)*(T[i+1][j][k] +
					T[i-1][j][k]+ T[i][j+1][k] + T[i][j-1][k] + T[i][j][k+1]
					+ T[i][j][k-1]) + 1/(3*C + 1)*Told[i][j][k];
				}
			}
		}
		if (error(Tnew, T)<=threshhold)
			break;
		setvalue(T, Tnew);
	}
	setvalue(T, Tnew);
	if (n==MAX_ITER)
		printf("WARNING: iteration failed to converge");
	printf("Iterations: %d\n\n", n);
	free_f3tensor(Told, 1, nx+1, 1, ny+1, 1, nz+1);
	free_f3tensor(Tnew, 1, nx+1, 1, ny+1, 1, nz+1);
}

void gaussSeidel(float ***T){
	float C = Cx;
	int i, j, k, n;
	int MAX_ITER = 1000;
	float ***Told = f3tensor(1, nx+1, 1, ny+1, 1, nz+1);
	float ***Tlast = f3tensor(1, nx+1, 1, ny+1, 1, nz+1);
	setvalue(Told, T);
	for (n=1;n<=MAX_ITER;n++)
	{
		setvalue(Tlast, T);
		for (i=2;i<=nx;i++){
			for (j=2;j<=ny;j++){
				for (k=2;k<=nz;k++){
					T[i][j][k] = C/2/(3*C + 1)*(T[i+1][j][k] +
					T[i-1][j][k]+ T[i][j+1][k] + T[i][j-1][k] + T[i][j][k+1]
					+ T[i][j][k-1]) + 1/(3*C + 1)*Told[i][j][k];
				}
			}
		}
		if (error(Tlast, T)<=threshhold)
			break;
	}
	if (n==MAX_ITER)
		printf("WARNING: iteration failed to converge");
	printf("Iterations: %d\n\n", n);
	free_f3tensor(Told, 1, nx+1, 1, ny+1, 1, nz+1);
	free_f3tensor(Tlast, 1, nx+1, 1, ny+1, 1, nz+1);
}

void SOR(float ***T){
	float C = Cx;
	int i, j, k, n;
	int MAX_ITER = 1000;
	float ***Told = f3tensor(1, nx+1, 1, ny+1, 1, nz+1);
	float ***Tlast = f3tensor(1, nx+1, 1, ny+1, 1, nz+1);
	setvalue(Told, T);
	for (n=1;n<=MAX_ITER;n++)
	{
		setvalue(Tlast, T);
		for (i=2;i<=nx;i++){
			for (j=2;j<=ny;j++){
				for (k=2;k<=nz;k++){
					T[i][j][k] = (1 - w)*T[i][j][k] + w*C/2/(3*C + 1)*(T[i+1][j][k]
					+ T[i-1][j][k]+ T[i][j+1][k] + T[i][j-1][k] + T[i][j][k+1] + 
					T[i][j][k-1]) + w/(3*C + 1)*Told[i][j][k];
				}
			}
		}
		if (error(Tlast, T)<=threshhold)
			break;
	}
	if (n==MAX_ITER)
		printf("WARNING: iteration failed to converge");
	printf("Iterations: %d\n\n", n);
	free_f3tensor(Told, 1, nx+1, 1, ny+1, 1, nz+1);
	free_f3tensor(Tlast, 1, nx+1, 1, ny+1, 1, nz+1);
}

float error(float ***T1, float ***T2){
	int i, j, k;
	float error = 0;
	for (i=1;i<=nx+1;i++){
		for (j=1;j<=ny+1;j++){
			for (k=1;k<=nz+1;k++){
				if (error<=fabs(T1[i][j][k]-T2[i][j][k])){
					error = fabs(T1[i][j][k] - T2[i][j][k]);
				}
			}
		}
	}
	return error;
}

void setvalue(float ***T1, float ***T2){
	int i, j, k;
	for (i=1;i<=nx+1;i++)
		for (j=1;j<=ny+1;j++)
			for (k=1;k<=nz+1;k++)
				T1[i][j][k] = T2[i][j][k];
}

/* function that sets the initial condition (Gaussian+noise) */
void gaussnoise(float ***T){
	int i, j, k;
	float xx, yy, zz;
	for (i=1;i<=nx+1;i++){
		for (j=1;j<=ny+1;j++){
			for (k=1;k<=nz+1;k++){
				/* T[i][j][k] = e^(-(5*dx*(i-1)-2.5)^2-
				  (5*dy*(j-1)-2.5)^2-(5*dz*(k-1)-2.5)^2) */
				xx = pow(5*dx*(i-1)-2.5, 2);
				yy = pow(5*dy*(j-1)-2.5, 2);
				zz = pow(5*dz*(k-1)-2.5, 2);
				T[i][j][k] = pow((float)e, (float)(0-xx-yy-zz));
			}
		}
	}
}

/* dirichlet (constant) conditions */
void dirichletBC(float ***T, float t){
	int i, j, k;
	for (j=1;j<=ny+1;j++){
		for (k=1;k<=nz+1;k++){
			T[1][j][k] = t;
			T[nx+1][j][k] = t;
		}
	}
	for (i=1;i<=nx+1;i++){
		for (k=1;k<=nz+1;k++){
			T[i][1][k] = t;
			T[i][ny+1][k] = t;
		}
	}
	for (i=1;i<=nx+1;i++){
		for (j=1;j<=ny+1;j++){
			T[i][j][1] = t;
			T[i][j][nz+1] = t;
		}
	}
}

/* periodic boundary conditions */
void periodicBC(float ***T, float x, float y, float z){
	int i, j, k;
	for (j=1;j<=ny+1;j++){
		for (k=1;k<=nz+1;k++){
			T[1][j][k] = x;
			T[nx+1][j][k] = x;
		}
	}
	for (i=1;i<=nx+1;i++){
		for (k=1;k<=nz+1;k++){
			T[i][1][k] = y;
			T[i][ny+1][k] = y;
		}
	}
	for (i=1;i<=nx+1;i++){
		for (j=1;j<=ny+1;j++){
			T[i][j][1] = z;
			T[i][j][nz+1] = z;
		}
	}
}

/* arbitrary time-independent "source term" */
void sourceterm(float ***T){
	int i, j, k;
	/* here uses the ratio of the distance from (0,0,0) to
	   (i,j,k) to the distance from (0,0,0) to (X,Y,Z), and
	   multiplied by dt*/
	for (i=1;i<=nx+1;i++)
		for (j=1;j<=ny+1;j++)
			for (k=1;k<=nz+1;k++)
				T[i][j][k] += sqrt((float)(k*k + j*j + i*i))/
				sqrt((float)((nx+1)*(nx+1) + (ny+1)*(ny+1) + (nz+1)*(nz+1)))*dt;
}

/* T[1,...dx+1, 1,...dy+1, dz/2+1] should be the most typical slice of
   the 3d matirx, write it into the file upon request */
void printToFile(float ***T, char *filename){
	int i, j, k;
	FILE *fp;
    if ((fp=fopen(filename, "a"))==NULL){
         printf("Can't open file.\n");
		 exit(0);
	}
	for (i=1;i<=nx+1;i++){
		for (j=1;j<=ny+1;j++){
			fprintf(fp, "%f ", T[i][j][nz/2+1]);
		}
		fprintf(fp, "%\n");
	}
    fclose(fp);
}