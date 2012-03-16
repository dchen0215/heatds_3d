
#define NPRE  2
#define NPOST 2
#define NGMAX 15

void addint(float ***uf, float ***uc, float ***res, int nf);
void copy(float ***aout, float ***ain, int n);
void fill0(float ***u, int n);
void interp(float ***uf, float ***uc, int nf);
void relax(float ***u, float ***rhs, int n);
void resid(float ***res, float ***u, float ***rhs, int n);
void rstrct(float ***uc, float ***uf, int nc);
void slvsml(float ***u, float ***rhs);
void printToFile(float ***u, int n);
void initialization(float ***u, int n);
void gaussSeidel_RB(float ***u, float ***rhs, int n, float C);
void gaussSeidel(float ***u, float ***rhs, int n, float C);
void jacobi(float ***u, float ***rhs, int n, float C);
