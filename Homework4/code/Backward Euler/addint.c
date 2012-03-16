void addint(float ***uf, float ***uc, float ***res, int nf)
/*
  Does coarse-to-fine interpolation and adds result to uf. nf is the fine-grid dimension. The
  coarse-grid solution is input as uc[1..nc][1..nc][1..nc], where nc = nf/2+1. The fine-grid solution
  is returned in uf[1..nf][1..nf][1..nf]. res[1..nf][1..nf][1..nf] is used for temporary storage.
*/
{
  void interp(float ***uf, float ***uc, int nf);
  int i, j, k;
  interp(res, uc, nf);
  for (i=1;i<=nf;i++)
	  for (j=1;j<=nf;j++)
		  for (k=1;k<=nf;k++)
			  uf[i][j][k] += res[i][j][k];
}
