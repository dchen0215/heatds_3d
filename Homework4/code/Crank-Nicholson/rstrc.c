void rstrct(float ***uc, float ***uf, int nc) {
	/* 
	Half-weighting restriction. nc is the coarse-grid dimension. The fine-grid solution is input in
	uf[1..2*nc-1][1..2*nc-1][1...2*nc-1], the coarse-grid solution is returned in uc[1..nc][1..nc][1...nc].
	*/

	int ic, iif, jc, jf, kc, kf, ncc = 2*nc - 1;
	/* Interior points.*/
	for (iif=3,ic=2;ic<nc;ic++,iif+=2) {
		for (jf=3,jc=2;jc<nc;jc++,jf+=2) {
			for (kf=3,kc=2;kf<nc;kc++,kf+=2) {
				uc[ic][jc][kc] = 0.5*uf[iif][jf][kf] + (uf[iif+1][jf][kf] + uf[iif-1][jf][kf]
				+ uf[iif][jf+1][kf] + uf[iif][jf-1][kf] + uf[iif][jf][kf+1] + uf[iif][jf][kf-1])/12.0;
			}
		}
	}
	/* Boundary points. */
	for (iif=1,ic=1;ic<=nc;ic++,iif+=2) {
		for (jf=1,jc=1;jc<=nc;jc++,jf+=2) {
			uc[ic][jc][1] = uf[iif][jf][1];
			uc[ic][jc][nc] = uf[iif][jf][ncc];
		}
	}
	for (jf=1,jc=1;jc<=nc;jc++,jf+=2) {
		for (kf=1,kc=1;kc<=nc;kc++,kf+=2) {
			uc[1][jc][kc] = uf[1][jf][kf];
			uc[nc][jc][kc] = uf[ncc][jf][kf];
		}
	}
	for (iif=1,ic=1;ic<=nc;ic++,iif+=2) {
		for (kf=1,kc=1;kc<=nc;kc++,kf+=2) {
			uc[ic][1][kc] = uf[iif][1][kf];
			uc[ic][nc][kc] = uf[iif][ncc][kf];
		}
	}
}