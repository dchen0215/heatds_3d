void slvsml(float ***u, float ***rhs){
	/* 
	Solution of the model problem on the coarsest grid, where h = 1
	2 . The right-hand side is input
	in rhs[1..3][1..3][1..3] and the solution is returned in 
	u[1..3][1..3][1..3].
	*/

	void fill0(float ***u, int n);
	float dx = 0.5;
	float dt = 0.001;
	float alpha = 1.0;
	float C = alpha*dt/(dx*dx);
	fill0(u, 3);
	u[2][2][2] = -(1 - 3*C)*rhs[2][2][2]/(1 + 3*C);
}
