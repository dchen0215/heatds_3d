1. Make change to the physical domain, partitioning, time step, etc as needed. Note that 
nx & dx, ny & dy, and nz & dz should be changed at the same time.
2. In line 112, body of the Crank-Nicolson, is where the iterative method is called. By
default, SOR is used to solve the linear system. You can change that to Jacobi or Gauss-
Seidel manually.
3. In the main function, Crank-Nicolson is called and the program will print out the
runtime per algorithmic timestep of the current iterative method that is being used.
4. The assumption made to the iterative methods is that it is a cubic domain and has the
same partitioning for each dimension, so the C is fixed for x, y, and z. You will find that
at the first line of each method, Cx is assigned to C.
5. The iterative method will give the times of iteration after each time step for the sake
of comparison and checking if the iterative methods are working right.
6. Every 10 time steps, the temperature of one slice of the cubic will be writen into the
text file once.
7. error() function is to calculate the error after each iteration to control the iterative
procedure.