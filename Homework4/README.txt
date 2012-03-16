1. There are two folders in "code." One is for Backward Euler, the other one is for
Crank-Nicholson.
2. I didn't implement these two methods in one program, becausse I think this way
may be more convenient to run and compare.
3. In relax.c, three different relaxation schemes are implemented. You can simply
change the function called in relax.c to use either one of them.
4. Change the problem size in main.c.
5. After you run the program, the runtime will be displayed on the terminal.
(By default, it is using Gauss-Seidel with R&B with the problem size of 257.)
6. And an output file named soln.txt will be created, in which you could find the result
of the matrix.
7. Move that file to the same path with Analysis.m, then run the Matlab program, a plot
will be created. (Note that: the problem size should be consistent, and this has to be
done manually.)
8. The answers of all the questions could be found in Performance Analysis.pdf.