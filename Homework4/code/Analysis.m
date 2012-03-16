clear all;
load 'soln.txt';
n = 257;
x = soln(1:n,1:n);
surfc(x);