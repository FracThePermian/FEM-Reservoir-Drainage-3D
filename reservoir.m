function [para] = reservoir
 
para.NX = 100;                 % number of blocks
para.NY = 100;

para.N = para.NX *para.NY;
para.k = ones(para.N,1)*25;
para.phi = ones(para.N,1)*0.2;
para.h = 10;
para.W = 20000;
para.L = 20000;
para.ct = 1e-6;
para.mu = 1;
para.Bw = 1;
para.dx = para.L/para.NX;
para.dy = para.W/para.NY;
para.Ax = para.dy*para.h;
para.Ay = para.dx*para.h


global k phi
k = ones(para.N,1)*0.2;
phi = ones(para.N,1)*25
