%% 2D distribution?
clear all;
r = linspace(1,100,100);
rd = r;

[R, Rd] = meshgrid(r, rd);
N = Rd.^(R-1).*exp(-Rd);

surf(R, Rd, N)