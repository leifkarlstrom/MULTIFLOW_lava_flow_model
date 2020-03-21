% make synthethic cone 

nx=100
ny=100
dx=1 
ncone=1
H=0.5;
pad=0;
periodic=1 ;
N=1;
var=546;
beta=1+2*H;
window=0;

[M Pm fm Pv fv] = synthspecNEW(nx,ny,dx,H,pad,window,var);

[P, locations]= makecone(nx, ny, dx, ncone)

Z = M + P ;

ShadeMap(Z, dx, 'cone')