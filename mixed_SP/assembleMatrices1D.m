function [M,L]=assembleMatrices1D(pg)

hx=pg.get_h(1);
hx=hx(:);

hxz=[hx;0];
zhx=[0;hx];
nvx=length(hx)+1;

%mx=spdiags( [1/6*hxz 1/3*(hxz+zhx) 1/6*zhx], [-1 0 1], nvx,nvx);
M=spdiags( [0*hxz 1/2*(hxz+zhx) 0*zhx], [-1 0 1], nvx,nvx);

hxz=[1./hx;0];
zhx=[0;1./hx];

L=spdiags( [-hxz (hxz+zhx) -zhx], [-1 0 1], nvx,nvx);


