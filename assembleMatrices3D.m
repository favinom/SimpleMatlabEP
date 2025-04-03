function [M,L]=assembleMatrices3D(hx,hy,hz)

hxz=[hx;0];
zhx=[0;hx];
nvx=length(hx)+1;

%mx=spdiags( [1/6*hxz 1/3*(hxz+zhx) 1/6*zhx], [-1 0 1], nvx,nvx);
mx=spdiags( [0*hxz 1/2*(hxz+zhx) 0*zhx], [-1 0 1], nvx,nvx);

hyz=[hy;0];
zhy=[0;hy];
nvy=length(hy)+1;

%my=spdiags( [1/6*hyz 1/3*(hyz+zhy) 1/6*zhy], [-1 0 1], nvy,nvy);
my=spdiags( [0*hyz 1/2*(hyz+zhy) 0*zhy], [-1 0 1], nvy,nvy);

hzz=[hz;0];
zhz=[0;hz];
nvz=length(hz)+1;

%my=spdiags( [1/6*hyz 1/3*(hyz+zhy) 1/6*zhy], [-1 0 1], nvy,nvy);
mz=spdiags( [0*hzz 1/2*(hzz+zhz) 0*zhz], [-1 0 1], nvz,nvz);


hxz=[1./hx;0];
zhx=[0;1./hx];

ax=spdiags( [-hxz (hxz+zhx) -zhx], [-1 0 1], nvx,nvx);

hyz=[1./hy;0];
zhy=[0;1./hy];

ay=spdiags( [-hyz (hyz+zhy) -zhy], [-1 0 1], nvy,nvy);

hzz=[1./hz;0];
zhz=[0;1./hz];

az=spdiags( [-hzz (hzz+zhz) -zhz], [-1 0 1], nvz,nvz);

M=kron(mz,kron(my,mx));
L=kron(mz,kron(my,ax))+kron(mz,kron(ay,mx))+kron(az,kron(my,mx));
