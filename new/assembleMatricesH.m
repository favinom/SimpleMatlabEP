function A=assembleMatricesH(pg,type,c)

dim=pg.dim;
ne=pg.get_ne;
conn=pg.get_connectivity;

if ~exist('c')
    c=ones(1,ne);
end

c=c(:)';

hx=pg.h(1);
hy=pg.h(2);
hz=pg.h(3);

Aref=[1 -1; -1 1];
Mref=[1/2 0; 0 1/2];

Ax=Aref/hx;
Ay=Aref/hy;
Az=Aref/hz;

Mx=Mref*hx;
My=Mref*hy;
Mz=Mref*hz;

if strcmp(type,'diff')

    if dim==1
        Aloc=Ax;
    elseif dim==2
        A2x=kron(My,Ax);
        A2y=kron(Ay,Mx);
        Aloc=A2x+A2y;
    elseif dim==3
        A3x=kron(Mz,kron(My,Ax));
        A3y=kron(Mz,kron(Ay,Mx));
        A3z=kron(Az,kron(My,Mx));
        Aloc=A3x+A3y+A3z;
    end

elseif strcmp(type,'mass')

    if dim==1
        Aloc=Mx;
    elseif dim==2
        Aloc=kron(My,Mx);
    elseif dim==3
        Aloc=kron(Mz,kron(My,Mx));
    end

else
    error('wrong string with matrix type')
end

ii=(1:size(Aloc,1))';
ii=repmat(ii,[1 size(Aloc,2)]);
jj=ii';

C=repmat(c,[length(Aloc(:)) 1]);
V=repmat(Aloc(:),[1,ne]);

V=C.*V;
I=conn(ii(:),:);
J=conn(jj(:),:);

A=sparse(I,J,V);
