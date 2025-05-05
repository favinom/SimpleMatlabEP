clear all
close all

dim=2;
ionicModelType=1; % 1 Paci % 2 Botti
factorize=0;

Tf=1; %4;
dt=1e-4;
nt=Tf/dt;
%nt=5000;
I_stim= 280; 
start_stim=5*1e-3;
stop_stim=5.3*1e-3;
diff=1e0;

% 0 Botti 
% 1 paci

%U_rest0 = -0.0908810000000000; 
U_rest0 = -0.081651; 
U_rest1 = -0.0734525804324366; 


T=linspace(0,Tf,nt+1);

if dim==0
    ne=[];
    h=[];
    VstoreStep=1;
end

if dim==1
    Xf=1;
    nex=100;
    hx=Xf/nex;
    ne=[nex];
    h=hx;
    VstoreStep=1;
    nsd=1; 
end

if dim==2
    Xf=1;
    nex=100;
    hx=Xf/nex;
    Yf=1;
    ney=100;
    hy=Yf/ney;
    ne=[nex ney];
    h=[hx hy];
    VstoreStep=10;
    nsd=2; 
end

if dim==3
    Xf=1;
    nex=100;
    hx=Xf/nex;
    Yf=1;
    ney=100;
    hy=Yf/ney;
    Zf=1;
    nez=2;
    hz=Zf/nez;
    ne=[nex ney nez];
    h=[hx hy hz];
    VstoreStep=1e16;
    nsd=1; 
end


pg=PointGrid(ne+1,h);
[M,L]=assembleMatrices(pg);

%Scelta=zeros(pg.get_nv,1);
percentuale_zeri = 0.3; % <-- modifica questo valore tra 0 e 1
n_zeri = round(percentuale_zeri * pg.get_nv);
valori = [zeros(n_zeri,1); ones(pg.get_nv - n_zeri,1)];
Scelta = valori(randperm(pg.get_nv));

U_rest=Scelta .* U_rest1 + (1 - Scelta) .* U_rest0;

[X,Y,Z]=pg.getCoo;
Iapp=zeros(pg.get_nv,1);
which=find(X<0.1 & Y<0.1 & Z<0.1);
Iapp(which)=I_stim;
Scelta(which)=1;

md=Monodomain(pg,M,L,T,factorize,U_rest,diff,VstoreStep,Scelta,nsd);

h_piccolo = Xf/nex;
H_grande = Xf/nsd;

md.filename = sprintf('results_iter_%.4f_%.4f.csv', h_piccolo, H_grande);
fid = fopen(md.filename, 'w');

md.V(:,1)=U_rest;

md.IappStartTime=start_stim;
md.IappStopTime=stop_stim;
md.Iapp=Iapp;

md.exportStep=100;

md.run;

fclose(fid);

md.plotAndAnalyze(); 

return
