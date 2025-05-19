clear all
close all

dim=2;
ionicModelType=5; % 1 HH % 2 TT % 3 Paci % 4 Botti % 5 Amin
factorize=0;

if ionicModelType==1
    U_rest = -54.387;
    Tf=30;
    nt=5000;
    I_stim=280; 
    start_stim=5;
    stop_stim=5.3;
    Di=3.15e-4;   % S/cm
    De=1.35e-3;   % S/cm            
end
if ionicModelType==2
    U_rest = -85.23; 
    Tf=500;
    nt=5000;
    I_stim=280;
    start_stim=5;
    stop_stim=5.3;
    Di=3.15e-4;   % S/cm
    De=1.35e-3;   % S/cm
end
if ionicModelType==3
    U_rest = -0.0734525804324366; 
    Tf=1.0;
    nt=5000;
    I_stim= 280; 
    start_stim=5*1e-3;
    stop_stim=5.3*1e-3;
    Di=3.15e-1;
    De=1.35;
end
if ionicModelType==4
    U_rest = -0.0908810000000000; 
    Tf=1 ;
    nt=5000;
    I_stim= 280; 
    start_stim=5*1e-3;
    stop_stim=5.3*1e-3;
    Di=3.15e-1;
    De=1.35;
end
if ionicModelType==5
    U_rest = -0.0734525804324366; 
    Tf=0.1;
    nt=5000;
    I_stim= 280; 
    start_stim=5*1e-3;
    stop_stim=5.3*1e-3;
    Di=3.15e-1;
    De=1.35;
end

T=linspace(0,Tf,nt+1);

if dim==0
    ne=[];
    h=[];
    VstoreStep=1;
end

if dim==1
    Xf=1;
    nex=100;
    hx=Xf/100;
    ne=[nex];
    h=hx;
    VstoreStep=1;
end

if dim==2
    Xf=0.2;
    nex=100;
    hx=Xf/100;
    Yf=0.2;
    ney=100;
    hy=Yf/100;
    ne=[nex ney];
    h=[hx hy];
    VstoreStep=1;
end

if dim==3
    Xf=1;
    nex=100;
    hx=Xf/nex;
    Yf=1;
    ney=100;
    hy=Yf/ney;
    Zf=1;
    nez=100;
    hz=Zf/nez;
    ne=[nex ney nez];
    h=[hx hy hz];
    VstoreStep=1e16;
end

pg=PointGrid(ne+1,h);
[Xc,Yc,Zc]=pg.getCooCenters;

which=find(0.4<Xc & Xc<0.6 & 0.4<Yc);
o=ones(size(Xc));

o(which)=1.0;

L=assembleMatricesH(pg,'diff',o);
M=assembleMatricesH(pg,'mass');


%%%%%%%%%%%%%%%%%%%%%%%%%% 
% ELECTRODES %
%%%%%%%%%%%%%%%%%%%%%%%%%%
nk = 1;
el{1} = [0.1; 0.1];
for i = 1:nk
    Mk{i} = zeros(size(Xc));  % Inizializza matrice dell'elettrodo
    dist_centri = sqrt((Xc - el{i}(1)).^2 + (Yc - el{i}(2)).^2);
    [~, which_k] = min(dist_centri(:));
    Mk{i}(which_k) = 1.0;
end

Mtot = zeros(size(Mk{1}));
for k = 1:nk
    Mtot = Mtot + Mk{k};
end
exportVTK_MEA(pg, Mtot)






[X,Y,Z]=pg.getCoo;
Iapp=zeros(pg.get_nv,1);
which=find(X<0.02 & Y<0.02 & Z<0.02);
Iapp(which)=I_stim;

bd=Bidomain(pg,M,L,T,ionicModelType,factorize,U_rest,Di,De,VstoreStep);

% bd.V(:,1)=U_rest;

bd.IappStartTime=start_stim;
bd.IappStopTime=stop_stim;
bd.Iapp=Iapp;

bd.exportStep=50;

bd.run;

bd.plotAndAnalyze(); 

return
