clear all
close all

dim=2;
ionicModelType=4; % 1 HH % 2 TT
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
    Xf=1;
    nex=100;
    hx=Xf/100;
    Yf=1;
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
[M,L]=assembleMatrices(pg);

[X,Y,Z]=pg.getCoo;
Iapp=zeros(pg.get_nv,1);
which=find(X<0.2 & Y<0.2 & Z<0.2);
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
