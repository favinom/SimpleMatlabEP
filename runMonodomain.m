clear all
close all

dim=2;
ionicModelType=3; % 1 HH % 2 TT % 3 Paci
factorize=0;

if ionicModelType==1
    U_rest = -54.387;
    Tf=30*10;
    nt=5000*10;
    I_stim=280; 
    start_stim=5;
    stop_stim=5.3;
    diff=1e-3;
end
if ionicModelType==2
    U_rest = -85.23; 
    Tf=500;
    nt=5000;
    I_stim=280;
    start_stim=5;
    stop_stim=5.3;
    diff=1e-3;
end
if ionicModelType==3
    U_rest = -0.0734525804324366; 
    Tf=1.0;
    nt=5000;
    I_stim= 280; 
    start_stim=5*1e-3;
    stop_stim=5.3*1e-3;
    diff=1e0;
end

T=linspace(0,Tf,nt+1);

if dim==0
    ne=[];
    h=[];
end

if dim==1
    Xf=1;
    nex=100;
    hx=Xf/100;
    ne=[nex];
    h=hx;
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
end

pg=PointGrid(ne+1,h);
[M,L]=assembleMatrices(pg);

[X,Y,Z]=pg.getCoo;
Iapp=zeros(pg.get_nv,1);
which=find(X<0.2 & Y<0.2 & Z<0.2);
Iapp(which)=I_stim;

md=Monodomain(pg,M,L,T,ionicModelType,factorize,U_rest,diff);

md.V(:,1)=U_rest;

md.IappStartTime=start_stim;
md.IappStopTime=stop_stim;
md.Iapp=Iapp;

md.exportStep=100;

md.run;

return
