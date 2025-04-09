clear all
close all

dim=2;
ionicModelType=2; % 1 HH % 2 TT
factorize=0;

if ionicModelType==1
    U_rest = -54.387;
    Tf=30;%*10;
    nt=5000;%*10;
end
if ionicModelType==2
    U_rest = -85.23; 
    Tf=500;
    nt=5000;
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
Iapp(which)=1;

bd=Bidomain(pg,M,L,T,ionicModelType,factorize,U_rest);

bd.V(:,1)=U_rest;

bd.IappStartTime=5;
bd.IappStopTime=5.3;
bd.Iapp=Iapp;

bd.exportStep=100;

bd.run;

return
