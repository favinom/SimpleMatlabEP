clear all
close all

%% 1D

%pg1=PointGrid([21],[0.1]);

%tic
%[M1,A1]=assembleMatrices(pg1);
%toc;

%tic
%A1N=assembleMatricesH(pg1,'diff');
%M1N=assembleMatricesH(pg1,'mass');
%toc;

%max(max(abs(M1-M1N)));
%max(max(abs(A1-A1N)));

%% 2D

%pg2=PointGrid([21,11],[0.1 0.4]);

%tic
%[M2,A2]=assembleMatrices(pg2);
%toc

%tic
%A2N=assembleMatricesH(pg2,'diff');
%M2N=assembleMatricesH(pg2,'mass');
%toc

%max(max(abs(M2-M2N)));
%max(max(abs(A2-A2N)));

%% 3D

%pg3=PointGrid([4,4,4],[0.1 0.2 0.1]);

%tic
%[M3,A3]=assembleMatrices(pg3);
%toc

%tic
%A3N=assembleMatricesH(pg3,'diff');
%M3N=assembleMatricesH(pg3,'mass');
%toc

%max(max(abs(M3-M3N)))
%max(max(abs(A3-A3N)))


%% how to use it
pg3=PointGrid([11,11,11],[0.1 0.1 0.1]);

[Xc,Yc,Zc]=pg3.getCooCenters;

c=ones(size(Xc));

which=find(Xc>0.5 & Yc>0.5);
c(which)=4;

A3N=assembleMatricesH(pg3,'diff',c);
M3N=assembleMatricesH(pg3,'mass',c);
