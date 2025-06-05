clear all
close all

addpath('../.')

dim=2;
ionicModelType=5; % 3 Paci % 4 Botti % 5 Amin

%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%% 
% SET IONIC MODELS, TIME STEP AND DIFFUSION PARAMETERS %
%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%% 

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
    Tf=0.6;
    nt=2000;
    I_stim= 280; 
    start_stim=5*1e-3;
    stop_stim=5.3*1e-3;
    %Di=3.15e-1; % mS/cm
    %De=1.35; % mS/cm
    %Di=1.3514;
    %De=0.31525;
    Di = 0.01;
    De = 0.1;
    %Di_test = [0.01, 0.05, 0.31525];
    %De_test = [0.1, 0.5, 1.3514];
end

T=linspace(0,Tf,nt+1);


%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%% 
% SET DIMENSIONS AND GEOMETRY %
%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%% 

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
    Xf=0.15;
    nex=150;
    hx=Xf/nex;
    Yf=0.15;
    ney=150;
    hy=Yf/ney;
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

which_isch=find(0.1<Xc & Xc<0.3 & 0.1<Yc & Yc<0.3);
isch=ones(size(Xc));

isch(which_isch)=1;

L=assembleMatricesH(pg,'diff',isch);
M=assembleMatricesH(pg,'mass');


%%%%%%%%%%%%%%%%%%%%%%%%%% 
% ELECTRODES %
%%%%%%%%%%%%%%%%%%%%%%%%%%
nk = 9;
el = cell(1, nk);

% center = [0.1; 0.1];
d = 200e-4;  % 0.02 cm
if nk==1
    center = [Xf/2; Yf/2];
    offsets = [0];
elseif nk==9
    center = [Xf/2; Yf/2];
    offsets = [-1 0 1] * d;
elseif nk == 64
    center = [0.2-d/2; 0.2-d/2];
    offsets = [-3 -2 -1 0 1 2 3 4] * d;
end

k = 1;

for i = 1:sqrt(nk)  % righe (y)
    for j = 1:sqrt(nk)  % colonne (x)
        dx = offsets(j);
        dy = offsets(i);  
        el{k} = center + [dx; dy];
        k = k + 1;
    end
end

for i = 1:nk
    dist_centri = sqrt((Xc - el{i}(1)).^2 + (Yc - el{i}(2)).^2);
    [~, which_k] = min(dist_centri(:));
    chi_el{i} = zeros(size(Xc));  % Inizializza matrice dell'elettrodo
    chi_el{i}(which_k) = 1.0;
    Mk{i}=assembleMatricesH(pg,'mass', chi_el{i});
    fk{i} = sum(Mk{i},2);
end

chi_tot = zeros(size(chi_el{1}));
for k = 1:nk
    chi_tot = chi_tot + chi_el{k};
end
exportVTK_MEA(pg, chi_tot)



%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%% 
% BOUNDARY NODES (3 edges) AND APPLIED CURRENT 
%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%% 

boundary_nodes = pg.getBoundary;

% plot_Boundary_and_electrodes
% pause

[X,Y,Z]=pg.getCoo;
Iapp=zeros(pg.get_nv,1);
which=find(X<0.02 & Y<0.02 & Z<0.02);
Iapp(which)=I_stim;


%%%%%%%%%%%%%%%%%%%%%%%%%%  
% RUN SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%  

mea=MEA(pg,M,L,T,ionicModelType,U_rest,Di,De,VstoreStep,fk,boundary_nodes);

mea.IappStartTime=start_stim;
mea.IappStopTime=stop_stim;
mea.Iapp=Iapp;

mea.exportStep=10;

mea.run;

plotAndAnalyze(mea,el); 

return
