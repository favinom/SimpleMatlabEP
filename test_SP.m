clear all
close all

Xf=1;
nex=512;
x=linspace(0,Xf,nex+1)';
hx=diff(x);

Yf=1;
ney=512;
y=linspace(0,Yf,ney+1)';
hy=diff(y);

[X,Y]=ndgrid(x,y);

[M,L]=assembleMatrices_old(hx,hy);
L=L+0.001*M; 

node_id = 1:(nex+1)*(ney+1);

id{1} = find((X<Xf/2) & (Y<Yf/2));
id{2} = find((X>Xf/2) & (Y<Yf/2));
id{3} = find((X<Xf/2) & (Y>Yf/2));
id{4} = find((X>Xf/2) & (Y>Yf/2));
id_S=node_id;

for j = 1:length(id)
    id_S=setdiff(id_S,id{j});
end

L_S = L(id_S,id_S);
for ids = 1:length(id)
    Lii{ids} = L(id{ids},id{ids});
    H{ids}=chol(Lii{ids});
    LSi{ids} = L(id_S,id{ids});
    LiS{ids} = L(id{ids},id_S);
end


V=ones(nex+1,ney+1);
b=V(id_S)';

tic
%fcn =@(x) app_SP(x,L_S,Lii,LSi,LiS);
fcn_fact =@(x) app_SP_fact(x,L_S,H,LSi,LiS);
x = pcg(fcn_fact,b,1e-7,1000);
toc


bfull=V(:);

tic
xfull = pcg(L,bfull,1e-6,1000);
toc
