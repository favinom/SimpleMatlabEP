function [ids,idi,vertical_lines,horizontal_lines]=buildSubdomainIds(pg,nex,ney)

[X,Y]=pg.getCoo;

x=unique(X(:));
y=unique(Y(:));

xmin=min(X(:));
xmax=max(X(:));

ymin=min(Y(:));
ymax=max(Y(:));

xs=linspace(xmin,xmax,nex+1);
ys=linspace(ymin,ymax,ney+1);

xs=xs(2:end-1);
ys=ys(2:end-1);


for i=1:length(xs)
    [~,w]=min(abs(x-xs(i)));
    xs(i)=x(w);
end

for j=1:length(ys)
    [~,w]=min(abs(y-ys(j)));
    ys(j)=y(w);
end

xs=[xmin-1 xs xmax+1];
ys=[ymin-1 ys ymax+1];

counter=0;
for i=1:length(xs)-1
    for j=1:length(ys)-1
        counter=counter+1;
        xxmin=xs(i);
        xxmax=xs(i+1);
        yymin=ys(j);
        yymax=ys(j+1);

        idi{counter}=find(xxmin<X & X<xxmax & yymin<Y & Y<yymax);
    end
end

ids=1:numel(X);
for i=1:length(idi)
    ids=setdiff(ids,idi{i});
end


% Calcolo dei vettori dei nodi sulle interfacce verticali
for i = 1:nex-1
    xval = xs(i+1); % posizione della linea verticale
    vertical_lines{i} = find(abs(X - xval) < 1e-12); % tolleranza per confronto numerico
end

% Calcolo dei vettori dei nodi sulle interfacce orizzontali
for j = 1:ney-1
    yval = ys(j+1); % posizione della linea orizzontale
    horizontal_lines{j} = find(abs(Y - yval) < 1e-12); % tolleranza per confronto numerico
end

