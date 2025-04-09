function plotSubid(pg,ids,idi)

[X,Y]=pg.getCoo;

plot(X(ids),Y(ids),'r*','MarkerSize',2);
hold on

for i=1:length(idi)
    plot(X(idi{i}),Y(idi{i}),'b*','MarkerSize',0.25);
end
