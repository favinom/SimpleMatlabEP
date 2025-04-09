function plotSubid(pg,ids,idi)

[X,Y]=pg.getCoo;

plot(X(ids),Y(ids),'r*');
hold on

for i=1:length(idi)
    plot(X(idi{i}),Y(idi{i}),'b*');
end
