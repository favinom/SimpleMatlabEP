% no clear all
close all

[m,n,l]=size(V);

for i=1:10:l
    surf(X(:,:),Y(:,:),V(:,:,i))
    grid off
    %ma=max(max(V(:,:,1,i)));
    %mi=min(min(V(:,:,1,i)));
    xlabel('x')
    ylabel('y')
    %format long
    %ma-mi
    view(0,90)
    zlim([-80 120])
    title([num2str((i-1)*dt)])
    pause(0.000001)
end

