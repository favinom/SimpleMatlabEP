acti = -200*ones(size(V(:,:,1)));
repo = -200*ones(size(acti));

v_up   = -50;            % threshold for acti
v_down = -45;          % 0.9*(-84);    %threshold for repo

i=1;
id_acti=find((V(:,:,1)>v_up));
time_interp=(T(i)).*(v_up-U_rest)./(V(:,:,i)-U_rest);
acti(id_acti)=time_interp(id_acti);

for i=2:length(T)
    id_acti=find((V(:,:,i)>v_up) & (V(:,:,i-1)<v_up));
    time_interp=T(i-1)+(T(i)-T(i-1)).*(v_up-V(:,:,i-1))./(V(:,:,i)-V(:,:,i-1));
    acti(id_acti)=time_interp(id_acti);
    id_repo=find((V(:,:,i)<v_down) & (V(:,:,i-1)>v_down));
    time_interp_repo=T(i-1)+(T(i)-T(i-1)).*(v_down-V(:,:,i-1))./(V(:,:,i)-V(:,:,i-1));
    repo(id_repo)=time_interp_repo(id_repo);
end

figure
surf(X(:,:),Y(:,:),acti(:,:))
contourf(X, Y, acti, 20) % 20 isochrones
colorbar
xlabel('x')
ylabel('y')
view(0,90)
title('acti')

figure
surf(X(:,:),Y(:,:),repo(:,:))
contourf(X, Y, repo, 20) % 20 isochrones
colorbar
xlabel('x')
ylabel('y')
view(0,90)
title('repo')

