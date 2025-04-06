function y = app_SP(x,L_S,Lii,LSi,LiS)
    somma=zeros(length(x),length(Lii)+1);
    for ids = 1:length(Lii)
        temp=LiS{ids}*x;
        temp=Lii{ids}\temp;
        somma(:,ids)=-LSi{ids}*temp;
    end
    somma(:,end)=L_S*x;
    y=sum(somma,2);

end