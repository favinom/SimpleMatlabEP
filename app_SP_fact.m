function y = app_SP_fact(x,L_S,H,LSi,LiS)
    somma=zeros(length(x),length(H)+1);
    for ids = 1:length(H)
        temp=LiS{ids}*x;
        temp2=H{ids}'\temp;
        temp=H{ids}\temp2;
        somma(:,ids)=-LSi{ids}*temp;
    end
    somma(:,end)=L_S*x;
    y=sum(somma,2);

end