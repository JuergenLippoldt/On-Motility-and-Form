function idx_num=idx_vec2num(idx_vec,s)
if ismatrix(idx_vec)
    idx_num=zeros(1,size(idx_vec,2));
else
    idx_num=0;
end
hyparea=ones(1,length(idx_num));
for s1=1:length(s)
    if s1>1
        hyparea=hyparea*s(s1-1);
        idx_num(:)=idx_num(1,:)+(idx_vec(s1,:)-1).*hyparea; 
    else
        idx_num(:)=idx_num(1,:)+idx_vec(s1,:).*hyparea;       
    end
end
end