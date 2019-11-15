function idx_vec=idx_num2vec(idx_num,s)
idx_vec=squeeze(zeros([ndims(s),length(idx_num)]));
rest=idx_num;
for s1=length(s):-1:2
    hyparea=ones(length(idx_num),1);
    for s2=s1-1:-1:1
        hyparea=hyparea*s(s2);
    end
    idx_vec(s1,:)=ceil(rest./hyparea)';
    rest=mod(rest,hyparea);
    rest(rest==0)=hyparea(rest==0);
end
idx_vec(1,:)=rest;
end


