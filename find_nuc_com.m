function [nuc_com,bw_nuc_com]=find_nuc_com(bw_nuc)

CC = bwconncomp(bw_nuc);
%[y_grid,x_grid]=meshgrid(1:size(bw_nuc,2),1:size(bw_nuc,1));
disty=bwdist(~bw_nuc);

nuc_com=zeros(length(size(bw_nuc)),length(CC.PixelIdxList));
for k=1:length(CC.PixelIdxList)
    [~,ind]=max(disty(CC.PixelIdxList{1,k}));
    [nuc_com(:,k)] = idx_num2vec(CC.PixelIdxList{1,k}(ind),size(bw_nuc));
%     nuc_com(1,k)=mean(x_grid(CC.PixelIdxList{1,k}));
%     nuc_com(2,k)=mean(y_grid(CC.PixelIdxList{1,k}));
end

    bw_nuc_com=zeros(size(bw_nuc));
    bw_nuc_com(idx_vec2num(nuc_com,size(bw_nuc)))=1;
    se = strel('disk',2);
    bw_nuc_com = imdilate(bw_nuc_com, se); 
    