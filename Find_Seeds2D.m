
function [seeds,seg_nuc]=Find_Seeds2D(nuc_org,Enmask,disk_radius)

want_to_seg=true;

% good parameter for Pauls measurements: diskradius=10, ws=20, bwopen=>0.02,0.05* ..

%disk_radius=15;   % typical feature size => mean parameter to play around with 
lower_bound_cell_area=5.5*disk_radius^2;  % estimate of smallest cell area in picture

nuc=nuc_org;%rgb2gray(nuc_org); % just copy

%% some filters to make the image nicer to work with

nuc=double(nuc)/max(double(nuc(:))); %normalize
nuc = adpmedian(nuc, 2*floor(disk_radius/2)+1);     % adaptiver Medianfilter gegen Salt&Pepper  
nuc=nuc-imopen(nuc,strel('disk',2*disk_radius));    % background correction - shouldnt be important for nucs     
nuc=double(nuc)/max(double(nuc(:))); % normalize


 ws=2*disk_radius;   % something to play with - size of average filter 
IM=averagefilter(nuc, [ws ws], 'replicate');
IM=imadjust((IM-nuc).*nuc); 

%% nucleus thresholding 
%empirically found (for stained nuclei)
seedsy=imopen(imclose(nuc>.35*graythresh(nuc) & adaptivethreshold(nuc,disk_radius*2)  & adaptivethreshold(nuc,(disk_radius*5)^2),strel('disk',2)) & IM ==0,strel('disk',round(disk_radius/10))) & Enmask; % & En < levelEn 

%% morphological cleaning
seedsy=imerode(seedsy,strel('disk',1));  % sorts out dirt 
seedsy=bwareaopen(seedsy,floor((1.7-sum(Enmask(:))/numel(Enmask))*0.04*disk_radius^2));
seedsy=imdilate(seedsy,strel('disk',floor(1)));
seedsy=bwareaopen(seedsy,floor(0.20*disk_radius^2));

holes=ones(size(seedsy))-seedsy;  
bigholes = bwareaopen(holes, floor(lower_bound_cell_area));
smallholes = holes & ~bigholes;
seedsy=seedsy | smallholes;

bw_nuc=seedsy;

%% center of nuclei 
% found as maxima in a distance matrix - good for separeting closeby nuclei

disty=bwdist(~seedsy);   % distance matric to find "center of mass" for nucs
disty=imgaussfilt(disty,disk_radius/4);  % an important parameter to play around with

msk = true(3,3,3);
msk(2,2,2) = false;
%# assign, to every voxel, the maximum of its neighbors
maxy = imdilate(disty,msk);
seeds = disty > maxy;
seeds = imdilate(seeds, strel('disk',floor(disk_radius/8)));   %% centers that are very close should fuse because they are one cell 
seeds = imdilate(seeds, strel('disk',floor(disk_radius/(5-4*(1-sum(Enmask(:))/numel(Enmask))))));
 [seeds_com,~]=find_nuc_com(seeds);
 seeds_com=idx_vec2num(seeds_com,[size(nuc_org,1),size(nuc_org,2)]); 
 seeds=zeros(size(seeds));
 seeds(seeds_com)=1;
 seeds = imdilate(seeds, strel('disk',floor(disk_radius/8)));

%figure; imshow(seeds)

%% segmentation of nuclei by watersheding around previously found seeds
if want_to_seg
    hy = fspecial('sobel');  %%typical water shedding https://de.mathworks.com/help/images/marker-controlled-watershed-segmentation.html
    hx = hy';
    Iy = imfilter(double(nuc)/256, hy, 'replicate');
    Ix = imfilter(double(nuc)/256, hx, 'replicate');
    gradmag = sqrt(Ix.^2 + Iy.^2);
    gradmag=gradmag/max(gradmag(:));
    gradmag2 = imimposemin(gradmag, bwmorph(~(Enmask & bw_nuc),'shrink','Inf')|  seeds);
    L = watershed(gradmag2);
I = zeros(size(gradmag2));
I(imdilate(L == 0, ones(1, 1))) = 1;

outlines= im2bw(I, 0.5);
seg_nuc=~outlines & bw_nuc;
 
%clean-up
CC_shapes=bwconncomp(seg_nuc);
shapes=zeros(size(nuc_org));
for c=2:length(CC_shapes.PixelIdxList)
    if ismember(1,ismember(seeds_com, CC_shapes.PixelIdxList{1,c}))==1 && length(CC_shapes.PixelIdxList{1,c}) < 100*disk_radius;
        shapes(CC_shapes.PixelIdxList{1,c})=1;
    end
end
seg_nuc=shapes;
end

%end
end