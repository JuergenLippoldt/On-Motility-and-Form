function [shapes]=watershed_seeds2D(phaco_org,Enmask,seedsy,seg_nuc,disk_radius,dispersed,nuc)

% phaco_org can be an actin signal - needs to be boundary information

%disk_radius=10;              % radius of a disk used to blur. The number of pixels should roughly equal 10µm. (whole number)
smoothy=1;                    % how often do you want to delate and erode to smooth the resulting grid (whole number)
lower_bound_cell_area=floor(.7*disk_radius^2);    % cells are larger than .. in square pixel. (whole number)
for_real_phaco=false;
%verstreut=true;

%% adjust seeds
 [seeds_com,~]=find_nuc_com(seedsy);
 seeds_com=idx_vec2num(seeds_com,[size(phaco_org,1),size(phaco_org,2)]);
 seedsy=imerode(seg_nuc,strel('disk',1));
 seedsy_points=zeros(size(phaco_org));
 seedsy_points(seeds_com)=1;
 seedsy_points=imdilate(seedsy_points,strel('disk',1));
  seedsy=seedsy+seedsy_points;
  seedsy(seedsy>1)=1;


%% actin - watershed
%typical water shedding https://de.mathworks.com/help/images/marker-controlled-watershed-segmentation.html

act=phaco_org;%rgb2gray(phaco_org);
act=act-imopen(act,strel('disk',3*disk_radius));
act = imadjust(act);
act = adpmedian(act, 2*floor(.75*disk_radius)+1); 

hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(act-.5*nuc)/256, hy, 'replicate');
Ix = imfilter(double(act-.5*nuc)/256, hx, 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
gradmag=gradmag/max(gradmag(:));

    gradmag = (imgaussfilt(gradmag,1.5)+.2).*(act+.5);
    gradmag=gradmag/max(gradmag(:));

clear hy hx Iy Ix

% adjust Enmask if it is calculated from phaco 
if for_real_phaco
if sum(Enmask(:))/(size(Enmask,1)*size(Enmask,2)) < .9 && dispersed
    gradmag2 = imimposemin(gradmag, bwmorph(~Enmask,'shrink','Inf')|  seedsy); % use this if you estimate the background
else
    gradmag2 = imimposemin(gradmag, seedsy); % use this if you estimate the background
end
else
if sum(Enmask(:))/(size(Enmask,1)*size(Enmask,2)) < .9 && dispersed
    gradmag2 = imimposemin(act, bwmorph(~Enmask,'shrink','Inf')|  seedsy); % use this if you estimate the background
else
    gradmag2 = imimposemin(act, seedsy); % use this if you estimate the background
end
end


%gradmag2 = imimposemin(gradmag,  seedsy);
L = watershed(gradmag2);
I = zeros(size(gradmag2));
%I(imdilate(L == 0, ones(3, 3))) = 1;
I(L == 0)=1;

outlines= im2bw(I, 0.5);
if dispersed
    mask=imbinarize(imgaussfilt(act,disk_radius/8),0.05);
    mask=bwmorph(mask,'dilate');
    mask=bwmorph(mask,'erode');
    mask=bwareaopen(mask,20);
    mask=~bwareaopen(~mask,20);
    mask=mask & Enmask;
    outlines(~mask)=0;
    disty=bwdist(~mask);
    outlines(disty==1)=1;
end

%% smooth and adjust
outlines_smooth=outlines;
for i=1:smoothy
outlines_smooth=bwmorph(outlines_smooth,'dilate');
end
for i=1:smoothy
outlines_smooth=bwmorph(outlines_smooth,'erode');
end
%imshow(outlines_smooth)

holes=ones(size(outlines_smooth))-outlines_smooth;
bigholes = bwareaopen(holes, lower_bound_cell_area);
smallholes = holes & ~bigholes;
outlines_smooth=outlines_smooth | smallholes;

%imshow(outlines_smooth)

clear holes smallholes bigholes i I I2 L

I2 = zeros(size(gradmag2));
I2(outlines_smooth == 1) = 1;
 
% discard unrealistic shapes (better have less info than faulty one)
CC_shapes=bwconncomp(~I2,4);
stats = regionprops(CC_shapes,'Area','Perimeter');
for i=1:length(stats)
    shape_para(i)=stats(i).Perimeter/stats(i).Area^.5;
end

shapes=zeros(size(act));
for c=2:length(CC_shapes.PixelIdxList)
    if ismember(1,ismember(seeds_com, CC_shapes.PixelIdxList{1,c}))==1 && length(CC_shapes.PixelIdxList{1,c}) > (.6*disk_radius)^2  && length(CC_shapes.PixelIdxList{1,c}) < (3.5*disk_radius)^2 && shape_para(c) < 8
        
        shapes(CC_shapes.PixelIdxList{1,c})=1;
    end
end

end