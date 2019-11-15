
% folder with images
pathy_main='E:\Friedl_19_02_13\shapes\_for elongation factor\NT shRNA\6 mg_ml';
oldp=cd(pathy_main);
filies=dir('*tif');

disk_radius=15;  % estimation of typical feature size in pixel - adjust for your images!


for lauf=1:length(filies)
    clearvars -except lauf filies pathy_main oldp disk_radius
    disp(filies(lauf).name);  % Just to get an estimation how long the script will be running

    pathy=[pathy_main,'\',filies(lauf).name];

%% load images

FileTif=pathy;
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
 
FinalImage=zeros(nImage,mImage,NumberImages,'uint16');
for i=1:NumberImages
   FinalImage(:,:,i)=imread(FileTif,'Index',i);
end

%% initial noise reduction and preparetion

    nuc_org=squeeze(FinalImage(:,:,1));%imread([pathy,'\',nuc_dir(f).name]);
    nuc_org=double(nuc_org)/double(max(nuc_org(:)));
%    nuc_org=imresize(imread([pathy,'\',nuc_dir(f).name]),.5);
%     demix_org=double(nuc_org)/double(max(nuc_org(:)));
    act_org=squeeze(FinalImage(:,:,2));%imread([pathy,'\',nuc_dir(f).name]);
    act_org=double(act_org)/double(max(act_org(:)));
    
    nuc=nuc_org;%rgb2gray(nuc_org); % just copy
    
    intensity_adjustment=imgaussfilt(nuc_org+act_org,floor(.7*disk_radius));
    intensity_adjustment=1./(intensity_adjustment/(max(intensity_adjustment(:))));
    
    nuc=double(nuc)/max(double(nuc(:))); %normalize
    nuc = adpmedian(nuc, 2*floor(disk_radius/2)+1);     % adaptiver Medianfilter gegen Salt&Pepper  
    nuc=nuc+.7*nuc.*intensity_adjustment;    % background correction - shouldnt be important for nucs     
    nuc=double(nuc)/max(double(nuc(:))); % normalize

    act=act_org;%rgb2gray(nuc_org); % just copy    
    act=double(act)/max(double(act(:))); %normalize
    act= adpmedian(act, 2*floor(disk_radius/2)+1);     % adaptiver Medianfilter gegen Salt&Pepper  
    act=act+.7*act.*intensity_adjustment;  %-imopen(act,strel('disk',2*disk_radius));    % background correction - shouldnt be important for nucs     
    act=double(act)/max(double(act(:))); % normalize

    %% finding cell free regions low in intensity and entropy
    %intensity
    int_mask=imgaussfilt(nuc+act,floor(.5*disk_radius));

    nuc2=nuc+act;
    fgm= ones(size(nuc));
    for ws=1:10
        fgm = min(fgm,adaptivethreshold(nuc2,2*ws^2*disk_radius));%im2bw(nuc(:,:,i),blub);%imbinarize(nuc(:,:,i));
    end
    
int_mask=int_mask/max(int_mask(:));
%int_mask = reshape(double(int_mask(:))/max(double(int_mask(:))),size(int_mask));
int_mask = imdilate(int_mask >graythresh(.65*int_mask(int_mask>0)),strel('disk',2));
int_mask = imerode(int_mask,strel('disk',2));
int_mask = bwareaopen(int_mask,floor((.5*disk_radius)^3));

% entropy

Eim = entropyfilt((act + nuc)/2);
%G = fspecial('gaussian',[10 10],3);
Enmask=zeros(size(Eim));

    Eim=Eim/max(Eim(:));
Enmask= im2bw(Eim,graythresh(Eim(mean(Eim)>0.05))); 
%Enmask=imfill(Enmask);
holes=imfill(Enmask,'holes')-Enmask;
Enmask=Enmask+ holes;

Enmask_without_int=Enmask;
Enmask=Enmask.*int_mask;
    
%% find nucleus positions 
    [seeds,seg_nuc]=Find_Seeds2D(nuc_org,Enmask,disk_radius );
    overlay_nuc=bsxfun(@plus,nuc_org,double(seeds));
 %   imshow(overlay_nuc)
    
 %% save nuc positions
    if true
    	imwrite(255*uint8(seeds),[filies(lauf).name(1:end-4),'seeds_dr15.png'] );
%    	imwrite(255*uint8(seg_nuc),fullfile(pathy,'seg_nuc.png') );
    	imwrite(overlay_nuc,[filies(lauf).name(1:end-4),'overlay_nuc_dr15.png'] );
    end
    
 %% segmentate cells by watershed
	if  all(size(act)==size(nuc))
          
      
        shapes=watershed_seeds2D(act,Enmask,seeds,seg_nuc,disk_radius,true,nuc);
        
        disty=bwdist(~shapes);
        outlines=disty==1;
        
  %% compose overlay and save segmentation
        overlay=zeros(size(act,1),size(act,2),3);
%overlay_nuc=zeros(size(act,1),size(act,2),size(act,3),3);
overlay(:,:,1)=nuc;
overlay(:,:,2)=act;%/255;
 imwrite(overlay,[filies(lauf).name(1:end-4),'act_nuc.png'] );
overlay(:,:,3)=outlines;%/255;
imwrite(overlay,[filies(lauf).name(1:end-4),'_overlay_shapes.png'] );
% figure;
% imshow(overlay);

%        figure; imshow(overlay)
imwrite(shapes,[filies(lauf).name(1:end-4),'shapes.png'] );
%         imwrite(act + nuc,[filies(lauf).name(1:end-4),'act_nuc.png'] );
%  	    imwrite(overlay,[filies(lauf).name(1:end-4),'_overlay_shapes.png'] );
	end	
% end
%waitbar( lauf/length(filies));
end
% end
%close(waitbar(0));

%% scripts for plotting 
% commented out for copy + paste
if false
h=[.25, .5, .25; .5, 1, .5; .25, .5, .25]/2;

for f=1:length(nuc_dir)
    tracksy_org=imread([pathy,'\',tracksy_dir(f).name]);
    tracks=zeros(size(tracksy_org));
    for c=1:2
        tracks(:,:,c)=imfilter(tracksy_org(:,:,c),h);
        tracks(:,:,c)=tracks(:,:,c)/max(max(tracks(:,:,c)));        
    end
    imwrite(tracks,fullfile(pathy,'tracks_on_black.png') );
	nuc_org=imread([pathy,'\',nuc_dir(f).name]);%ControlDemixNucleus(f).cdata;%%squeeze(Nuc_video(1:end-60,:,:,f));     
    phaco_org=imread([pathy,'\',phaco_dir(f).name]);
    imwrite(tracks+double(nuc_org(100:100+500-1,1600:1600+1800-1,:)+phaco_org(100:100+500-1,1600:1600+1800-1,:))/255,fullfile(pathy,'tracks_on nuc.png') );
end
end
% end

% resave
% FileTif='de_bw.tif';
% InfoImage=imfinfo(FileTif);
% mImage=InfoImage(1).Width;
% nImage=InfoImage(1).Height;
% NumberImages=length(InfoImage);
% 
% for i=1:NumberImages
%    FinalImage=imread(FileTif,'Index',i);
%    de_bw(:,:,i)=logical(FinalImage);
% %   imwrite(uint8(FinalImage)*255,'seeds2.tiff' );
% end