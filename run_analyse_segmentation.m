
%% A script that calculates, organizes and plots parameters of logical maps of cell segmentations
% Statistical functions and plots are commented out (if false) at the
% bottom of the script for easy adjustments + copy and paste

folder='E:\Friedl_19_02_13\shapes\analysis\all';
cd(folder)
disk_radius=15;

filies=dir('*shapes*');
shape_measures_all=[];

for f=1:length(filies)
    shapes_org=logical(imread(filies(f).name));
    boarder=true(size(shapes_org));
    boarder(2:end-1,2:end-1)=0;
    
    outside=imdilate(shapes_org+boarder,strel('disk',round(1.3*disk_radius)));
    outside=imerode(outside,strel('disk',round(1.4*disk_radius)));
    
    outer_region=imdilate(~outside,strel('disk',round(12*disk_radius)));
    
    CC_shapes=bwconncomp(shapes_org,6);
    
    shape_inside=zeros(length(CC_shapes.PixelIdxList) ,1);
    for c=1:length(CC_shapes.PixelIdxList) 
        if sum( outer_region(CC_shapes.PixelIdxList{c}))==0 
            shape_inside(c)=1;
        end        
    end
    
    shape_measures=shape_inside;
    shape_m=regionprops(CC_shapes,'Area','Perimeter','MajorAxisLength','MinorAxisLength','Orientation');
    
    neighbour_struck={};
    neighbour_struck_further={};
    for c=1:length(CC_shapes.PixelIdxList) 
        nn=zeros(length(CC_shapes.PixelIdxList) ,1);
        nnf=zeros(length(CC_shapes.PixelIdxList) ,1);
        large_cell=zeros(size(shapes_org));
        large_cell(CC_shapes.PixelIdxList{c})=1;
        large_cell=imdilate(large_cell,strel('disk',7));
        larger_cell=zeros(size(shapes_org));
        larger_cell(CC_shapes.PixelIdxList{c})=1;
        larger_cell=imdilate(larger_cell,strel('disk',round(2*disk_radius)));
        for z=1:length(CC_shapes.PixelIdxList) 
            if  sum( large_cell(CC_shapes.PixelIdxList{z}))>0
                nn(z)=1;
            end        
            if  sum( larger_cell(CC_shapes.PixelIdxList{z}))>0
                nnf(z)=1;
            end 
        end
        nn(c)=0;
        nnf(c)=0;
        neighbour_struck{c}=find(nn);     
        neighbour_struck_further{c}=find(nnf);   
    end
    
    
     for c=1:length(CC_shapes.PixelIdxList) 
         shape_measures(c,2)=shape_m(c).Area;
         shape_measures(c,3)=shape_m(c).Perimeter/shape_m(c).Area^.5;
         shape_measures(c,4)=shape_m(c).MajorAxisLength/shape_m(c).MinorAxisLength;
         if length(neighbour_struck{c})<4
            shape_measures(c,5)=NaN;
         else
            mor=mean(arrayfun(@(x) shape_m(x).Orientation,neighbour_struck{c}));
            shape_measures(c,5)=(2*((cosd(shape_m(c).Orientation-mor))^2)-1);
         end
         if length(neighbour_struck_further{c})<8
            shape_measures(c,6)=NaN;
         else
            mor=mean(arrayfun(@(x) shape_m(x).Orientation,neighbour_struck_further{c}));
            shape_measures(c,6)=(2*((cosd(shape_m(c).Orientation-mor))^2)-1);
         end
     end
    shape_measures_all{f}=shape_measures;
    
end


if false 
    
    
	relevant_shapes=[];
    for i=1:length(shape_measures_all)

        not_secound=logical([0,1,0,1,1,1,1,1,0,1,0,1,1,1,1,1,1,1,0,1,0,1,1,1,0,1,0,1,1,1,1,1]); % there where multiple maps for some experiments because of a manual correction to check the validity of the segmentation
        if contains(filies(i).name,'CDH1') & contains(filies(i).name,'bottom') & contains(filies(i).name,'6_mg') & not_secound(i)
            relevant_shapes=[relevant_shapes; shape_measures_all{i}];
        end
        
    end
    if false
        relevant_mean_AR=nanmean(relevant_shapes_CDH1_up(:,4));
        [h,p]=kstest2(relevant_shapes_NT_bottom(:,4),relevant_shapes_CDH1_bottom(:,4))
    end
    
    relevant_mean_shape=nanmean(relevant_shapes(:,3));
    relevant_mean_shape_para_outside=nanmean(relevant_shapes(relevant_shapes(:,1)==0,3));
    relevant_mean_shape_para_inside=nanmean(relevant_shapes(relevant_shapes(:,1)==1,3));
    relevant_std_shape_para_outside=nanstd(relevant_shapes(relevant_shapes(:,1)==0,3));
    relevant_std_shape_para_inside=nanstd(relevant_shapes(relevant_shapes(:,1)==1,3));
    relevant_std_shape_para=nanstd(relevant_shapes(:,3)); 
    relevant_mean_AR_outside=nanmean(relevant_shapes(relevant_shapes(:,1)==0,4));
    relevant_mean_AR_inside=nanmean(relevant_shapes(relevant_shapes(:,1)==1,4));
    relevant_mean_AR=nanmean(relevant_shapes(:,4));
    relevant_std_AR_outside=nanstd(relevant_shapes(relevant_shapes(:,1)==0,4));
    relevant_std_AR_inside=nanstd(relevant_shapes(relevant_shapes(:,1)==1,4));
    relevant_std_AR=nanstd(relevant_shapes(:,4));
    relevant_mean_Orientation_outside=nanmean(relevant_shapes(relevant_shapes(:,1)==0,5));
    relevant_mean_Orientation_inside=nanmean(relevant_shapes(relevant_shapes(:,1)==1,5));
    relevant_mean_Orientation=nanmean(relevant_shapes(:,5));
    relevant_std_Orientation_outside=nanstd(relevant_shapes(relevant_shapes(:,1)==0,5));
    relevant_std_Orientation_inside=nanstd(relevant_shapes(relevant_shapes(:,1)==1,5));
    relevant_std_Orientation=nanstd(relevant_shapes(:,5));


    
    for i=1:length(shape_measures_all)
        mean_shape_para_outside(i)=nanmean(shape_measures_all{i}(shape_measures_all{i}(:,1)==0,3));
        mean_shape_para_inside(i)=nanmean(shape_measures_all{i}(shape_measures_all{i}(:,1)==1,3));
        mean_shape_para(i)=nanmean(shape_measures_all{i}(:,3));
        std_shape_para_outside(i)=nanstd(shape_measures_all{i}(shape_measures_all{i}(:,1)==0,3));
        std_shape_para_inside(i)=nanstd(shape_measures_all{i}(shape_measures_all{i}(:,1)==1,3));
        std_shape_para(i)=nanstd(shape_measures_all{i}(:,3));
        median_shape_para_outside(i)=nanmedian(shape_measures_all{i}(shape_measures_all{i}(:,1)==0,3));
        median_shape_para_inside(i)=nanmedian(shape_measures_all{i}(shape_measures_all{i}(:,1)==1,3));
        median_shape_para(i)=nanmedian(shape_measures_all{i}(:,3));
        mad_shape_para_outside(i)=mad(shape_measures_all{i}(shape_measures_all{i}(:,1)==0,3));
        mad_shape_para_inside(i)=mad(shape_measures_all{i}(shape_measures_all{i}(:,1)==1,3));
        mad_shape_para(i)=mad(shape_measures_all{i}(:,3));
        
        mean_AR_outside(i)=nanmean(shape_measures_all{i}(shape_measures_all{i}(:,1)==0,4));
        mean_AR_inside(i)=nanmean(shape_measures_all{i}(shape_measures_all{i}(:,1)==1,4));
        mean_AR(i)=nanmean(shape_measures_all{i}(:,4));
        std_AR_outside(i)=nanstd(shape_measures_all{i}(shape_measures_all{i}(:,1)==0,4));
        std_AR_inside(i)=nanstd(shape_measures_all{i}(shape_measures_all{i}(:,1)==1,4));
        std_AR(i)=nanstd(shape_measures_all{i}(:,4));
        median_AR_outside(i)=nanmedian(shape_measures_all{i}(shape_measures_all{i}(:,1)==0,4));
        median_AR_inside(i)=nanmedian(shape_measures_all{i}(shape_measures_all{i}(:,1)==1,4));
        median_AR(i)=nanmedian(shape_measures_all{i}(:,4));
        mad_AR_outside(i)=mad(shape_measures_all{i}(shape_measures_all{i}(:,1)==0,4));
        mad_AR_inside(i)=mad(shape_measures_all{i}(shape_measures_all{i}(:,1)==1,4));
        mad_AR(i)=mad(shape_measures_all{i}(:,4));
        
        mean_Orientation_outside(i)=nanmean(shape_measures_all{i}(shape_measures_all{i}(:,1)==0,5));
        mean_Orientation_inside(i)=nanmean(shape_measures_all{i}(shape_measures_all{i}(:,1)==1,5));
        mean_Orientation(i)=nanmean(shape_measures_all{i}(:,5));
        std_Orientation_outside(i)=nanstd(shape_measures_all{i}(shape_measures_all{i}(:,1)==0,5));
        std_Orientation_inside(i)=nanstd(shape_measures_all{i}(shape_measures_all{i}(:,1)==1,5));
        std_Orientation(i)=nanstd(shape_measures_all{i}(:,5));
        median_Orientation_outside(i)=nanmedian(shape_measures_all{i}(shape_measures_all{i}(:,1)==0,5));
        median_Orientation_inside(i)=nanmedian(shape_measures_all{i}(shape_measures_all{i}(:,1)==1,5));
        median_Orientation(i)=nanmedian(shape_measures_all{i}(:,5));
        mad_Orientation_outside(i)=mad(shape_measures_all{i}(shape_measures_all{i}(:,1)==0,5));
        mad_Orientation_inside(i)=mad(shape_measures_all{i}(shape_measures_all{i}(:,1)==1,5));
        mad_Orientation(i)=mad(shape_measures_all{i}(:,5));
        
        mean_Orientationlarge_outside(i)=nanmean(shape_measures_all{i}(shape_measures_all{i}(:,1)==0,6));
        mean_Orientationlarge_inside(i)=nanmean(shape_measures_all{i}(shape_measures_all{i}(:,1)==1,6));
        mean_Orientationlarge(i)=nanmean(shape_measures_all{i}(:,6));
        std_Orientationlarge_outside(i)=nanstd(shape_measures_all{i}(shape_measures_all{i}(:,1)==0,6));
        std_Orientationlarge_inside(i)=nanstd(shape_measures_all{i}(shape_measures_all{i}(:,1)==1,6));
        std_Orientationlarge(i)=nanstd(shape_measures_all{i}(:,6));
        median_Orientationlarge_outside(i)=nanmedian(shape_measures_all{i}(shape_measures_all{i}(:,1)==0,6));
        median_Orientationlarge_inside(i)=nanmedian(shape_measures_all{i}(shape_measures_all{i}(:,1)==1,6));
        median_Orientationlarge(i)=nanmedian(shape_measures_all{i}(:,6));
        mad_Orientationlarge_outside(i)=mad(shape_measures_all{i}(shape_measures_all{i}(:,1)==0,6));
        mad_Orientationlarge_inside(i)=mad(shape_measures_all{i}(shape_measures_all{i}(:,1)==1,6));
        mad_Orientationlarge(i)=mad(shape_measures_all{i}(:,6));
        
        not_secound=logical([0,1,0,1,1,1,1,1,0,1,0,1,1,1,1,1,1,1,0,1,0,1,1,1,0,1,0,1,1,1,1,1]);
        
%         
% p{1}=plot(((0:69)/4), CD6_core,'LineWidth',3,'Color',[8,48,107]/255);
% hold on
% p{2}=plot(((0:69)/4), CD2_core,'LineWidth',3,'Color',[33,113,181]/255);
% hold on
% p{3}=plot(((0:69)/4), NG6_core,'LineWidth',3,'Color',[0,109,44]/255);
    end
    
    
if false 
    AR_NT=[];
	AR_CDH1=[];
for f=1:length(not_secound)
	NTorCDH1(f)=contains(filies(f).name,'NT');
    bottomorup(f)=contains(filies(f).name,'bottom');
    densly(f)=contains(filies(f).name,'6_');
    if densly(f)
    if bottomorup(f)
        if NTorCDH1(f)
            AR_NT=[AR_NT; shape_measures_all{f}(shape_measures_all{f}(:,1)==0,4)];
        else
            AR_CDH1=[AR_CDH1; shape_measures_all{f}(shape_measures_all{f}(:,1)==0,4)];            
        end
    end
    end
end
  
        s1=scatter(mean_AR_outside(not_secound & densly & NTorCDH1 & bottomorup),std_AR_outside((not_secound)& densly & NTorCDH1 & bottomorup));
        s1.LineWidth = 0.6;
        s1.MarkerEdgeColor = [0,109,44]/255;
        s1.MarkerFaceColor = [0,109,44]/255;
            
%         s11=scatter(mean_AR_outside(not_secound & NTorCDH1 & ~bottomorup),std_AR_outside((not_secound)& NTorCDH1 & ~bottomorup));
%         s11.LineWidth = 0.6;
%         s11.MarkerEdgeColor = [80,189,88]/255;
%         s11.MarkerFaceColor = [80,189,88]/255;

        hold on
        s2=scatter(mean_AR_outside(not_secound & densly & ~NTorCDH1 & bottomorup),std_AR_outside((not_secound)& densly & ~NTorCDH1 & bottomorup));
        s2.LineWidth = 0.6;       
        s2.MarkerEdgeColor = [8,48,107]/255;
        s2.MarkerFaceColor = [8,48,107]/255;
 
%         s21=scatter(mean_AR_outside(not_secound & ~NTorCDH1 & ~bottomorup),std_AR_outside((not_secound)& ~NTorCDH1 & ~bottomorup));
%         s21.LineWidth = 0.6;    
%         s21.MarkerEdgeColor = [48,88,187]/255;
%         s21.MarkerFaceColor = [48,88,187]/255;

%         hold on
%         s3=scatter(mean_AR_inside(not_secound & NTorCDH1 & bottomorup),std_AR_inside((not_secound) & NTorCDH1 & bottomorup));
%         s3.LineWidth = 0.6;
%         s3.MarkerEdgeColor = [0,109,44]/255;
%         s3.MarkerFaceColor = [0,109,44]/255;
% 
%         s31=scatter(mean_AR_inside(not_secound & NTorCDH1 & ~bottomorup),std_AR_inside((not_secound) & NTorCDH1 & ~bottomorup));
%         s31.LineWidth = 0.6;    
%         s31.MarkerEdgeColor = [80,189,88]/255;
%         s31.MarkerFaceColor = [80,189,88]/255;
% 
%         hold on
%         s4=scatter(mean_AR_inside(not_secound & ~NTorCDH1 & bottomorup),std_AR_inside((not_secound)& ~NTorCDH1 & bottomorup));
%         s4.LineWidth = 0.6;
%         s4.MarkerEdgeColor = [8,48,107]/255;
%         s4.MarkerFaceColor = [8,48,107]/255;
%         
%         s41=scatter(mean_AR_inside(not_secound & ~NTorCDH1 &~bottomorup),std_AR_inside((not_secound)& ~NTorCDH1 &~bottomorup));
%         s41.LineWidth = 0.6;
%         s41.MarkerEdgeColor = [48,88,187]/255;
%         s41.MarkerFaceColor = [48,88,187]/255;

         hold on 
        p5=plot(1.9:0.05:2.85,0.808*(1.9:0.05:2.85)- 0.85,'Color','black','LineWidth',2 ); % s.d.(AR) = 0.808AR? 0.85 Atia et al
        legend([s1, s2, p5],{'shNT ','shCDH1 ','predicted relation'})
        %legend([s3, s4, p5],{'shNT - bottom layer','shCDH1 - bottom layer','predicted relation'})
        


end
     my_boxplot(shape_measures_org(shape_measures_org(:,1)==0,5),shape_measures_org(shape_measures_org(:,1)==1,5),shape_measures_corrected(shape_measures_corrected(:,1)==0,5),shape_measures_corrected(shape_measures_corrected(:,1)==1,5),'NT_2mg/ml_outside raw','NT_2mg/ml_inside raw', 'NT_2mg/ml_outside corrected','NT_2mg/ml_inside corrected');

    
	my_boxplot(shape_measures_all_CDH1_2_bottom(shape_measures_all_CDH1_2_bottom(:,1)==0,4),shape_measures_all_CDH1_2_bottom(shape_measures_all_CDH1_2_bottom(:,1)==1,4),shape_measures_all_CDH1_6_bottom(shape_measures_all_CDH1_6_bottom(:,1)==0,4),shape_measures_all_CDH1_6_bottom(shape_measures_all_CDH1_6_bottom(:,1)==1,4),shape_measures_all_NT2_bottom(shape_measures_all_NT2_bottom(:,1)==0,4),shape_measures_all_NT2_bottom(shape_measures_all_NT2_bottom(:,1)==1,4),shape_measures_all_NT6_bottom(shape_measures_all_NT6_bottom(:,1)==0,4),shape_measures_all_NT6_bottom(shape_measures_all_NT6_bottom(:,1)==1,4),'CDH1_2mg/ml_outside','CDH1_2mg/ml_inside','CDH1_6mg/ml_outside','CDH1_6mg/ml_inside','NT_2mg/ml_outside','NT_2mg/ml_inside','NT_6mg/ml_outside','NT_6mg/ml_inside')
     
	my_boxplot(shape_measures_all_CDH1_2_bottom(shape_measures_all_CDH1_2_bottom(:,1)==0,3),shape_measures_all_CDH1_2_bottom(shape_measures_all_CDH1_2_bottom(:,1)==1,3),shape_measures_all_CDH1_6_bottom(shape_measures_all_CDH1_6_bottom(:,1)==0,3),shape_measures_all_CDH1_6_bottom(shape_measures_all_CDH1_6_bottom(:,1)==1,3),shape_measures_all_NT2_bottom(shape_measures_all_NT2_bottom(:,1)==0,3),shape_measures_all_NT2_bottom(shape_measures_all_NT2_bottom(:,1)==1,3),shape_measures_all_NT6_bottom(shape_measures_all_NT6_bottom(:,1)==0,3),shape_measures_all_NT6_bottom(shape_measures_all_NT6_bottom(:,1)==1,3),'CDH1_2mg/ml_outside','CDH1_2mg/ml_inside','CDH1_6mg/ml_outside','CDH1_6mg/ml_inside','NT_2mg/ml_outside','NT_2mg/ml_inside','NT_6mg/ml_outside','NT_6mg/ml_inside')
 
	my_boxplot(shape_measures_all_CDH1_2_bottom(shape_measures_all_CDH1_2_bottom(:,1)==0,2),shape_measures_all_CDH1_2_bottom(shape_measures_all_CDH1_2_bottom(:,1)==1,2),shape_measures_all_CDH1_6_bottom(shape_measures_all_CDH1_6_bottom(:,1)==0,2),shape_measures_all_CDH1_6_bottom(shape_measures_all_CDH1_6_bottom(:,1)==1,2),shape_measures_all_NT2_bottom(shape_measures_all_NT2_bottom(:,1)==0,2),shape_measures_all_NT2_bottom(shape_measures_all_NT2_bottom(:,1)==1,2),shape_measures_all_NT6_bottom(shape_measures_all_NT6_bottom(:,1)==0,2),shape_measures_all_NT6_bottom(shape_measures_all_NT6_bottom(:,1)==1,2),'CDH1_2mg/ml_outside','CDH1_2mg/ml_inside','CDH1_6mg/ml_outside','CDH1_6mg/ml_inside','NT_2mg/ml_outside','NT_2mg/ml_inside','NT_6mg/ml_outside','NT_6mg/ml_inside')
 
 	my_boxplot(shape_measures_all_CDH1_2_bottom(shape_measures_all_CDH1_2_bottom(:,1)==0,6),shape_measures_all_CDH1_2_bottom(shape_measures_all_CDH1_2_bottom(:,1)==1,6),shape_measures_all_CDH1_6_bottom(shape_measures_all_CDH1_6_bottom(:,1)==0,6),shape_measures_all_CDH1_6_bottom(shape_measures_all_CDH1_6_bottom(:,1)==1,6),shape_measures_all_NT2_bottom(shape_measures_all_NT2_bottom(:,1)==0,6),shape_measures_all_NT2_bottom(shape_measures_all_NT2_bottom(:,1)==1,6),shape_measures_all_NT6_bottom(shape_measures_all_NT6_bottom(:,1)==0,6),shape_measures_all_NT6_bottom(shape_measures_all_NT6_bottom(:,1)==1,6),'CDH1_2mg/ml_outside','CDH1_2mg/ml_inside','CDH1_6mg/ml_outside','CDH1_6mg/ml_inside','NT_2mg/ml_outside','NT_2mg/ml_inside','NT_6mg/ml_outside','NT_6mg/ml_inside')
    
    my_boxplot(shape_measures_all_CDH1_2_bottom(shape_measures_all_CDH1_2_bottom(:,1)==0,5),shape_measures_all_CDH1_2_bottom(shape_measures_all_CDH1_2_bottom(:,1)==1,5),shape_measures_all_CDH1_6_bottom(shape_measures_all_CDH1_6_bottom(:,1)==0,5),shape_measures_all_CDH1_6_bottom(shape_measures_all_CDH1_6_bottom(:,1)==1,5),shape_measures_all_NT2_bottom(shape_measures_all_NT2_bottom(:,1)==0,5),shape_measures_all_NT2_bottom(shape_measures_all_NT2_bottom(:,1)==1,5),shape_measures_all_NT6_bottom(shape_measures_all_NT6_bottom(:,1)==0,5),shape_measures_all_NT6_bottom(shape_measures_all_NT6_bottom(:,1)==1,5),'CDH1_2mg/ml_outside','CDH1_2mg/ml_inside','CDH1_6mg/ml_outside','CDH1_6mg/ml_inside','NT_2mg/ml_outside','NT_2mg/ml_inside','NT_6mg/ml_outside','NT_6mg/ml_inside')
  
    
	my_boxplot(shape_measures_all_CDH1_2_bottom(:,3),shape_measures_all_CDH1_2_up(:,3),shape_measures_all_CDH1_6_bottom(:,3),shape_measures_all_CDH1_6_up(:,3),shape_measures_all_NT2_bottom(:,3),shape_measures_all_NT2_up(:,3),shape_measures_all_NT6_bottom(:,3),shape_measures_all_NT6_up(:,3),'CDH1_2mg/ml_bottom','CDH1_2mg/ml_up','CDH1_6mg/ml_bottom','CDH1_6mg/ml_up','NT_2mg/ml_bottom','NT_2mg/ml_up','NT_6mg/ml_bottom','NT_6mg/ml_up')



nanmedian(shape_measures_all_CDH1_2_bottom(shape_measures_all_CDH1_2_bottom(:,1)==0,6))
nanmedian(shape_measures_all_CDH1_2_bottom(shape_measures_all_CDH1_2_bottom(:,1)==1,6))
nanmedian(shape_measures_all_CDH1_6_bottom(shape_measures_all_CDH1_6_bottom(:,1)==0,6))
nanmedian(shape_measures_all_CDH1_6_bottom(shape_measures_all_CDH1_6_bottom(:,1)==1,6))

nanmedian(shape_measures_all_NT2_bottom(shape_measures_all_NT2_bottom(:,1)==0,6))
nanmedian(shape_measures_all_NT2_bottom(shape_measures_all_NT2_bottom(:,1)==1,6))
nanmedian(shape_measures_all_NT6_bottom(shape_measures_all_NT6_bottom(:,1)==0,6))
nanmedian(shape_measures_all_NT6_bottom(shape_measures_all_NT6_bottom(:,1)==1,6))
end






if false
    shape_measures_all_CDH1_2_bottom(:,2)=shape_measures_all_CDH1_2_bottom(:,2)*0.6210004^2;    
    shape_measures_all_CDH1_6_bottom(:,2)=shape_measures_all_CDH1_6_bottom(:,2)*0.6210004^2;    
    shape_measures_all_CDH1_2_up(:,2)=shape_measures_all_CDH1_2_up(:,2)*0.6210004^2;    
    shape_measures_all_CDH1_6_up(:,2)=shape_measures_all_CDH1_6_up(:,2)*0.6210004^2;
    shape_measures_all_NT2_bottom(:,2)=shape_measures_all_NT2_bottom(:,2)*0.6210004^2;    
    shape_measures_all_NT6_bottom(:,2)=shape_measures_all_NT6_bottom(:,2)*0.6210004^2;    
    shape_measures_all_NT2_up(:,2)=shape_measures_all_NT2_up(:,2)*0.6210004^2;    
    shape_measures_all_NT6_up(:,2)=shape_measures_all_NT6_up(:,2)*0.6210004^2;
    
end