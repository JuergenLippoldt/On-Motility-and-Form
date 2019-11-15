% extract some track properties 
% the calculations are trivial - the problem is mainly not to mix up the data structure 

pathy='E:\Friedl_19_02_13\tracks to analyse';

cd(pathy)

track_files=dir('*.xml');

for trackys=1:length(track_files)

trackname=track_files(trackys).name;

%% load and catagorize tracks

[Tracks, ~] = importTrackMateTracks([pathy, '\' ,trackname]);

mueh_per_pixel=1.05;  %1/0.95234

Tracks(cellfun('length',Tracks) <11)=[];

% the third-dimension 0 seems to bother the code
for j=1:length(Tracks)
    if size(Tracks{j},2)>3
    Tracks{j}(:,4)=[]; %I am to dumb for cellfun :P
    end
end

% the code assumes that the tracks have a position at each time step between start and finish 
% if there is a gap it is estimate as the mean between the points between the gap
% does not happen often and should not change anything significant
Tracks=fill_gaps(Tracks);

for j=1:length(Tracks)
   
	Tracks{j}=round(Tracks{j});

end

for f=length(Tracks):-1:1  %sort out tracks of faulty pixel
    if sum(sum(abs(Tracks{f}(:,2:3)-Tracks{f}(1,2:3))))<10
        Tracks(f)=[];
    end    
end

T = Tracks;
% This tells for each frame the tracks of which it is made up and makes a
% structure with it for later use
[FrameTracks, FrameTrackCoordinates] = FindAllTracksInFrames(T);
clear T
FrameTrackCoordinates = EliminateDuplicateDataPoints( FrameTracks, FrameTrackCoordinates, Tracks );
FrameTrackCoordinates = EliminateDuplicateDataPoints( FrameTracks, FrameTrackCoordinates, Tracks );

cellnumber{trackys}=arrayfun(@(blub) length(FrameTracks{blub}),1:length( FrameTracks),'UniformOutput',false);

% these strangly high numbers are used, because this script is copied and
% pasted from a script that is looking at much more involved stuff
for j=1:length(Tracks)
    Tracks{j}(:,7)=NaN;
    Tracks{j}(:,11:12)=NaN;
    Tracks{j}(:,17:28)=NaN;
    Tracks{j}(:,29)=0;
    Tracks{j}(:,31)=NaN;
	Tracks{j}(:,33)=NaN;
    Tracks{j}(:,34)=0;
    Tracks{j}(:,39)=2;
end

for j=1:length(Tracks)
    Tracks{j}(2:end,29)=Tracks{j}(2:end,2)-Tracks{j}(1:end-1,2);
	Tracks{j}(2:end,30)=Tracks{j}(2:end,3)-Tracks{j}(1:end-1,3);
    Tracks{j}(1,29:30)=NaN;
    Tracks{j}(4:end-3,32)=mueh_per_pixel*((Tracks{j}(7:end,2)-Tracks{j}(1:end-6,2)).^2+(Tracks{j}(7:end,3)-Tracks{j}(1:end-6,3)).^2+(Tracks{j}(7:end,4)-Tracks{j}(1:end-6,4)).^2).^0.5;
    Tracks{j}(1:min(end,3),32)=NaN; Tracks{j}(max(end-3,1):end,32)=NaN; 
end    

for f=2:length(FrameTracks)
    FrameTrack_with_vel=FrameTracks{f}(ismember(FrameTracks{f},FrameTracks{f-1}));
    FrameTrack_with_vel=FrameTrack_with_vel(~logical(sum(isnan(cell2mat(arrayfun(@(tiq) Tracks{tiq}(Tracks{tiq}(:,1)==f-1,29:30)',FrameTrack_with_vel,'UniformOutput',false))))));
 	distance_vec_all=ipdm(mueh_per_pixel*cell2mat(arrayfun(@(tiq) Tracks{tiq}(Tracks{tiq}(:,1)==f-1,2:3)',FrameTrack_with_vel,'UniformOutput',false))','Subset','Maximum', 'Limit',25,'Result','Structure');
    
    for j=1:length(Tracks)
        if sum(distance_vec_all.columnindex == j)>2
            neighbourtracks=FrameTrack_with_vel(distance_vec_all.rowindex(distance_vec_all.columnindex == j & distance_vec_all.rowindex ~=j));
            Tracks{j}(Tracks{j}(:,1)==f-1,35)=Tracks{j}(Tracks{j}(:,1)==f-1,29)-nanmean(arrayfun(@(tiq) Tracks{tiq}(Tracks{tiq}(:,1)==f-1,29),neighbourtracks));
            Tracks{j}(Tracks{j}(:,1)==f-1,36)=Tracks{j}(Tracks{j}(:,1)==f-1,30)-nanmean(arrayfun(@(tiq) Tracks{tiq}(Tracks{tiq}(:,1)==f-1,30),neighbourtracks));
        else
            Tracks{j}(Tracks{j}(:,1)==f-1,35)=NaN;
            Tracks{j}(Tracks{j}(:,1)==f-1,36)=NaN;
        end

    end        
end

for j=1:length(Tracks)
	Tracks{j}(1,37:38)=Tracks{j}(1,2:3);
	for t=2:size(Tracks{j},1)
        Tracks{j}(t,37:38)=Tracks{j}(t-1,37:38)+Tracks{j}(t,35:36);
	end
end


%% categorize regions for tracks
% distance to center < 100 pixel => core 
% 10% of furthest outside cells => edge 
for t=1:length(FrameTracks)
    frame_dist=[];
    	center_x(t)=mean( arrayfun(@(tiq) Tracks{tiq}(Tracks{tiq}(:,1)==t-1,2),FrameTracks{t}));
        center_y(t)=mean( arrayfun(@(tiq) Tracks{tiq}(Tracks{tiq}(:,1)==t-1,3),FrameTracks{t}));
    for f=1:length(FrameTracks{t})
        frame_dist(f)=((Tracks{FrameTracks{t}(f)}(Tracks{FrameTracks{t}(f)}(:,1)==t-1,2)-center_x(t))^2+(Tracks{FrameTracks{t}(f)}(Tracks{FrameTracks{t}(f)}(:,1)==t-1,3)-center_y(t))^2)^.5;
        if frame_dist(f)<100
            Tracks{FrameTracks{t}(f)}(Tracks{FrameTracks{t}(f)}(:,1)==t-1,39)=1;
        end
    end
    [~,Idx] = sort(frame_dist);
    for i=round(.9*length(FrameTracks{t})):length(FrameTracks{t})
        Tracks{FrameTracks{t}(Idx(i))}(Tracks{FrameTracks{t}(Idx(i))}(:,1)==t-1,39)=3;
    end
end

for i=1:length(FrameTracks)
    FrameTracks_1{i}=FrameTracks{i}(arrayfun(@(tiq) Tracks{tiq}(Tracks{tiq}(:,1)==i-1,39)==1,FrameTracks{i}));
    FrameTracks_2{i}=FrameTracks{i}(arrayfun(@(tiq) Tracks{tiq}(Tracks{tiq}(:,1)==i-1,39)==2,FrameTracks{i}));
    FrameTracks_3{i}=FrameTracks{i}(arrayfun(@(tiq) Tracks{tiq}(Tracks{tiq}(:,1)==i-1,39)==3,FrameTracks{i}));
end

%% calculate MSD and spatial velocity correlation
% delta t - distance pairs
time_windows(1)=5;
time_windows(2)=length(FrameTracks)-5;
time_frame_of_interest=time_windows(2)-time_windows(1)-20;
dt_dist_sq=cell(min(381,time_frame_of_interest),2);

for i=1:length(Tracks)
    trq=Tracks{i}(Tracks{i}(:,1)>time_windows(1) & Tracks{i}(:,1)<time_windows(2),[2:3,39]);
    if size(trq,1)>2
        for tau=0:min(length(trq)-1,min(380,time_frame_of_interest-1))
            dt_dist_sq{tau+1}(end+1:end+length(trq)-tau,:)=[(trq(1:end-tau,1)-trq(1+tau:end,1)).^2+(trq(1:end-tau,2)-trq(1+tau:end,2)).^2, trq(1:end-tau,3) ];
        end
    end
    
end

for t=1:length(dt_dist_sq)
    MSD_indi(t,1)=mean(dt_dist_sq{t}(:,1));
    MSD_indi(t,2)=mean(dt_dist_sq{t}(dt_dist_sq{t}(:,2)==1,1));    
    MSD_indi(t,3)=mean(dt_dist_sq{t}(dt_dist_sq{t}(:,2)==2,1));    
    MSD_indi(t,4)=mean(dt_dist_sq{t}(dt_dist_sq{t}(:,2)==3,1));
%     MSD_normed_indi(t)=nanmean(dt_dist_normed_sq{t});
    non_gauss_indi(t,1)=mean(dt_dist_sq{t}(:,1).^2)/(3*(mean(dt_dist_sq{t}(:,1)).^2))-1;
end
 MSD{trackys}=MSD_indi;
%  MSD_normed{trackys}=MSD_normed_indi;
 non_gauss{trackys}=non_gauss_indi;

[spa_vel_corr_1,spa_or_corr_1]=spatial_corr_tracks(Tracks,FrameTracks,time_windows,FrameTracks_1,mueh_per_pixel);
[spa_vel_corr_2,spa_or_corr_2]=spatial_corr_tracks(Tracks,FrameTracks,time_windows,FrameTracks_2,mueh_per_pixel);
[spa_vel_corr_3,spa_or_corr_3]=spatial_corr_tracks(Tracks,FrameTracks,time_windows,FrameTracks_3,mueh_per_pixel);
[spa_vel_corr_all,spa_or_corr_all]=spatial_corr_tracks(Tracks,FrameTracks,time_windows,FrameTracks,mueh_per_pixel);

spatial_velocity_corr{trackys,1}=spa_vel_corr_all;
spatial_velocity_corr{trackys,2}=spa_or_corr_all;
spatial_velocity_corr{trackys,3}=spa_vel_corr_1;
spatial_velocity_corr{trackys,4}=spa_or_corr_1;
spatial_velocity_corr{trackys,5}=spa_vel_corr_2;
spatial_velocity_corr{trackys,6}=spa_or_corr_2;
spatial_velocity_corr{trackys,7}=spa_vel_corr_3;
spatial_velocity_corr{trackys,8}=spa_or_corr_3;

%% temporal correlation

xvelocity_all=[];
yvelocity_all=[];
xor_all=[];
yor_all=[];
    
for f=time_windows(1):time_windows(2)
    xvelocity_all=[xvelocity_all, arrayfun(@(NT) (Tracks{NT}(Tracks{NT}(:,1)==f-1,29)),FrameTracks{f})];
	yvelocity_all=[yvelocity_all, arrayfun(@(NT) (Tracks{NT}(Tracks{NT}(:,1)==f-1,30)),FrameTracks{f})];  
    xor_all=[xor_all, arrayfun(@(NT) (Tracks{NT}(Tracks{NT}(:,1)==f-1,29)/((Tracks{NT}(Tracks{NT}(:,1)==f-1,29))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,30))^2)^.5),FrameTracks{f})];
	yor_all=[yor_all, arrayfun(@(NT) (Tracks{NT}(Tracks{NT}(:,1)==f-1,30)/((Tracks{NT}(Tracks{NT}(:,1)==f-1,29))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,30))^2)^.5),FrameTracks{f})];   
end
    xvelocity_all(isnan(xvelocity_all))=[];
	yvelocity_all(isnan(yvelocity_all))=[]; 
    xor_all(isnan(xor_all))=[];
	yor_all(isnan(yor_all))=[]; 
    
    xvelocity_m_all=nanmean(xvelocity_all);
 	yvelocity_m_all=nanmean(yvelocity_all);    
    xor_m_all=nanmean(xor_all);
 	yor_m_all=nanmean(yor_all);  
    
	var_of_vel_vec_all=(nanvar(xvelocity_all)+nanvar(yvelocity_all));
    var_of_or_all=(nanvar(xor_all)+nanvar(yor_all));
   

Autocorr_i=cell(min(61,time_windows(2)-time_windows(1)),2);
Autocorr_or_i=cell(min(61,time_windows(2)-time_windows(1)),2);

for i=1:length(Tracks)
	trq=Tracks{i}(Tracks{i}(:,1)>time_windows(1) & Tracks{i}(:,1)<time_windows(2),[2:3,39]);
    trq_vel=Tracks{i}(Tracks{i}(:,1)>time_windows(1) & Tracks{i}(:,1)<time_windows(2),29:30);
    trq_or=Tracks{i}(Tracks{i}(:,1)>time_windows(1) & Tracks{i}(:,1)<time_windows(2),29:30)./(Tracks{i}(Tracks{i}(:,1)>time_windows(1) & Tracks{i}(:,1)<time_windows(2),29).^2+Tracks{i}(Tracks{i}(:,1)>time_windows(1) & Tracks{i}(:,1)<time_windows(2),30).^2).^.5;
    if size(trq,1)>3
        for tau=0:min(length(trq)-1,60)
        	Autocorr_i{tau+1}(end+1:end+length(trq)-tau,:)=[((trq_vel(1:end-tau,1)-xvelocity_m_all).*(trq_vel(1+tau:end,1)-xvelocity_m_all)+(trq_vel(1:end-tau,2)-yvelocity_m_all).*(trq_vel(1+tau:end,2)-yvelocity_m_all))/var_of_vel_vec_all,trq(1:end-tau,3)];
            Autocorr_or_i{tau+1}(end+1:end+length(trq)-tau,:)=[((trq_or(1:end-tau,1)-xor_m_all).*(trq_or(1+tau:end,1)-xor_m_all)+(trq_or(1:end-tau,2)-yor_m_all).*(trq_or(1+tau:end,2)-yor_m_all))/var_of_or_all,trq(1:end-tau,3)];
        end
    end
end
for i=1:min(61,time_windows(2)-time_windows(1))
    Autocorr{trackys}(i,1)=nanmean(Autocorr_i{i}(:,1));
    Autocorr{trackys}(i,3)=nanmean(Autocorr_i{i}(Autocorr_i{i}(:,2)==1,1));
    Autocorr{trackys}(i,5)=nanmean(Autocorr_i{i}(Autocorr_i{i}(:,2)==2,1));    
    Autocorr{trackys}(i,7)=nanmean(Autocorr_i{i}(Autocorr_i{i}(:,2)==3,1));
    Autocorr{trackys}(i,2)=nanmean(Autocorr_or_i{i}(:,1));
    Autocorr{trackys}(i,4)=nanmean(Autocorr_or_i{i}(Autocorr_or_i{i}(:,2)==1,1));
    Autocorr{trackys}(i,6)=nanmean(Autocorr_or_i{i}(Autocorr_or_i{i}(:,2)==2,1));    
    Autocorr{trackys}(i,8)=nanmean(Autocorr_or_i{i}(Autocorr_or_i{i}(:,2)==3,1));
end



clearvars -except Autocorr MSD corr_all non_gauss trackys track_files pathy  MSD_normed spatial_velocity_corr

end

%% plot scripts to copy + paste

if false
    
for trackys=1:length(track_files)

trackname=track_files(trackys).name;

if contains( trackname,'NT')
    if contains( trackname,'2_mg_ml') || contains( trackname,'2 mg')
        color=[65,171,93]/255;
    elseif contains( trackname,'6_mg_ml') || contains( trackname,'6 mg')
        color=[0,109,44]/255;
    end
    
elseif  contains( trackname,'CDH')
	if contains( trackname,'2_mg_ml') || contains( trackname,'2 mg')
        color=[33,113,181]/255;
	elseif contains( trackname,'6_mg_ml') || contains( trackname,'6 mg')
        color=[8,48,107]/255;
    end
end


% hold on
% [xDataOut, yDataOut, yCI] = SteffenSmoothCI(spatial_velocity_corr{trackys,2}(2,:),spatial_velocity_corr{trackys,2}(1,:),3,50);
% p{trackys}=plot(xDataOut(xDataOut>15), yDataOut(xDataOut>15),'LineWidth',3,'Color',color);
% 
 hold on
% ciplot_steffen(xDataOut(xDataOut>7), yDataOut(xDataOut>7)-yCI(xDataOut>7),yDataOut(xDataOut>7)+yCI(xDataOut>7),'red')
 p{trackys}=plot((10:length(cellnumber{trackys})-10)/4,cell2mat(cellnumber{trackys}(10:end-10)),'LineWidth',3,'Color',color);

end
l = findobj(gca,'Type','Line');
legend([p{1} p{6} p{11} p{16}],{'CDH1 2 mg/ml','CDH1 6 mg/ml','NT 2 mg/ml','NT 6 mg/ml'})



figure;
plot((0:70)/4,MSD_shNT(1:71),'LineWidth',3)
hold on
plot((0:70)/4,MSD_CDH1(1:71),'LineWidth',3)

figure;
plot((1:60)/4,Autocorr_shNT(2:61),'LineWidth',3)
hold on
plot((1:60)/4,Autocorr_CDH1(2:61),'LineWidth',3)

figure;
plot((0:70)/4,non_gauss_shNT(1:71),'LineWidth',3)
hold on
plot((0:70)/4,non_gauss_CDH1(1:71),'LineWidth',3)

hold on
[xDataOut, yDataOut, yCI] = SteffenSmoothCI(corr_all_shNT(2,:),corr_all_shNT(1,:),3,50);
plot(xDataOut(xDataOut>7), yDataOut(xDataOut>7),'LineWidth',3)
hold on
ciplot_steffen(xDataOut(xDataOut>7), yDataOut(xDataOut>7)-yCI(xDataOut>7),yDataOut(xDataOut>7)+yCI(xDataOut>7),'blue')

hold on
[xDataOut, yDataOut, yCI] = SteffenSmoothCI(corr_all_CDH1(2,:),corr_all_CDH1(1,:),3,50);
plot(xDataOut(xDataOut>7), yDataOut(xDataOut>7),'LineWidth',3)
hold on
ciplot_steffen(xDataOut(xDataOut>7), yDataOut(xDataOut>7)-yCI(xDataOut>7),yDataOut(xDataOut>7)+yCI(xDataOut>7),'red')

legend([p1 p3],{'shNT-2 mg movie 1','CDH1-2 mg movie 1'})

 NG2_all =[];
 NG2_core=[];
 NG6_core=[];
 CD2_core=[];
 CD6_core=[];
 NG2_expanse=[];
 NG6_expanse=[];
 CD2_expanse=[];
 CD6_expanse=[];
 NG2_edge=[];
 NG6_edge=[];
 CD2_edge=[];
 CD6_edge=[];
 
for trackys=1:length(track_files) 
trackname=track_files(trackys).name;
if contains( trackname,'NT')
    if contains( trackname,'2_mg_ml') || contains( trackname,'2 mg')
        color=[65,171,93]/255;
        NG2_all    =[NG2_all,     spatial_velocity_corr{trackys,2}];
        NG2_core   =[NG2_core,    spatial_velocity_corr{trackys,4}];
        NG2_expanse=[NG2_expanse, spatial_velocity_corr{trackys,6}];
        NG2_edge   =[NG2_edge,    spatial_velocity_corr{trackys,8}];     
    elseif contains( trackname,'6_mg_ml') || contains( trackname,'6 mg')
        color=[0,109,44]/255;
        NG6_core   =[NG6_core,    spatial_velocity_corr{trackys,4}];
        NG6_expanse=[NG6_expanse, spatial_velocity_corr{trackys,6}];
        NG6_edge   =[NG6_edge,    spatial_velocity_corr{trackys,8}];   
    end
    
elseif  contains( trackname,'CDH')
	if contains( trackname,'2_mg_ml') || contains( trackname,'2 mg')
        color=[33,113,181]/255;
        CD2_core   =[CD2_core,    spatial_velocity_corr{trackys,4}];
        CD2_expanse=[CD2_expanse, spatial_velocity_corr{trackys,6}];
        CD2_edge   =[CD2_edge,    spatial_velocity_corr{trackys,8}]; 
	elseif contains( trackname,'6_mg_ml') || contains( trackname,'6 mg')
        color=[8,48,107]/255;
        CD6_core   =[CD6_core,    spatial_velocity_corr{trackys,4}];
        CD6_expanse=[CD6_expanse, spatial_velocity_corr{trackys,6}];
        CD6_edge   =[CD6_edge,    spatial_velocity_corr{trackys,8}]; 
    end
end
end
if false 
    [xDataOut, yDataOut, zDataOut, zCI,  TWeight] = SteffenSmoothCI3D( NG2_all(2,NG2_all(2,:)>15),NG2_all(3,NG2_all(2,:)>15) ,NG2_all(1,NG2_all(2,:)>15), 2, 100, 100 );
    zDataOut2=zDataOut;
    zDataOut2(TWeight<5)=NaN;
%zDataOut2(zDataOut2==0)=0.001;
%zDataOut2=1./zDataOut2;
figure;
s=surf(yDataOut,xDataOut,zDataOut2);
caxis([0, 1])
s.EdgeColor = 'black';
end

if false %% statistic 
    [xDataOut, yDataOut_CD2_edge, weight_CD2_edge] = SteffenSmoothCI_total_weight(CD2_edge(2,:),CD2_edge(1,:),3,50);
    [xDataOut, yDataOut_CD6_edge, weight_CD6_edge] = SteffenSmoothCI_total_weight(CD6_edge(2,:),CD6_edge(1,:),3,50);
    [xDataOut, yDataOut_CD2_expanse, weight_CD2_expanse] = SteffenSmoothCI_total_weight(CD2_expanse(2,:),CD2_expanse(1,:),3,50);
    [xDataOut, yDataOut_CD6_expanse, weight_CD6_expanse] = SteffenSmoothCI_total_weight(CD6_expanse(2,:),CD6_expanse(1,:),3,50);
    [xDataOut, yDataOut_NG6_edge, weight_NG6_edge] = SteffenSmoothCI_total_weight(NG6_edge(2,:),NG6_edge(1,:),3,50);
    [xDataOut, yDataOut_NG6_expanse, weight_NG6_expanse] = SteffenSmoothCI_total_weight(NG6_expanse(2,:),NG6_expanse(1,:),3,50);
    [xDataOut, yDataOut_NG2_edge, weight_NG2_edge] = SteffenSmoothCI_total_weight(NG2_edge(2,:),NG2_edge(1,:),3,50);
    [xDataOut, yDataOut_NG2_expanse, weight_NG2_expanse] = SteffenSmoothCI_total_weight(NG2_expanse(2,:),NG2_expanse(1,:),3,50);
    
    for i=9:50
        p_CD2_edge_CD6_edge(i-8) = compare_correlation_coefficients(yDataOut_CD2_edge(i),yDataOut_CD6_edge(i),weight_CD2_edge(i)/10,weight_CD6_edge(i)/10);
        p_NT6_edge_CD6_edge(i-8) = compare_correlation_coefficients(yDataOut_NG6_edge(i),yDataOut_CD6_edge(i),weight_NG6_edge(i)/10,weight_CD6_edge(i)/10);
        p_CD2_mid_CD6_mid(i-8) = compare_correlation_coefficients(yDataOut_CD2_expanse(i),yDataOut_CD6_expanse(i),weight_CD2_expanse(i)/10,weight_CD6_expanse(i)/10);
        p_NT6_mid_CD6_mid(i-8) = compare_correlation_coefficients(yDataOut_NG6_expanse(i),yDataOut_CD6_expanse(i),weight_NG6_expanse(i)/10,weight_CD6_expanse(i)/10);
        p_CD2_edge_CD6_mid(i-8) = compare_correlation_coefficients(yDataOut_CD2_edge(i),yDataOut_CD6_expanse(i),weight_CD2_edge(i)/10,weight_CD6_expanse(i)/10);
        p_NT6_edge_CD6_mid(i-8) = compare_correlation_coefficients(yDataOut_NG6_edge(i),yDataOut_CD6_expanse(i),weight_NG6_edge(i)/10,weight_CD6_expanse(i)/10);
        p_NT6_edge_NT2_mid(i-8) = compare_correlation_coefficients(yDataOut_NG6_edge(i),yDataOut_NG2_expanse(i),weight_NG6_edge(i)/10,weight_NG2_expanse(i)/10);
        p_NT6_edge_NT2_edge(i-8) = compare_correlation_coefficients(yDataOut_NG6_edge(i),yDataOut_NG2_edge(i),weight_NG6_edge(i)/10,weight_NG2_edge(i)/10);
        p_NT6_mid_NT2_mid(i-8) = compare_correlation_coefficients(yDataOut_NG6_expanse(i),yDataOut_NG2_expanse(i),weight_NG6_expanse(i)/10,weight_NG2_expanse(i)/10);
        p_CD2_mid_CD6_mid(i-8) = compare_correlation_coefficients(yDataOut_CD6_expanse(i),yDataOut_CD2_expanse(i),weight_CD6_expanse(i)/10,weight_CD2_expanse(i)/10);
        
    end 
end

[xDataOut, yDataOut, yCI] = SteffenSmoothCI(CD6_core(2,:),CD6_core(1,:),3,50);
p{1}=plot(xDataOut(xDataOut>15), yDataOut(xDataOut>15),'LineWidth',3,'Color',[8,48,107]/255);
hold on
[xDataOut, yDataOut, yCI] = SteffenSmoothCI(CD2_core(2,:),CD2_core(1,:),3,50);
p{2}=plot(xDataOut(xDataOut>15), yDataOut(xDataOut>15),'LineWidth',3,'Color',[33,113,181]/255);
hold on
[xDataOut, yDataOut, yCI] = SteffenSmoothCI(NG6_core(2,:),NG6_core(1,:),3,50);
p{3}=plot(xDataOut(xDataOut>15), yDataOut(xDataOut>15),'LineWidth',3,'Color',[0,109,44]/255);
hold on
[xDataOut, yDataOut, yCI] = SteffenSmoothCI(NG2_core(2,:),NG2_core(1,:),3,50);
p{4}=plot(xDataOut(xDataOut>15), yDataOut(xDataOut>15),'LineWidth',3,'Color',[65,171,93]/255);
hold on
[xDataOut, yDataOut, yCI] = SteffenSmoothCI(CD6_expanse(2,:),CD6_expanse(1,:),3,50);
p{5}=plot(xDataOut(xDataOut>15), yDataOut(xDataOut>15),'LineWidth',3,'Color',[88,48,107]/255);
hold on
[xDataOut, yDataOut, yCI] = SteffenSmoothCI(CD2_expanse(2,:),CD2_expanse(1,:),3,50);
p{6}=plot(xDataOut(xDataOut>15), yDataOut(xDataOut>15),'LineWidth',3,'Color',[113,113,181]/255);
hold on
[xDataOut, yDataOut, yCI] = SteffenSmoothCI(NG6_expanse(2,:),NG6_expanse(1,:),3,50);
p{7}=plot(xDataOut(xDataOut>15), yDataOut(xDataOut>15),'LineWidth',3,'Color',[80,109,44]/255);
hold on
[xDataOut, yDataOut, yCI] = SteffenSmoothCI(NG2_expanse(2,:),NG2_expanse(1,:),3,50);
p{8}=plot(xDataOut(xDataOut>15), yDataOut(xDataOut>15),'LineWidth',3,'Color',[145,171,93]/255);
hold on
[xDataOut, yDataOut, yCI] = SteffenSmoothCI(CD6_edge(2,:),CD6_edge(1,:),3,50);
p{9}=plot(xDataOut(xDataOut>15), yDataOut(xDataOut>15),'LineWidth',3,'Color',[188,48,107]/255);
hold on
[xDataOut, yDataOut, yCI] = SteffenSmoothCI(CD2_edge(2,:),CD2_edge(1,:),3,50);
p{10}=plot(xDataOut(xDataOut>15), yDataOut(xDataOut>15),'LineWidth',3,'Color',[213,113,181]/255);
hold on
[xDataOut, yDataOut, yCI] = SteffenSmoothCI(NG6_edge(2,:),NG6_edge(1,:),3,50);
p{11}=plot(xDataOut(xDataOut>15), yDataOut(xDataOut>15),'LineWidth',3,'Color',[180,109,44]/255);
hold on
[xDataOut, yDataOut, yCI] = SteffenSmoothCI(NG2_edge(2,:),NG2_edge(1,:),3,50);
p{12}=plot(xDataOut(xDataOut>15), yDataOut(xDataOut>15),'LineWidth',3,'Color',[245,171,93]/255);

legend([p{2} p{6} p{10} p{1} p{5} p{9} p{4} p{8} p{12} p{3} p{7} p{11}],{'CDH1 2 mg/ml - core','CDH1 2 mg/ml - expanse','CDH1 2 mg/ml - edge','CDH1 6 mg/ml - core','CDH1 6 mg/ml - expanse','CDH1 6 mg/ml - edge','NT 2 mg/ml - core','NT 2 mg/ml - expanse','NT 2 mg/ml - edge','NT 6 mg/ml - core','NT 6 mg/ml - expanse','NT 6 mg/ml - edge'})



 NG2_core=[];
 NG6_core=[];
 CD2_core=[];
 CD6_core=[];
 NG2_expanse=[];
 NG6_expanse=[];
 CD2_expanse=[];
 CD6_expanse=[];
 NG2_edge=[];
 NG6_edge=[];
 CD2_edge=[];
 CD6_edge=[];
for trackys=1:length(track_files) 
trackname=track_files(trackys).name;
if contains( trackname,'NT')
    if contains( trackname,'2_mg_ml') || contains( trackname,'2 mg')
        color=[65,171,93]/255;
        NG2_core   =[NG2_core,    MSD{trackys}(1:70,2)];
        NG2_expanse=[NG2_expanse, MSD{trackys}(1:70,3)];
        NG2_edge   =[NG2_edge,    MSD{trackys}(1:70,4)];     
    elseif contains( trackname,'6_mg_ml') || contains( trackname,'6 mg')
        color=[0,109,44]/255;
        NG6_core   =[NG6_core,    MSD{trackys}(1:70,2)];
        NG6_expanse=[NG6_expanse, MSD{trackys}(1:70,3)];
        NG6_edge   =[NG6_edge,    MSD{trackys}(1:70,4)];   
    end
    
elseif  contains( trackname,'CDH')
	if contains( trackname,'2_mg_ml') || contains( trackname,'2 mg')
        color=[33,113,181]/255;
        CD2_core   =[CD2_core,    MSD{trackys}(1:70,2)];
        CD2_expanse=[CD2_expanse, MSD{trackys}(1:70,3)];
        CD2_edge   =[CD2_edge,    MSD{trackys}(1:70,4)]; 
	elseif contains( trackname,'6_mg_ml') || contains( trackname,'6 mg')
        color=[8,48,107]/255;
        CD6_core   =[CD6_core,    MSD{trackys}(1:70,2)];
        CD6_expanse=[CD6_expanse, MSD{trackys}(1:70,3)];
        CD6_edge   =[CD6_edge,    MSD{trackys}(1:70,4)]; 
    end
end
end
NG2_core=mean(NG2_core,2);
NG2_expanse=mean(NG2_expanse,2);
NG2_edge=mean(NG2_edge,2);
NG6_core=mean(NG6_core,2);
NG6_expanse=mean(NG6_expanse,2);
NG6_edge=mean(NG6_edge,2);
CD2_core=mean(CD2_core,2);
CD2_expanse=mean(CD2_expanse,2);
CD2_edge=mean(CD2_edge,2);
CD6_core=mean(CD6_core,2);
CD6_expanse=mean(CD6_expanse,2);
CD6_edge=mean(CD6_edge,2);

p{1}=plot(((0:69)/4), CD6_core,'LineWidth',3,'Color',[8,48,107]/255);
hold on
p{2}=plot(((0:69)/4), CD2_core,'LineWidth',3,'Color',[33,113,181]/255);
hold on
p{3}=plot(((0:69)/4), NG6_core,'LineWidth',3,'Color',[0,109,44]/255);
hold on
p{4}=plot(((0:69)/4), NG2_core,'LineWidth',3,'Color',[65,171,93]/255);
hold on
p{5}=plot(((0:69)/4), CD6_expanse,'LineWidth',3,'Color',[88,48,107]/255);
hold on
p{6}=plot(((0:69)/4), CD2_expanse,'LineWidth',3,'Color',[113,113,181]/255);
hold on
p{7}=plot(((0:69)/4), NG6_expanse,'LineWidth',3,'Color',[80,109,44]/255);
hold on
p{8}=plot(((0:69)/4), NG2_expanse,'LineWidth',3,'Color',[145,171,93]/255);
hold on
p{9}=plot(((0:69)/4), CD6_edge,'LineWidth',3,'Color',[188,48,107]/255);
hold on
p{10}=plot((0:69)/4, CD2_edge,'LineWidth',3,'Color',[213,113,181]/255);
hold on
p{11}=plot((0:69)/4, NG6_edge,'LineWidth',3,'Color',[180,109,44]/255);
hold on
p{12}=plot((0:69)/4, NG2_edge,'LineWidth',3,'Color',[245,171,93]/255);

hold on
plot (exp((-1:.5:2.5)),400*exp((-1:.5:2.5)),'color','red','linewidth',2)

hold on
plot (exp((-1:.5:2.5)),25*exp((-1:.5:2.5)).^1.20,'color','red','linewidth',2)
legend([p{2} p{6} p{10} p{1} p{5} p{9} p{4} p{8} p{12} p{3} p{7} p{11}],{'CDH1 2 mg/ml - core','CDH1 2 mg/ml - expanse','CDH1 2 mg/ml - edge','CDH1 6 mg/ml - core','CDH1 6 mg/ml - expanse','CDH1 6 mg/ml - edge','NT 2 mg/ml - core','NT 2 mg/ml - expanse','NT 2 mg/ml - edge','NT 6 mg/ml - core','NT 6 mg/ml - expanse','NT 6 mg/ml - edge'})
hold on
plot (exp((-1:.5:2.5)),200*exp((-1:.5:2.5)),'color','red','linewidth',3)

hold on
plot (exp((-1:.5:2.5)),32*exp((-1:.5:2.5)).^1.30,'color','red','linewidth',3)
end
