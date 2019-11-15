function [spa_vel_corr,spa_or_corr]=spatial_corr_tracks3D(Tracks,FrameTracks,time_windows,selection,mueh_per_pixel)

xvelocity_all=[];
yvelocity_all=[];
zvelocity_all=[];
xor_all=[];
yor_all=[];
zor_all=[];
    
for f=time_windows(1):time_windows(2)
    xvelocity_all=[xvelocity_all, arrayfun(@(NT) (Tracks{NT}(Tracks{NT}(:,1)==f-1,29)),selection{f})];
	yvelocity_all=[yvelocity_all, arrayfun(@(NT) (Tracks{NT}(Tracks{NT}(:,1)==f-1,30)),selection{f})];  
	zvelocity_all=[zvelocity_all, arrayfun(@(NT) (Tracks{NT}(Tracks{NT}(:,1)==f-1,31)),selection{f})]; %%% AHHH
    xor_all=[xor_all, arrayfun(@(NT) (Tracks{NT}(Tracks{NT}(:,1)==f-1,29)/((Tracks{NT}(Tracks{NT}(:,1)==f-1,29))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,30))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,31))^2)^.5),selection{f})];
	yor_all=[yor_all, arrayfun(@(NT) (Tracks{NT}(Tracks{NT}(:,1)==f-1,30)/((Tracks{NT}(Tracks{NT}(:,1)==f-1,29))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,30))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,31))^2)^.5),selection{f})];   
	zor_all=[zor_all, arrayfun(@(NT) (Tracks{NT}(Tracks{NT}(:,1)==f-1,31)/((Tracks{NT}(Tracks{NT}(:,1)==f-1,29))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,30))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,31))^2)^.5),selection{f})];   
end
    xvelocity_all(isnan(xvelocity_all))=[];
	yvelocity_all(isnan(yvelocity_all))=[]; 
    zvelocity_all(isnan(zvelocity_all))=[]; 
    xor_all(isnan(xor_all))=[];
	yor_all(isnan(yor_all))=[]; 
    zor_all(isnan(zor_all))=[]; 
    
    xvelocity_m_all=nanmean(xvelocity_all);
 	yvelocity_m_all=nanmean(yvelocity_all);    
    zvelocity_m_all=nanmean(zvelocity_all); 
    xor_m_all=nanmean(xor_all);
 	yor_m_all=nanmean(yor_all);    
    zor_m_all=nanmean(zor_all);
   
%    abs_of_mean_vel_vec_all=(nanmean(xvelocity_all)^2+nanmean(yvelocity_all)^2)^.5;
	var_of_vel_vec_all=nanvar(xvelocity_all)+nanvar(yvelocity_all);%+nanvar(zvelocity_all);
    var_of_or_all=nanvar(xor_all)+nanvar(yor_all);%+nanvar(zor_all);
  
    
    spa_vel_corr=[];
    spa_or_corr=[];

    
for f=time_windows(1):time_windows(2)
    FrameTrack_with_vel=FrameTracks{f}(ismember(FrameTracks{f},FrameTracks{f-1}));    
    selection_with_vel=selection{f}(ismember(selection{f},FrameTracks{f-1}));
    FrameTrack_with_vel=FrameTrack_with_vel(~logical(sum(isnan(cell2mat(arrayfun(@(tiq) Tracks{tiq}(Tracks{tiq}(:,1)==f-1,29:30)',FrameTrack_with_vel,'UniformOutput',false))))));
 	distance_vec_all=ipdm(mueh_per_pixel*cell2mat(arrayfun(@(tiq) Tracks{tiq}(Tracks{tiq}(:,1)==f-1,2:4)',selection_with_vel,'UniformOutput',false))',mueh_per_pixel*cell2mat(arrayfun(@(tiq) Tracks{tiq}(Tracks{tiq}(:,1)==f-1,2:4)',FrameTrack_with_vel,'UniformOutput',false))','Subset','Maximum', 'Limit',100,'Result','Structure');
    spa_vel_corr=[spa_vel_corr , [1/(var_of_vel_vec_all)*((arrayfun(@(NT) Tracks{NT}(Tracks{NT}(:,1)==f-1,29),selection_with_vel(distance_vec_all.rowindex))-xvelocity_m_all).*(arrayfun(@(NT) Tracks{NT}(Tracks{NT}(:,1)==f-1,29),FrameTrack_with_vel(distance_vec_all.columnindex))-xvelocity_m_all)+(arrayfun(@(NT) Tracks{NT}(Tracks{NT}(:,1)==f-1,30),selection_with_vel(distance_vec_all.rowindex))-yvelocity_m_all).*(arrayfun(@(NT) Tracks{NT}(Tracks{NT}(:,1)==f-1,30),FrameTrack_with_vel(distance_vec_all.columnindex))-yvelocity_m_all));distance_vec_all.distance']];
    spa_or_corr=[spa_or_corr , [1/(var_of_or_all)*((arrayfun(@(NT) Tracks{NT}(Tracks{NT}(:,1)==f-1,29)/(((Tracks{NT}(Tracks{NT}(:,1)==f-1,29))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,30))^2)^.5),selection_with_vel(distance_vec_all.rowindex))-xor_m_all).*(arrayfun(@(NT) Tracks{NT}(Tracks{NT}(:,1)==f-1,29)/(((Tracks{NT}(Tracks{NT}(:,1)==f-1,29))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,30))^2)^.5),FrameTrack_with_vel(distance_vec_all.columnindex))-xor_m_all)+(arrayfun(@(NT) Tracks{NT}(Tracks{NT}(:,1)==f-1,30)/(((Tracks{NT}(Tracks{NT}(:,1)==f-1,29))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,30))^2)^.5),selection_with_vel(distance_vec_all.rowindex))-yor_m_all).*(arrayfun(@(NT) Tracks{NT}(Tracks{NT}(:,1)==f-1,30)/(((Tracks{NT}(Tracks{NT}(:,1)==f-1,29))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,30))^2)^.5),FrameTrack_with_vel(distance_vec_all.columnindex))-yor_m_all));distance_vec_all.distance';mean([arrayfun(@(NT) (((Tracks{NT}(Tracks{NT}(:,1)==f-1,29))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,30))^2)^.5),selection_with_vel(distance_vec_all.rowindex));arrayfun(@(NT) (((Tracks{NT}(Tracks{NT}(:,1)==f-1,29))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,30))^2)^.5),FrameTrack_with_vel(distance_vec_all.columnindex))],1)]];

    %3D corr below - z possition to random
%     spa_vel_corr=[spa_vel_corr , [1/(var_of_vel_vec_all)*((arrayfun(@(NT) Tracks{NT}(Tracks{NT}(:,1)==f-1,29),selection_with_vel(distance_vec_all.rowindex))-xvelocity_m_all).*(arrayfun(@(NT) Tracks{NT}(Tracks{NT}(:,1)==f-1,29),FrameTrack_with_vel(distance_vec_all.columnindex))-xvelocity_m_all)+(arrayfun(@(NT) Tracks{NT}(Tracks{NT}(:,1)==f-1,30),selection_with_vel(distance_vec_all.rowindex))-yvelocity_m_all).*(arrayfun(@(NT) Tracks{NT}(Tracks{NT}(:,1)==f-1,30),FrameTrack_with_vel(distance_vec_all.columnindex))-yvelocity_m_all)+(arrayfun(@(NT) Tracks{NT}(Tracks{NT}(:,1)==f-1,31),selection_with_vel(distance_vec_all.rowindex))-zvelocity_m_all).*(arrayfun(@(NT) Tracks{NT}(Tracks{NT}(:,1)==f-1,31),FrameTrack_with_vel(distance_vec_all.columnindex))-zvelocity_m_all));distance_vec_all.distance']];
%     spa_or_corr=[spa_or_corr , [1/(var_of_or_all)*((arrayfun(@(NT) Tracks{NT}(Tracks{NT}(:,1)==f-1,29)/(((Tracks{NT}(Tracks{NT}(:,1)==f-1,29))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,30))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,31))^2)^.5),selection_with_vel(distance_vec_all.rowindex))-xor_m_all).*(arrayfun(@(NT) Tracks{NT}(Tracks{NT}(:,1)==f-1,29)/(((Tracks{NT}(Tracks{NT}(:,1)==f-1,29))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,30))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,31))^2)^.5),FrameTrack_with_vel(distance_vec_all.columnindex))-xor_m_all)+(arrayfun(@(NT) Tracks{NT}(Tracks{NT}(:,1)==f-1,30)/(((Tracks{NT}(Tracks{NT}(:,1)==f-1,29))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,30))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,31))^2)^.5),selection_with_vel(distance_vec_all.rowindex))-yor_m_all).*(arrayfun(@(NT) Tracks{NT}(Tracks{NT}(:,1)==f-1,30)/(((Tracks{NT}(Tracks{NT}(:,1)==f-1,29))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,30))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,31))^2)^.5),FrameTrack_with_vel(distance_vec_all.columnindex))-yor_m_all)+(arrayfun(@(NT) Tracks{NT}(Tracks{NT}(:,1)==f-1,31)/(((Tracks{NT}(Tracks{NT}(:,1)==f-1,29))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,30))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,31))^2)^.5),selection_with_vel(distance_vec_all.rowindex))-zor_m_all).*(arrayfun(@(NT) Tracks{NT}(Tracks{NT}(:,1)==f-1,31)/(((Tracks{NT}(Tracks{NT}(:,1)==f-1,29))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,30))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,31))^2)^.5),FrameTrack_with_vel(distance_vec_all.columnindex))-zor_m_all));distance_vec_all.distance';mean([arrayfun(@(NT) (((Tracks{NT}(Tracks{NT}(:,1)==f-1,29))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,30))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,31))^2)^.5),selection_with_vel(distance_vec_all.rowindex));arrayfun(@(NT) (((Tracks{NT}(Tracks{NT}(:,1)==f-1,29))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,30))^2+(Tracks{NT}(Tracks{NT}(:,1)==f-1,31))^2)^.5),FrameTrack_with_vel(distance_vec_all.columnindex))],1)]];

    
%    mobility=arrayfun(@(tiq) Tracks{tiq}(Tracks{tiq}(:,1)==f-1,32),FrameTracks{f});
%    mob_corr_all=[mob_corr_all,[1/(var_of_mel_vec_abs_all)*((arrayfun(@(NT) Tracks{NT}(Tracks{NT}(:,1)==f-1,32),FrameTrack_with_vel(distance_vec_all.rowindex))-velocity_m_abs_all).*(arrayfun(@(NT) Tracks{NT}(Tracks{NT}(:,1)==f-1,32),FrameTrack_with_vel(distance_vec_all.columnindex))-velocity_m_abs_all));distance_vec_all.distance']];
%    or_corr_all=[or_corr_all,1/mean_orr_sqrt*(arrayfun(@(NT) Tracks{NT}(Tracks{NT}(:,1)==f-1,33),FrameTrack_with_vel(distance_vec_all.rowindex)).*arrayfun(@(NT) conj(Tracks{NT}(Tracks{NT}(:,1)==f-1,33)),FrameTrack_with_vel(distance_vec_all.columnindex)))];
end

%mob_corr_all(:,isnan(mob_corr_all(1,:)))=[];    
%corr_all_indi(:,(isnan(corr_all_indi(1,:))))=[];
%or_corr_all(:,(isnan(or_corr_all(1,:))))=[];
%corr_all{trackys}=corr_all_indi;

% 
% [xDataOut, yDataOut, yCI] = SteffenSmoothCI(spa_or_corr(2,:),spa_or_corr(1,:),3,50);
% plot(xDataOut(yCI<50), yDataOut(yCI<50),'LineWidth',3)

end