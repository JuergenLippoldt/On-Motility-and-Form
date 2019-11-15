% extract some track properties 
% the calculations are trivial - the problem is mainly not to mix up the data structure 

pathy='E:\Friedl_19_02_13\tracks to analyse';
cd(pathy)

track_files=dir('*.xml');

NT_vel_6=cell(1,3);
NT_vel_2=cell(1,3);
CDH1_vel_2=cell(1,3);
CDH1_vel_6=cell(1,3);
for trackys=1:length(track_files)
    trackname=track_files(trackys).name;
    disp(trackname)
    
    %% load and clean up
   [Tracks, ~] = importTrackMateTracks([pathy, '\' ,trackname]);

mueh_per_pixel=1.05;  %1/0.95234#

Tracks(cellfun('length',Tracks) <11)=[];

% the third-dimension 0 seems to bother the code
for j=1:length(Tracks)
    if size(Tracks{j},2)>3
    Tracks{j}(:,4)=[]; %I was to dumb for cellfun :P
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

%% calculate the parameter of interest and save it in the suitable data structure 
for j=1:length(Tracks)
     Tracks{j}(4:end-3,32)=mueh_per_pixel*((Tracks{j}(6:end-1,2)-Tracks{j}(2:end-5,2)).^2+(Tracks{j}(6:end-1,3)-Tracks{j}(2:end-5,3)).^2+(Tracks{j}(6:end-1,4)-Tracks{j}(2:end-5,4)).^2).^0.5;
end   

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

for j=1:length(Tracks)
    if contains( trackname,'NT')
        if  contains( trackname,'2 mg') || contains( trackname,'2_mg') 
            for i= 4:4:size(Tracks{j},1)-4
                NT_vel_2{Tracks{j}(i,39)}=[NT_vel_2{Tracks{j}(i,39)}, Tracks{j}(i,32)];
            end
        elseif contains( trackname,'6 mg') || contains( trackname,'6_mg') 
            for i= 4:4:size(Tracks{j},1)-4
                NT_vel_6{Tracks{j}(i,39)}=[NT_vel_6{Tracks{j}(i,39)}, Tracks{j}(i,32)];
            end
        end
    elseif contains( trackname,'CDH')
        if  contains( trackname,'2 mg') || contains( trackname,'2_mg') 
            for i= 4:4:size(Tracks{j},1)-4
                CDH1_vel_2{Tracks{j}(i,39)}=[CDH1_vel_2{Tracks{j}(i,39)}, Tracks{j}(i,32)];
            end
        elseif contains( trackname,'6 mg') || contains( trackname,'6_mg') 
            for i= 4:4:size(Tracks{j},1)-4
                CDH1_vel_6{Tracks{j}(i,39)}=[CDH1_vel_6{Tracks{j}(i,39)}, Tracks{j}(i,32)];
            end
        end
    end
    
end

end

my_boxplot(NT_vel_6{2}, CDH1_vel_6{2},'shNT 6 mg/ml', 'shCDH1 6 mg/ml')

my_boxplot(NT_vel_6{2}, CDH1_vel_6{2},NT_vel_6{3}, CDH1_vel_6{3},'shNT 6 mg/ml Mid', 'shCDH1 6 mg/ml Mid','shNT 6 mg/ml Edge', 'shCDH1 6 mg/ml Edge')

[h,p] = kstest2(NT_vel_6{3},CDH1_vel_6{3});



