function NewFrameTrackCoordinates = EliminateDuplicateDataPoints( FrameTracks, FrameTrackCoordinates, Tracks )
% for unknown reasons, there are duplicate positions in some frames ...
% lets try to replace them with an interpolation of the
% last/next track positions
%
% The positions are only replaced in FrameTrackCoordinates - original
% tracks untouched!!!!!!!!!!!!

	for f = 1:length(FrameTracks)
		
		% get UNIQUE cell positions - keep the order stable
		[C, ia, ic] = unique(FrameTrackCoordinates{f}, 'rows', 'stable');
		
		% check if unique positions are less than the original ones
		if length(C) ~= length(FrameTrackCoordinates{f})
			disp(['duplicate data in frame ' num2str(f)]);
			
			% find where the removed (duplicate) position is
			removed = setdiff(1:length(FrameTrackCoordinates{f})', ia );	

			% now find not only the removed one but also the remaining ones
			doublette = [];
			for r = 1:length(removed)
				doublette(end+1:end+2) = find( FrameTrackCoordinates{f}(:,1) == FrameTrackCoordinates{f}(removed(r),1) & ...
					FrameTrackCoordinates{f}(:,2) == FrameTrackCoordinates{f}(removed(r),2));
			end
			
			
			% replace the doublette positions by interpolated ones from the
			% corresponding tracks
			CorrespondingTracks = FrameTracks{f}(doublette);
			for j=1:length(CorrespondingTracks)
				
				% find the corresponding row in the track
				index	= find( Tracks{CorrespondingTracks(j)}(:,1) == f-1 );
				
				% interpolate new position
				if index>1 && index<length(Tracks{CorrespondingTracks(j)}(:,1))
					NewPosition = mean( Tracks{CorrespondingTracks(j)}([index-1 index+1],2:end) );
				elseif index==1 % cant interpolate, take position '2'
					NewPosition = Tracks{CorrespondingTracks(j)}(2,2:end);
				else % cant interpolate, take position 'end-1'
					NewPosition = Tracks{CorrespondingTracks(j)}(end-1,2:end);				
                end
                
                % if new position happens to be identical to the old one:
                % adjust it !! (try future positions of track first - 
                % if track is ending try older postions)
				lauf=1;
                no_endless_loops=true; % a switch to avoid an endless loop of +1 -1 +1 -1 +1 -1
                while NewPosition == FrameTrackCoordinates{f}(doublette(j),:)
                    if index+lauf<=size(Tracks{CorrespondingTracks(j)},1) && no_endless_loops
                        NewPosition=Tracks{CorrespondingTracks(j)}( index+lauf,2:end);
                        lauf=lauf+1;
                    else
                        no_endless_loops=false;
                        lauf=lauf-1;
                        NewPosition=Tracks{CorrespondingTracks(j)}( index+lauf,2:end);
                    end
                end
                
                % check if the other doublett happens to have the same
                % interpolatet position => if yes ... just keep the old one
                
                if mod(j,2)==1
                    otherNewPosition=NewPosition;
                else
                    if NewPosition==otherNewPosition
                        NewPosition=FrameTrackCoordinates{f}(doublette(j),:);
                    end
                end
                
				% replace the position in FrameTrack coordinates with this new
				% position
				disp( ['replace: ' num2str(FrameTrackCoordinates{f}(doublette(j),:)) ' -> ' num2str(NewPosition)]);
				%disp();
				%disp(NewPosition);
				FrameTrackCoordinates{f}(doublette(j),:) = NewPosition;
				
			end
			disp('positions replaced by interpolations.');
			
		end		
	end

	NewFrameTrackCoordinates = FrameTrackCoordinates;
end