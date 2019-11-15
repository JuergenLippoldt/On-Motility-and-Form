function [FrameTracks, FrameTrackCoordinates] = FindAllTracksInFrames( Tracks )
% find within Tracks all that lie in frame FrameNr	
	
	FrameTracks = cell( max(cellfun( @(x) max(x(:,1)), Tracks )) + 1, 1 );
	FrameTrackCoordinates = FrameTracks;
	
	progress = waitbar(0, 'rumgichteln...');
	
	for j=1:length( Tracks )
		ThisTrack	= Tracks{j};
		IsInFrames	= ThisTrack(:,1);
		
		if ~isnan(IsInFrames) % track might be "deleted" by replacing with NaN
			for f=1:length(IsInFrames)
				FrameTracks{IsInFrames(f)+1}(end+1) = j;
				FrameTrackCoordinates{IsInFrames(f)+1}(end+1,:) = ThisTrack(f,2:end);
			end
		end
		
		
		waitbar(j/length(Tracks), progress);
		
	end
	
	
	% we must call UNIQUE on the data or else we will get problems
	% later with voronoin().
	% not clear why it is not unique in the first place
	%for f=1:length(FrameTracks)
	%	unique(FrameTrackCoordinates, 'rows');

	
	try delete(progress); catch; end

end