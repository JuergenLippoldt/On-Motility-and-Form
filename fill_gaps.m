function [Tracks]=fill_gaps(Tracks)
waitbar(0, 'filling graps...');
for lauf = 1:length(Tracks)
    if (Tracks{lauf}(end,1)-Tracks{lauf}(1,1)) > length(Tracks{lauf})-1
        interim=Tracks{lauf};
        for  lauf2=2:length(Tracks{lauf})
            if (Tracks{lauf}(lauf2,1)-Tracks{lauf}(lauf2-1,1)) > 1
                blub2=interp1([Tracks{lauf}(lauf2-1,1),Tracks{lauf}(lauf2,1)],[Tracks{lauf}(lauf2-1,2),Tracks{lauf}(lauf2,2)],Tracks{lauf}(lauf2-1,1):1:Tracks{lauf}(lauf2,1));
                blub3=interp1([Tracks{lauf}(lauf2-1,1),Tracks{lauf}(lauf2,1)],[Tracks{lauf}(lauf2-1,3),Tracks{lauf}(lauf2,3)],Tracks{lauf}(lauf2-1,1):1:Tracks{lauf}(lauf2,1));
                interim=[interim(1:find(interim(:,1)==Tracks{lauf}(lauf2-1,1)),:);[Tracks{lauf}(lauf2-1,1)+1:1:Tracks{lauf}(lauf2,1)-1;blub2(2:end-1);blub3(2:end-1)]';interim(find(interim(:,1)==Tracks{lauf}(lauf2,1)):end,:)];%(Tracks{lauf}(lauf2,1)-Tracks{lauf}(1,1):Tracks{lauf}(end,1)-Tracks{lauf}(1,1)+Tracks{lauf}(lauf2,1)-Tracks{lauf}(lauf2-1,1)-1,:)=interim(Tracks{lauf}(lauf2-1,1)-Tracks{lauf}(1,1):Tracks{lauf}(end,1)-Tracks{lauf}(1,1)-1,:);
                
 %               interim (Tracks{lauf}(lauf2-1,1)+1:Tracks{lauf}(lauf2,1)-1,:)=[Tracks{lauf}(lauf2-1,1)+1:1:Tracks{lauf}(lauf2,1)-1,round(blub2(2:end-1)),round(blub3(2:end-1))];   
            end
 
        end
        for lauf3=2:length(interim)
            if interim(lauf3,1)~=interim(lauf3-1,1)+1
                error('gap-filling went wrong')
            end
        end
        Tracks{lauf}=interim;
    end

    waitbar( lauf/length(Tracks));
end
close(waitbar(0));
end