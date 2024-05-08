function [ neworder, inFocus ] = getNewOrder( inFocus )
%GETNEWORDER get the correct order to place consecutive planes together
%   Detailed explanation goes here
    list = [inFocus.zpos];
    list(9:16) = [];
    ordZpos = sort([list],'descend');%Descending order because the 
    %first plane to get into focus is the plane that is the highest in Z.

    for i = 1:length(list)
        x = ordZpos(i);
        Positions = find([inFocus.zpos] == x);
        for j = 1:length(Positions)
            z = Positions(j)
            inFocus(z).globalch = i
        end
    end
     
        % inFocus(i).globalch = x;
        % inFocus([inFocus.zpos] == ordZpos(i)).globalch = i; 
    
    focus = cat(1,inFocus.zpos);
    [~,neworder] = sort(focus,'descend');
    
end
