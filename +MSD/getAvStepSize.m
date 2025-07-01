function [AvStepSize] = getAvStepSize(Coord)
    dim = size(Coord, 2);
    for i = 1:size(Coord, 1)-1
        if dim == 1
            Step(i) = sqrt((Coord(i+1,1) - Coord(i,1)).^2);
        elseif dim == 2
            Step(i) = sqrt((Coord(i+1,1) - Coord(i,1)).^2 + (Coord(i+1,2) - Coord(i,2)).^2);
        elseif dim == 3
            Step(i) = sqrt((Coord(i+1,1) - Coord(i,1)).^2 + (Coord(i+1,2) - Coord(i,2)).^2 + (Coord(i+1,3) - Coord(i,3)).^2);
        end
    end
    AvStepSize = mean(Step);
end

