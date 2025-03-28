function [MSD,D] = calc(coord)

dim = size(coord,2);

switch dim
    case 1
        coord(:,2) = 0;
    case 2        
    otherwise
        error('unexpected dimension for the vector')
end

%%%Fist calculate angle difference between each frame

Dtheta = diff(coord(:,1).');
Dphi = diff(coord(:,2).');

Dtheta = mod(Dtheta + pi, 2*pi) - pi;
Dphi = mod(Dphi + pi, 2*pi) - pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MSD = zeros(size(coord,1)-1,1);
%Calculate mean-squere-displacement
for i = 1:size(coord,1)-1
    
    stp = i;
    cnt =  1;
    D1  = [];
    while cnt<=stp && cnt+stp<=size(coord,1)
        
        idx = cnt:stp:size(coord,1);
        DX  = diff(coord(idx,1).');
        DY  = diff(coord(idx,2).');
        DZ  = diff(coord(idx,3).');
        D1  = [D1 sqrt(DX.^2 + DY.^2 + DZ.^2)];
        cnt = cnt+1;
        
        if ~isempty(D1)
            
            D2=D1(~isnan(D1));
            
            if ~isempty(D2)
                
                MSD(i) = mean(D2.^2);
                
            else
                
                MSD(i) = NaN;
                
            end
        end
    end %while
end

D   = D(:);
MSD = MSD(:);

end

