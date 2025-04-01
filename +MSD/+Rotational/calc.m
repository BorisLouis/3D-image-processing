function [MSAD,D] = calc(coord)

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

MSAD = zeros(size(coord,1)-1,1);
%Calculate mean-square-displacement
for i = 1:size(coord,1)-1
    
    stp = i;
    cnt =  1;
    Dr  = [];
    while cnt<=stp && cnt+stp<=size(coord,1)
        
        idx = cnt:stp:size(coord,1);
        Dtheta  = diff(coord(idx,1).');
        Dphi  = diff(coord(idx,2).');
        Dr  = [Dr sqrt(Dtheta.^2 + Dphi.^2)];
        cnt = cnt+1;
        
        if ~isempty(Dr)
            
            D=Dr(~isnan(Dr));
            
            if ~isempty(D)
                
                MSAD(i) = mean(D.^2);
                
            else
                
                MSAD(i) = NaN;
                
            end
        end
    end %while
end

MSAD = MSAD(:);

end

