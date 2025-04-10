function [MSD] = calc2(coord, tau, expTime)

    dim = size(coord,2);
    
    switch dim
        case 1
            coord(:,2) = 0;
        case 2        
        otherwise
            error('unexpected dimension for the vector')
    end
    
    MSAD = zeros(size(coord,1)-1,1);
    %Calculate mean-square-displacement
    maxIdx = round((max(tau)-min(tau))/expTime);
    AngDisp = [];
    for dt = 1:size(coord,1)-1
        
        cnt =  1;
        
        while cnt<=dt && cnt+dt<=size(coord,1)
            
            idx = cnt:dt:size(coord,1);
            DTheta = diff(coord(idx,1)).';
            DPhi = diff(coord(idx,2)).';
            TimeLagcorr = abs(diff(tau(idx)) - dt*expTime) < 0.002;

            DTheta(TimeLagcorr == 0) = [];
            DPhi(TimeLagcorr == 0) = [];

            AngDisp= [AngDisp; (sqrt((DTheta).^2 + (DPhi).^2))'];

            cnt = cnt + 1;
        end

        if ~isempty(AngDisp)
            Dr = AngDisp(~isnan(AngDisp));
            
            if ~isempty(Dr)
                MSAD(dt) = mean(Dr.^2);
            else
                MSAD(dt) = NaN;
            end
        else
            MSAD(dt) = NaN;
        end
        TimeLag(dt) = dt*expTime;
    end

    MSD = [MSAD'; TimeLag];
end