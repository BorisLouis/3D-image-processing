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
            SelectTheta = coord(idx,1).';
            SelectPhi = coord(idx,2).';
            TimeLagcorr = abs(diff(tau(idx)) - dt*expTime) < 0.002;

            for z = 1:size(SelectTheta, 2)-1
                if TimeLagcorr(z) == 1
                    theta1 = SelectTheta(z);
                    phi1 = SelectPhi(z);
                    r1 = [cos(theta1)*cos(phi1) cos(theta1)*sin(phi1) sin(theta1)];
    
                    theta2 = SelectTheta(z+1);
                    phi2 = SelectPhi(z+1);
                    r2 = [cos(theta2)*cos(phi2) cos(theta2)*sin(phi2) sin(theta2)];

                    AngDisp= [AngDisp; real(acos(dot(r1, r2)))];
                else
                    continue
                end
            end
            cnt = cnt + 1;
        end

        if ~isempty(AngDisp)
            Dr = AngDisp(~isnan(AngDisp));
            Dr(Dr == 0) = [];
            
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
