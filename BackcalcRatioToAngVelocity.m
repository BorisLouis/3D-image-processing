for q = 1:size(CommonTraces, 1)
    try
        ratio(:,q) = CommonTraces.Ratio{q,1};
        time = [1:1:size(ratio,1)]*0.01;
        Height = 1.083133;
    
        TotalInt = mean(CommonTraces.Int1{q,1} + CommonTraces.Int2{q,1});  
        Amplitude = 0.0009372*TotalInt - 0.06379;
        % Amplitude = 0.1;
    
        for i = 1:size(ratio, 1)
            % if i ~= 1
            %     CandidateAngles(1,1) = real(asin((ratio(i,1) - Height)./Amplitude));
            %     CandidateAngles(1,2) = pi - real(asin((ratio(i,1) - Height)./Amplitude)); 
            %     [~, idx] = min(abs(CandidateAngles - AnglesTest(i-1,1)));
            %     AnglesTest(i,1) = CandidateAngles(idx);
            % elseif i == 1
                AnglesTest(i,q) = real(asin((ratio(i,q) - Height)./Amplitude)./1.90);
            % end
        end
    
        % Theoretical_Angles = [0:0.001:pi./2];
        % TheoreticalRatio = Height + Amplitude*cos(1.90*Theoretical_Angles);
        % 
        % for z = 1:size(ratio,1)
        %     [~, idx] = min(abs(TheoreticalRatio - ratio(z)));
        %     Angles(z,1) = Theoretical_Angles(1,idx);
        % end
        % 
        theta_raw = AnglesTest(:,q);
        
        tau_values = 1:floor(length(theta_raw)/2);
        % MSD = zeros(size(tau_values));
        
        for i = 1:length(tau_values)
            tau = tau_values(i);
            displacements = (theta_raw(1+tau:end) - theta_raw(1:end-tau)).^2;
            displacements(displacements == 0) = [];
            MSD(q,i) = mean(displacements);
        end
        
        time_lags = tau_values * mean(diff(time));
        p = fit(time_lags(1, 33:100).', MSD(q, 33:100).', 'a*x+b');
        coeff = coeffvalues(p);
        D_theta = coeff(1) / 2;
        angular_speed(q,1) = (real(sqrt(2 * D_theta)*(180/pi)))./2;
    catch
        angular_speed(q,1) = NaN;
    end
end