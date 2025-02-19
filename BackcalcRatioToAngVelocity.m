for q = 1:size(CommonTraces, 1)
    ratio = CommonTraces.Ratio{q,1};
    time = [1:1:size(ratio,1)]*0.01;
    Height = 1.083133;

    TotalInt = mean(CommonTraces.Int1{q,1} + CommonTraces.Int2{q,1});  
    Amplitude = 0.0009372*TotalInt - 0.06379;

    Theoretical_Angles = [0:0.001:pi./2];
    TheoreticalRatio = Height + Amplitude*cos(1.90*Theoretical_Angles);
    
    for z = 1:size(ratio,1)
        [~, idx] = min(abs(TheoreticalRatio - ratio(z)));
        Angles(z,1) = Theoretical_Angles(1,idx);
    end
    
    theta_raw = Angles(56:130)
    
    tau_values = 1:floor(length(theta_raw)/2);
    MSD = zeros(size(tau_values));
    
    for i = 1:length(tau_values)
        tau = tau_values(i);
        displacements = (theta_raw(1+tau:end) - theta_raw(1:end-tau)).^2;
        MSD(i) = mean(displacements);
    end
    
    time_lags = tau_values * mean(diff(time));
    p = fit(time_lags(1, 1:35).', MSD(1, 1:35).', 'a*x');
    coeff = coeffvalues(p);
    D_theta = coeff(1) / 2;
    angular_speed(q,1) = real(sqrt(2 * D_theta)*(180/pi));
end