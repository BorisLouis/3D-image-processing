maxIter = 5000;
Height = 0;
Amplitude = 0.777915;

for q = 1:size(CommonTraces, 1)
    ratio(:,q) = CommonTraces.Difference{q,1};
    time = [1:1:size(ratio,1)]*0.01;
    TotalInt = mean(CommonTraces.Int1{q,1} + CommonTraces.Int2{q,1});  
    
    for i = 1:size(ratio, 1)
        AnglesTest(i, q) = real(asin((ratio(i, q) - Height) ./ Amplitude) ./ 1.74);
    end
    theta_raw = AnglesTest(:, q);
    tau_values = 1:floor(length(theta_raw) / 2);

    for i = 1:length(tau_values)
        tau = tau_values(i);
        displacements = (theta_raw(1 + tau:end) - theta_raw(1:end - tau)).^2;
        displacements(displacements == 0) = [];
        MSD(q, i) = mean(displacements);
    end

    time_lags = tau_values * mean(diff(time));
    p = fit(time_lags(1, 33:100).', MSD(q, 33:100).', 'a*x+b');
    coeff = coeffvalues(p);
    D_theta = coeff(1) / 2;
    angular_speed(q, 1) = (real(sqrt(2 * D_theta) * (180 / pi))) / 2;

end
