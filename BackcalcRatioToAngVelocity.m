maxIter = 5000;
Height = 1.08;
tolX = 1e-6;     
tolFun = 1e-6; 
options = optimset('MaxIter', maxIter, 'TolX', tolX, 'TolFun', tolFun, 'Display', 'iter');

numTrials = 5;
best_sigma2 = Inf;
best_params = [];

for trial = 1:numTrials
    Slope0 = 0.0005 + 0.001 * rand;
    Intersect0 = -0.1 + 0.1 * rand;
    params0 = [Slope0; Intersect0];
    cost_func = @(params) cost_function(params, CommonTraces, Height);
    if exist('fminunc', 'file')
        options_fminunc = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'MaxIterations', maxIter, 'TolFun', tolFun);
        [params_opt, fval] = fminunc(cost_func, params0, options_fminunc);
    else
        [params_opt, fval] = fminsearch(cost_func, params0, options);
    end

    if fval < best_sigma2
        best_sigma2 = fval;
        best_params = params_opt;
    end
end

params_opt = best_params;
disp("Best parameters found:");
disp(params_opt);

Slope_opt = params_opt(1);
Intersect_opt = params_opt(2);

for q = 1:size(CommonTraces, 1)
    ratio(:,q) = CommonTraces.Ratio{q,1};
    time = [1:1:size(ratio,1)]*0.01;
    TotalInt = mean(CommonTraces.Int1{q,1} + CommonTraces.Int2{q,1});  
    Amplitude = Slope_opt(1) * TotalInt + Intersect_opt(1);
    
    for i = 1:size(ratio, 1)
        AnglesTest(i, q) = real(asin((ratio(i, q) - Height) ./ Amplitude) ./ 1.90);
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

function sigma2 = cost_function(params, CommonTraces, Height)
    for q = 1:size(CommonTraces, 1)
        ratio(:,q) = CommonTraces.Ratio{q,1};
        time = [1:1:size(ratio,1)]*0.01;
        TotalInt = mean(CommonTraces.Int1{q,1} + CommonTraces.Int2{q,1});  
        Amplitude = params(1) * TotalInt + params(2);
        
            for i = 1:size(ratio, 1)
                AnglesTest(i, q) = real(asin((ratio(i, q) - Height) ./ Amplitude) ./ 1.90);
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

    Theoretical_difference = abs(angular_speed - 25);
    sigma2 = mean(Theoretical_difference);  % The function must return a single scalar value
end
