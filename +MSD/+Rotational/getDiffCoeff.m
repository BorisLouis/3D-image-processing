function Dr = getDiffCoeff(msad,tau,fitRange,dim)
    
    switch dim
        case '1D'
            error('Cannot do rotational tracking in 1D')
        case '2D'
            div = 2;
        case '3D'
            div = 4;
        otherwise
            error('Unknown dim, dim needs to be provided as 1D 2D or 3D')
    end

    % assert(min(size(msad))==1,'MSD needs to be provided as a vector')
    % assert(and(fitRange<=1,isnumeric(fitRange)),'fit Range needs to be numerical between 0 and 1');
    
    tofit = msad(1, 1:fitRange);
    tau   = msad(2, 1:fitRange);

    % [f, gov]    = fit(msad(2,:)',msad(1,:)','a*(1-(1-b^2)*exp(-c*x))');
    try
        Lower = [0, 0, 0];
        Upper = [2*max(msad(1,:)), sqrt(msad(1,1)/(2*max(msad(1,:)))), 100];
        if div == 4
            [f, gov]    = fit(msad(2,:)',msad(1,:)','a*(1-(1-b^2)*exp(-(1.6*c*x)^0.95))', 'Lower', Lower, 'Upper', Upper);
        elseif div == 2
            [f, gov]    = fit(msad(2,:)',msad(1,:)','a*(1-(1-b^2)*exp(-(4*c*x)))');%, 'Lower', Lower, 'Upper', Upper);
        end
        figure()
        plot(f, msad(2,:)',msad(1,:)');
        % 
        % 
        if gov.rsquare > 0.75
            g = coeffvalues(f);
            Dr = g(3);
        else
            Dr = NaN;
        end
    catch
        Dr = NaN;
    end
end