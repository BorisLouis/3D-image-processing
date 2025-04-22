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
    [f, gov]    = fit(msad(2,:)',msad(1,:)','a*(1-exp(-c*x))');
   
    % figure()
    % plot(f, msad(2,:)',msad(1,:)');
    % % % plot(f, tau',tofit')
    
    if gov.rsquare > 0.85
        g = coeffvalues(f);
        Dr = g(2)/div;
    else
        Dr = NaN;
    end
end