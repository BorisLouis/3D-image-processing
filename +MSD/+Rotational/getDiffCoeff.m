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
    f     = fit(tau(:),tofit(:),'a*x+b');
    
    g = coeffvalues(f);
    Dr = g(1)/div;
end