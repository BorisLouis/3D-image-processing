function [Dr, correctionParams] = GetDiffusion(Gr, tau, Radius1, Temp, Dimension, ExpModel, Plot, eta_truth, partType, Angle)
    if strcmp(ExpModel, 'Single Exponential')
        expfunction = @(p,t) (exp(-t./p(1))).';
    elseif strcmp(Model, 'Bi Exponential')
        expfunction = @(p,t) (exp(-t./p(1)) + 1 * exp(-t./p(2))).';
    elseif strcmp(Model, 'Stretched Exponential')
        expfunction = @(p,t) (exp(-(t./p(1)).^p(2))).';
    else
        error('Please indicate right model for exponential fit')
    end
                    
    opts = optimset('Display','off');
    Nfit = 10;%round(length(tau)/3);
    tfit = tau(1:Nfit);
    yfit = Gr(1:Nfit);
    
    p0 = [tau(1, find(islocalmin(abs(Gr-(1./exp(1)))) == 1, 1, 'first'))];
    lb = [0,0]; ub=[2,Inf];
    p_single = lsqcurvefit(expfunction,p0,tfit,yfit,lb,ub,opts);

    if strcmp(ExpModel, 'Single Exponential')
        TauC = p_single(1);
    elseif strcmp(Model, 'Bi Exponential')
        TauC = max(p_single(1), p_single(2));
    elseif strcmp(Model, 'Stretched Exponential')
        TauC = p_single(1);
    else
        error('Please indicate right model for exponential fit')
    end

    alpha = 0.00300197 - exp(-0.0417047*eta_truth -3.08577);
    beta = 0.0239624 - 0.00649547*log(eta_truth);
    C = 6.54169 - exp(-0.0352001*eta_truth + 4.03419);
    correctionParams = [];

    TauC = TauC .* (alpha * TauC.^beta + C);
    % if nargin > 7 && ~isempty(eta_truth)   % only apply if ground truth provided
    %     % convert TauC -> viscosity (uncorrected)
    %     Dr_uncorr = 1 ./ (6 * TauC);  % simple relation (3D)
    %     eta_uncorr = MSD.Rotational.getViscosity(Dr_uncorr, Radius1, partType, Temp, Angle);
    % 
    %     % correction function
    %     corrFun = @(p, TauC) TauC .* (p(1) * TauC.^p(2) + p(3));
    % 
    %     % cost function vs ground truth viscosity
    %     costFun = @(p) sum((MSD.Rotational.getViscosity(1 ./ (6 * corrFun(p, TauC)), Radius1, partType, Temp, Angle) - eta_truth).^2);
    % 
    %     % fit correction parameters
    %     p0_corr = [0, 0, 1]; 
    %     bestP = fminsearch(costFun, p0_corr, opts);
    % 
    %     % apply correction
    %     TauC = corrFun(bestP, TauC);
    % 
    %     % store parameters
    %     correctionParams = [bestP(1),bestP(2), bestP(3), eta_uncorr, eta_truth];
    % else
    %     correctionParams = [];
    % end
        
    if strcmp(Dimension, '3D')
        Dr = (((((1./(2*pi*TauC))*2*pi))./8));
    elseif strcmp(Dimension, '2D')
        Dr = (((((1./(2*pi*TauC))*2*pi))./4));
    end
    %nR = (3*1.380649*10^-23*Temp*log(Radius1(1)./Radius1(2)))./(pi*Dr*(Radius1(1)*10^(-9))^3)*10^3;

    if Plot == 1
        figure;
        plot(tau',Gr,'ko','DisplayName','Autocorrelation'); hold on
        plot(tfit,expfunction(p_single,tfit),'r-','LineWidth',2,'DisplayName','Single exp fit');
        xlabel('\tau (s)'); ylabel('G_r(\tau)');
        legend; 
        title(sprintf('Single exponential fit (\\tau = %.3f s)', TauC));
    end
                    
end

