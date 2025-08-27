function [Model] = TestModels(Gr, tau)
    single_exp    = @(p,t) (p(1) * exp(-t./p(2))).';
    bi_exp        = @(p,t) (p(1) * exp(-t./p(2)) + (1-p(1)) * exp(-t./p(3))).';
    stretched_exp = @(p,t) (p(1) * exp(-(t./p(2)).^p(3))).';
    
    opts = optimset('Display','off');
    Nfit = round(length(tau)/3);
    tfit = tau(1:Nfit);
    yfit = Gr(1:Nfit);
    
    p0 = [1, mean(tfit)];
    lb = [0,0]; ub=[2,Inf];
    p_single = lsqcurvefit(single_exp,p0,tfit,yfit,lb,ub,opts);
    
    p0 = [0.5, mean(tfit)/2, mean(tfit)*2];
    lb = [0,0,0]; ub = [1,Inf,Inf];
    p_bi = lsqcurvefit(bi_exp,p0,tfit,yfit,lb,ub,opts);
    
    p0 = [1, mean(tfit), 0.8];
    lb = [0,0,0]; ub=[2,Inf,1];
    p_stretch = lsqcurvefit(stretched_exp,p0,tfit,yfit,lb,ub,opts);
    
    models = {'Single','Bi','Stretched'};
    params = {p_single,p_bi,p_stretch};
    funcs  = {single_exp,bi_exp,stretched_exp};
    
    N = numel(yfit);
    results = table('Size',[3 4],'VariableTypes',{'double','double','double','double'},...
    'VariableNames',{'Chi2red','AIC','BIC','RSS'},'RowNames',models);
    
    for i=1:3
        y_model = funcs{i}(params{i},tfit);
        resid   = yfit - y_model;
        k = numel(params{i});
        RSS = sum(resid.^2);
        chi2red = RSS/(N-k);
        AIC = 2*k + N*log(RSS/N);
        BIC = k*log(N) + N*log(RSS/N);
        results{i,:} = [chi2red,AIC,BIC,RSS];
    end

    [~, idx] = min(results.Chi2red);
    if idx == 1
        Model = 'Single Exponential';
    else
        if idx == 2 %Best idx for biexponential
            if abs(results.BIC(2)) - abs(results.BIC(1)) < 2
                Model = 'Single Exponential';
            else
                Model = 'Bi Exponential';
            end
        elseif idx == 3
            if abs(results.BIC(3)) - abs(results.BIC(1)) < 2
                Model = 'Single Exponential';
            else
                Model = 'Stretched Exponential';
            end
        end
    end
end

