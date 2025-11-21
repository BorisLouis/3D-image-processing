function [omega0] = fitSACF(r_map)
    [nx, ny] = size(r_map);
    cx = floor(nx/2) + 1;
    cy = floor(ny/2) + 1;
    
    [xi, eta] = meshgrid(1:ny, 1:nx);
    xi = xi - cy;     % shift to center at zero
    eta = eta - cx;
    
    rho = sqrt(xi.^2 + eta.^2);   % radial distance
    
    % -------------------------------------------------------------
    % 2. Flatten arrays for easier processing
    % -------------------------------------------------------------
    rho_vec = rho(:);
    r_vec   = r_map(:);
    
    % remove NaNs if present
    valid = ~isnan(r_vec);
    rho_vec = rho_vec(valid);
    r_vec   = r_vec(valid);
    
    % -------------------------------------------------------------
    % 3. Define radial bins
    % -------------------------------------------------------------
    dr = 1;   % bin width in pixels
    r_max = max(rho_vec);
    edges = 0:dr:r_max;
    bin_centers = edges(1:end-1) + dr/2;
    
    r_rad = zeros(size(bin_centers));
    N_bin = zeros(size(bin_centers));
    
    for k = 1:length(bin_centers)
        mask = (rho_vec >= edges(k)) & (rho_vec < edges(k+1));
        r_rad(k) = mean(r_vec(mask));
        N_bin(k) = sum(mask);
    end
    
    % Remove empty bins
    nonempty = N_bin > 0;
    rho_fit = bin_centers(nonempty);
    r_fit   = r_rad(nonempty);
    w       = sqrt(N_bin(nonempty));   % weights = sqrt(counts)
    
    % -------------------------------------------------------------
    % 4. Initial parameter estimates
    % -------------------------------------------------------------
    Ginf0 = mean(r_fit(end-3:end));             % asymptote guess
    G0_0  = r_fit(1) - Ginf0;                   % amplitude guess
    omega0_0 = rho_fit(find(r_fit <= Ginf0 + G0_0*exp(-1),1));
    if isempty(omega0_0)
        omega0_0 = 2;  % fallback guess
    end
    
    p0 = [G0_0, omega0_0, Ginf0];   % initial parameters
    
    % -------------------------------------------------------------
    % 5. Define model function
    % -------------------------------------------------------------
    model_fun = @(p, rho) p(1)*exp(-(rho.^2)/(p(2)^2)) + p(3);
    
    % -------------------------------------------------------------
    % 6. Perform nonlinear least squares fit (weighted)
    % -------------------------------------------------------------
    opts = optimoptions('lsqcurvefit','Display','off');
    
    lb = [-Inf, 0, -Inf];   % omega0 must be positive
    ub = [ Inf, Inf, Inf];
    
    p_fit = lsqcurvefit(@(p, rho) w .* model_fun(p,rho), ...
                        p0, rho_fit, w .* r_fit, lb, ub, opts);
    
    G0     = p_fit(1);
    omega0 = p_fit(2);
    Ginf   = p_fit(3);
    
    % -------------------------------------------------------------
    % 7. Plot results
    % -------------------------------------------------------------
    figure; hold on;
    plot(rho_fit, r_fit, 'ko', 'MarkerFaceColor','k','DisplayName','Radial Average');
    rho_plot = linspace(0, max(rho_fit), 200);
    plot(rho_plot, model_fun(p_fit, rho_plot), 'r-', 'LineWidth',2, ...
         'DisplayName','Fit');
    
    xlabel('\rho (pixels)');
    ylabel('r(\rho)');
    legend();
    title(sprintf('Fit Result:  \\omega_0 = %.3f pixels', omega0));
    grid on;
end

