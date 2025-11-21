function [sacf] = SACF(Frame)
    I = double(Frame);
    [ny, nx] = size(I);

    % Step 1: spatial mean
    mu = mean(I(:));

    % Step 2: fluctuations
    dI = I - mu;

    % Prepare output: autocorrelation for each spatial lag
    % lag_x = -(nx-1) : (nx-1)
    % lag_y = -(ny-1) : (ny-1)
    r = zeros(2*ny-1, 2*nx-1);

    % Compute all lags
    for eta = -(ny-1):(ny-1)
        for xi = -(nx-1):(nx-1)

            % Overlapping region for this shift
            x1 = max(1, 1+xi) : min(nx, nx+xi);
            x2 = x1 - xi;

            y1 = max(1, 1+eta): min(ny, ny+eta);
            y2 = y1 - eta;

            % numerator: average of product of fluctuations
            num = mean( dI(y1, x1) .* dI(y2, x2), 'all' );

            % denominator: mean^2
            sacf(eta+ny, xi+nx) = num / (mu * mu);
        end
    end
end

