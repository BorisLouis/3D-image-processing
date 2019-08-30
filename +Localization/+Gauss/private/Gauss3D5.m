function F = Gauss3D5(x,data)

    F = x(1)./sqrt(2*pi*x(2)).*...
        exp(-(((data(:,:,:,1) - x(5)) .^2) ./ (2 * x(2) .^2))-...
             (((data(:,:,:,2) - x(6)) .^2) ./ (2 * x(2) .^2))-...
             (((data(:,:,:,3) - x(7)) .^2) ./ (2 * x(3) .^2))) +...
        x(1)./sqrt(2*pi*x(2)).*...
        exp(-(((data(:,:,:,1) - x(8)) .^2) ./ (2 * x(2) .^2))-...
             (((data(:,:,:,2) - x(9)) .^2) ./ (2 * x(2) .^2))-...
             (((data(:,:,:,3) - x(10)) .^2) ./ (2 * x(3) .^2))) + ...
        x(1)./sqrt(2*pi*x(2)).*...
        exp(-(((data(:,:,:,1) - x(11)) .^2) ./ (2 * x(2) .^2))-...
             (((data(:,:,:,2) - x(12)) .^2) ./ (2 * x(2) .^2))-...
             (((data(:,:,:,3) - x(13)) .^2) ./ (2 * x(3) .^2))) +...
        x(1)./sqrt(2*pi*x(2)).*...
        exp(-(((data(:,:,:,1) - x(14)) .^2) ./ (2 * x(2) .^2))-...
             (((data(:,:,:,2) - x(15)) .^2) ./ (2 * x(2) .^2))-...
             (((data(:,:,:,3) - x(16)) .^2) ./ (2 * x(3) .^2))) + ...
        x(1)./sqrt(2*pi*x(2)).*...
        exp(-(((data(:,:,:,1) - x(17)) .^2) ./ (2 * x(2) .^2))-...
             (((data(:,:,:,2) - x(18)) .^2) ./ (2 * x(2) .^2))-...
             (((data(:,:,:,3) - x(19)) .^2) ./ (2 * x(3) .^2))) + x(4);


end
