
function [QP,mask] = getQPFullPhase(stack,s,mask)
    
    % mirror the data and compute adequate Fourier space grid
    [stackM,kx,kz] = QP_package.getMirroredStack(stack,s);
    
    if nargin < 3 % if no mask are provided
    % compute usefull stuff
    th = asin(s.optics.NA/s.optics.n);
    th_ill = asin(s.optics.NA_ill/s.optics.n);
    k0max = s.optics.n*2*pi/(s.optics.lambda - s.optics.dlambda/2);
    k0min = s.optics.n*2*pi/(s.optics.lambda + s.optics.dlambda/2);
    
    % compute Fourier space grid and the phase mask
    [Kx,Kz] = meshgrid(kx,kz);
    if isempty(s.optics.kzT)
        mask2D = Kz > k0max*(1-cos(th_ill));
    else
        mask2D = Kz > s.optics.kzT;
    end
    
    if strcmp(s.proc.applyFourierMask, 'true') % => compute the CTF mask for extra denoising
    % CTF theory 
    maskCTF = ((Kx-k0max*sin(th_ill)).^2 + (Kz-k0max*cos(th_ill)).^2 <= k0max^2) & ...
              ((Kx+k0min*sin(th_ill)).^2 + (Kz-k0min).^2 >= k0min.^2) & ...
              Kx>=0 & ...
              Kz < k0max*(1-cos(th));
    maskCTF = maskCTF | maskCTF(:,end:-1:1);
    mask2D = mask2D & maskCTF;
    end
    
    % since we assume a circular symetric CTF, we expand the 2Dmask in 3D
    mask = QP_package.map3D(mask2D);
     
    end
    
    % Cross-Spectral Density calculation
    Ik = fftshift(fftn(fftshift(stackM)));
    %EDIT BORIS based on : 
    % Main file for HBDTI
    % Version 1.0 -
    % Related Reference:
    % last modified on 06/17/2022
    % by  Linpeng Lu, and Chao Zuo (zuochao@njust.edu.cn)
    
    % Amplitude = abs(ifft2(fftshift(Freal)));
    % Phase= angle(ifft2(fftshift(Freal)));
    % mask = QP_package.cropXY(mask, max(size(Ik)));
    
    Gamma = Ik.*mask; % cross-spectral density
    
    csd = ifftshift(ifftn(ifftshift(Gamma)));
    csd = csd(1:size(stack,1),1:size(stack,2),1:size(stack,3));
    %csd = csd(1:size(mask,1),1:size(mask,2),1:size(stack,3)); % remove the mirrored input
    
    csd_complex = csd;
    csd_sum = sum(csd_complex, 3);
    QP = imag(csd_sum);
    % QP = angle(csd_sum + mean(stack(:))/s.optics.alpha);
    % Amplitude = abs(csd);
    % QP = angle(csd + mean(stack(:))/s.optics.alpha);
    close all



end