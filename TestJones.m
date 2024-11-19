theta_deg = 0:1:360; % HWP angles in degrees
theta = deg2rad(theta_deg); % Convert to radians

i = [0, 10, 20, 30, 40, 50];
% Ein = [1; 0]; %for 0° polarized light
% Ein = (1/sqrt(2))*[1; 1]; % for 45° polarized light
% Ein = [0; 1]; %for 90° polarized light

figure();
for j = 1:size(i,2)
    Ein = [cos(i(j)); sin(i(j))];
    
    Ip = zeros(size(theta)); % Intensity in p-polarized channel
    Is = zeros(size(theta)); % Intensity in s-polarized channel
    
    for k = 1:size(theta,2)
        HWP = [cos(2*theta(k)), sin(2*theta(k)); 
               sin(2*theta(k)), -cos(2*theta(k))];
        
        Eout = HWP * Ein;
        
        Ep = [1, 0; 0, 0] * Eout; % p-polarized
        Es = [0, 0; 0, 1] * Eout; % s-polarized
        
        Ip(k) = norm(Ep)^2;
        Is(k) = norm(Es)^2;
    end
    
    subplot(size(i,2), 1, j)
    plot(theta_deg, Ip, 'r', 'LineWidth', 2); 
    hold on;
    plot(theta_deg, Is, 'b', 'LineWidth', 2);
    xlabel('theta (degrees)');
    ylabel('Intensity');
    legend('channel1', 'channel2');
    text = append('Jones Formalism - input angle polarized light ', num2str(i(j)), '°');
    title(text);
    hold on
end