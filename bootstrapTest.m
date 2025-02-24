% 1. Laad je dataset in (40 datapunten)
data = [
1.358600589
0.926639296
1.116943829
0.952714456

1.24853433
1.244734353

0.938617336

1.01673387

1.22789444
0.97889192
0.997217835
1.066755146
1.134436628
1.167496775

1.092784894
1.063265867
0.981834555

0.999619675
1.115391658
0.948068484
1.082130371
1.029539348

0.967239698


1.149330757


1.167608814
1.263880745
1.038442828
0.886212509


0.981884546


1.256692719
1.095845509
1.164104353
1.083303925

];

data = abs(data);

% 2. Aantal bootstrap-simulaties
B = 1000; % Hoe meer, hoe nauwkeuriger

% 3. Bootstrap het gemiddelde
boot_means = bootstrp(B, @mean, data);

% 4. Bereken 95% betrouwbaarheidsinterval
CI = prctile(boot_means, [2.5, 97.5]); 

% 5. Resultaten weergeven
disp(['Steekproefgemiddelde: ', num2str(mean(data))]);
disp(['95% betrouwbaarheidsinterval: ', num2str(CI(1)), ' tot ', num2str(CI(2))]);
figure()
% 6. Visualisatie (optioneel, maar handig)
histogram(boot_means, 30);
hold on;
xline(CI(1), 'r', 'LineWidth', 2); % Ondergrens CI
xline(CI(2), 'r', 'LineWidth', 2); % Bovengrens CI
xline(mean(data), 'g', 'LineWidth', 2); % Steekproefgemiddelde
title('Bootstrap distribution of the mean');
xlabel('Mean');
ylabel('Frequence');
legend('Bootstrap-averages', '95% CI - low', '95% CI - high', 'sample average');
hold off;
mu0 = mean(data)
Relatieve_breedte = (CI(2) - CI(1))./mu0;



