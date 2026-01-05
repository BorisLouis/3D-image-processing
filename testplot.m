close all
rows = [70, 345];
cols = [144, 159];
Fig = figure();
for i = 1:2

Profile = squeeze(QPmap(rows(i), cols(i), 4, :));
plot(Profile)
LegendNames{i} = append('Particle ', num2str(i));
hold on
end
ylabel('Phase (°)')
xlabel('Frame z-stack')
title('1000nm brightfield sample1 plane 4')

legend(LegendNames)
saveas(Fig, append('D:\Data Hannah\202512117\Brightfield', filesep, '1000nm_s1.png'));
