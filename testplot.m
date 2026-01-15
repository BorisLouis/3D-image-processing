close all

FieldNames = fieldnames(BigResults_500);
Planes = 1:8;
for i = 1:size(FieldNames, 1)
    FieldName = FieldNames{i};
    List500 = nanmean(rmoutliers(BigResults_500.(FieldName),1));
    List1000 = nanmean(rmoutliers(BigResults_1000.(FieldName),1));
    List2000 = nanmean(rmoutliers(BigResults_2000.(FieldName),1));

    Fig = figure()
    hold on
    scatter(Planes, List500, 'filled');
    scatter(Planes, List1000, 'filled');
    scatter(Planes, List2000, 'filled');
    legend({'500 nm PS', '1000 nm PS', '2000 nm PS'})
    xlabel('Plane')
    ylabel('Phase (°)')
    title(FieldName)   

    if i == 1
        ylim([0 15])
    elseif i == 2
        ylim([0 80])
    elseif i == 3
        ylim([0 120])
    elseif i == 4
        ylim([-1.5 1.5])
    end
    saveas(Fig, append('D:\Data Hannah\202512117\Brightfield', filesep, FieldName, '.png'))
end