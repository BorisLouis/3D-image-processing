ResultsDiff = readtable('G:\multicolor_polarization\multicolor\20240906_polymerisation_dual_color\300nm_60nm_sample_1_HETEROGENEITY\Results.xlsx');

Times = {'3' '5' '7' '9' '11' '13'};
TimeColors = {[0.4980392156862745 0.6901960784313725 0.4117647058823529], [0.5725490196078431 0.11372549019607843 0.08235294117647059]};
TrackRes1 = table2array(ResultsDiff(2:92,["Var2", "Var3", "Var4", "Var5", "Var6", "Var7"]));
TrackRes2 = table2array(ResultsDiff(2:92,["Var10","Var11", "Var12", "Var13", "Var14", "Var15"]));
meanTrackRes1 = nanmean(TrackRes1, 1);
meanTrackRes2 = nanmean(TrackRes2, 1);

GroupedData = {TrackRes1 TrackRes2}; 

legendEntries = {'200 nm PS' '60 nm AuNP'};

N = numel(GroupedData);
delta = linspace(-.1,.1,N); %// define offsets to distinguish plots
width = .2; %// small width to avoid overlap
cmap = hsv(N); %// colormap
legWidth = 1.8; %// make room for legend

figure;
hold on;

for ii=1:N 
    labels = Times; 
    boxplot(GroupedData{ii},'Color', TimeColors{ii}, 'boxstyle','filled', ...
        'position',(1:numel(labels))+delta(ii), 'labels',labels, 'width', 0.1, 'outliersize', 0.000001)
    plot(NaN,1,'color',TimeColors{ii});
end
xlabel('Polymerization time (min)'); ylabel('Diffusion coefficient (µm^2s^-^1)');
xlim([1+2*delta(1) numel(labels)+2*delta(N)])
ylim([-0.5 7])

MeanTable = [nanmean(GroupedData{1,1}, 1); nanmean(GroupedData{1,2}, 1)];
for ii = 1:size(MeanTable, 2)
    Difference (1, ii) = (abs(MeanTable(1,ii)-MeanTable(2,ii)))./mean([MeanTable(1,ii); MeanTable(2,ii)])*100;
end

legend(legendEntries);