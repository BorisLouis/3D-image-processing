ResultsDiff = readtable('G:\multicolor_polarization\multicolor\20240813_oil_MC_NPs_200nm\Results_MC_diff.xlsx');

ResultsCh1 = table2array(ResultsDiff(1:end,["Var2", "Var3", "Var4", "Var5", "Var6"]));
ResultsCh2 = table2array(ResultsDiff(1:end,["Var9", "Var10", "Var11", "Var12", "Var13"]));

Results = [ResultsCh1(:), ResultsCh2(:)];


figure()
boxplot(Results)

mean1 = nanmean(ResultsCh1(:));
mean2 = nanmean(ResultsCh2(:));

std1 = std(ResultsCh1(:), "omitmissing");
std2 = std(ResultsCh2(:), "omitmissing");

error = abs((mean1-mean2))./((mean1+mean2)./2);

