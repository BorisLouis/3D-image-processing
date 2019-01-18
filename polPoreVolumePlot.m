%% ATTENTION !!!! NO CLEARING OF DATA HERE AS WE USE DATA GENERATED BY 
%FULLLOADJOHANNES.M

nStack = 5;
E1K = [exist('S1_1K','var'), exist('S2_1K','var')];
n1K = sum(E1K);
E5K = [exist('S1_5K','var'), exist('S2_5K','var')];
n5K = sum(E5K);
n2Plot = sum(E1K +E5K); %count number of variable present


Concentrations = [0.25,0.5,1];
%% Volume picture 1K
szMarker = 20;
assert(~all(E1K==0),'No Data found for 1K')

polRat025 = zeros(nStack*2,1);
polRat050 = zeros(nStack*2,1);
polRat100 = zeros(nStack*2,1);
poreRat025 = zeros(nStack*2,1);
poreRat050 = zeros(nStack*2,1);
poreRat100 = zeros(nStack*2,1);

for j = 1:nStack
    
    polRat025(1+2*(j-1)) = S1_1K.Table025(j).Data(1,:).ratioPol;
    polRat025(2+2*(j-1)) = S2_1K.Table025(j).Data(1,:).ratioPol;
    
    polRat050(1+2*(j-1)) = S1_1K.Table050(j).Data(1,:).ratioPol;
    polRat050(2+2*(j-1)) = S2_1K.Table050(j).Data(1,:).ratioPol;
    
    polRat100(1+2*(j-1)) = S1_1K.Table100(j).Data(1,:).ratioPol;
    polRat100(2+2*(j-1)) = S2_1K.Table100(j).Data(1,:).ratioPol;
    
    poreRat025(1+2*(j-1)) = S1_1K.Table025(j).Data(1,:).ratioPores;
    poreRat025(2+2*(j-1)) = S2_1K.Table025(j).Data(1,:).ratioPores;
    
    poreRat050(1+2*(j-1)) = S1_1K.Table050(j).Data(1,:).ratioPores;
    poreRat050(2+2*(j-1)) = S2_1K.Table050(j).Data(1,:).ratioPores;
    
    poreRat100(1+2*(j-1)) = S1_1K.Table100(j).Data(1,:).ratioPores;
    poreRat100(2+2*(j-1)) = S2_1K.Table100(j).Data(1,:).ratioPores;
       
end

polRat = [mean(polRat025),mean(polRat050),mean(polRat100)];
poreRat = [mean(poreRat025),mean(poreRat050),mean(poreRat100)];

errPol = [3*std(polRat025),3*std(polRat050),3*std(polRat100)];
errPore = [3*std(poreRat025),3*std(poreRat050),3*std(poreRat100)];

figure()
hold on

errorbar(Concentrations,polRat,errPol,'bo');
errorbar(Concentrations,poreRat,errPol,'ro');

xlim([0 1]);
ylim([0 1]);
xlabel('PIC concentration(mg/mL)');
ylabel('volume ratio');

legend({'polymer volume ratio','pore volume ratio'});



%% Volume picture 5 K

szMarker = 20;
assert(~all(E1K==0),'No Data found for 1K')

polRat025 = zeros(nStack*2,1);
polRat050 = zeros(nStack*2,1);
polRat100 = zeros(nStack*2,1);
poreRat025 = zeros(nStack*2,1);
poreRat050 = zeros(nStack*2,1);
poreRat100 = zeros(nStack*2,1);

for j = 1:nStack
    
    polRat025(1+2*(j-1)) = S1_5K.Table025(j).Data(1,:).ratioPol;
    polRat025(2+2*(j-1)) = S2_5K.Table025(j).Data(1,:).ratioPol;
    
    polRat050(1+2*(j-1)) = S1_5K.Table050(j).Data(1,:).ratioPol;
    polRat050(2+2*(j-1)) = S2_5K.Table050(j).Data(1,:).ratioPol;
    
    polRat100(1+2*(j-1)) = S1_5K.Table100(j).Data(1,:).ratioPol;
    polRat100(2+2*(j-1)) = S2_5K.Table100(j).Data(1,:).ratioPol;
    
    poreRat025(1+2*(j-1)) = S1_5K.Table025(j).Data(1,:).ratioPores;
    poreRat025(2+2*(j-1)) = S2_5K.Table025(j).Data(1,:).ratioPores;
    
    poreRat050(1+2*(j-1)) = S1_5K.Table050(j).Data(1,:).ratioPores;
    poreRat050(2+2*(j-1)) = S2_5K.Table050(j).Data(1,:).ratioPores;
    
    poreRat100(1+2*(j-1)) = S1_5K.Table100(j).Data(1,:).ratioPores;
    poreRat100(2+2*(j-1)) = S2_5K.Table100(j).Data(1,:).ratioPores;
       
end

polRat = [mean(polRat025),mean(polRat050),mean(polRat100)];
poreRat = [mean(poreRat025),mean(poreRat050),mean(poreRat100)];

errPol = [3*std(polRat025),3*std(polRat050),3*std(polRat100)];
errPore = [3*std(poreRat025),3*std(poreRat050),3*std(poreRat100)];

figure()
hold on

errorbar(Concentrations,polRat,errPol,'bo');
errorbar(Concentrations,poreRat,errPol,'ro');

xlim([0 1]);
ylim([0 1]);
xlabel('PIC concentration(mg/mL)');
ylabel('volume ratio');

legend({'polymer volume ratio','pore volume ratio'});