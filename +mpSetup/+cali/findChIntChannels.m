function [Icorrf2] = findChIntChannels(Icorrf1, Icorrf2, IntCh1, IntCh2)
    %Function to put the intensity of channel 2 the same as channel 1
    
    for i = 1:size(IntCh1)
        IntCh1(i,5) = IntCh1(i,5)*Icorrf1(IntCh1(i,1));
        IntCh1(i,6) = IntCh1(i,6)*Icorrf1(IntCh1(i,2));

        IntCh2(i,5) = IntCh2(i,5)*Icorrf2(IntCh2(i,1));
        IntCh2(i,6) = IntCh2(i,6)*Icorrf2(IntCh2(i,2));
    end

    Ratio = nanmean(IntCh1(:,5:6)./IntCh2(:,5:6), 'all');

    Icorrf2 = Icorrf2*Ratio;
end

