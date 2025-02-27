for i = 1:size(CommonTraces, 1)
    Trace1 = CommonTraces.Int1{i,1};
    Trace2 = CommonTraces.Int2{i,1};

    TotTrace = Trace1 + Trace2;
    MeanTot = median(TotTrace);
    
    CommonTraces.TotIntTrace{i,1} = TotTrace;
    CommonTraces.MeanTotInt{i,1} = MeanTot;
end