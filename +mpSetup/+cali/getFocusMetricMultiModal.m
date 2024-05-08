function [ focus_met, in_focus, fit ] = getFocusMetric( chData1c, chData2c, chData3c, chData4c, Z1, Z2)
%GETFOCUSMETRIC gets information about when the channels are in focus.
%Inspired in the work done in EPFL for 3D SOFI

    N = size(chData1c,4);
    %nChan = size(chData1c,3);
    nChan = 4;
    focus_met = zeros(N, nChan*2);
    fit = nan(N, nChan*2*2);
    in_focus(nChan*4).cam   =  [];
    in_focus(nChan*4).ch    =  [];
    in_focus(nChan*4).frame =  [];
    in_focus(nChan*4).zpos  =  [];
    
    for i = 1:nChan
        in_focus(i).cam   =  1;
        in_focus(i).ch    =  i;
        % get image sequence for 2 channels to compare
        %remove border pixel due to bug with one data set, should not
        %affect others
        tmp=squeeze(chData1c(10:end-10,10:end-10,i,:));
        
        % see when each channel is in focus
        %focus_met(:,i) = squeeze(mean(max(tmp)));
        focus_met(:,i) = squeeze(mean(max(imgradient3(tmp))));
        [ zFocus, fitTmp ] = mpSetup.cali.getSubResPlanePosition(focus_met(:,i),Z1);
        fit(1:length(fitTmp),2*i-1:2*i,:) = fitTmp;
        in_focus(i).zpos = zFocus;
        [~,in_focus(i).frame] = min(abs(Z1-zFocus));
    end
    
    for i = nChan+1:nChan*2
        in_focus(i).cam   =  2;
        in_focus(i).ch    =  i-nChan;
        % get image sequence for 2 channels to compare
        tmp=squeeze(chData2c(10:end-10,10:end-10,i-nChan,:));
        
        % see when each channel is in focus
        %focus_met(:,i) = squeeze(mean(max(tmp)));
        focus_met(:,i) = squeeze(mean(max(imgradient3(tmp))));
        
        [ zFocus, fitTmp ] = mpSetup.cali.getSubResPlanePosition(focus_met(:,i),Z2);
        fit(1:length(fitTmp),2*i-1:2*i,:) = fitTmp;
        in_focus(i).zpos = zFocus;
        [~,in_focus(i).frame] = min(abs(Z2-zFocus));
    end    

    for i = nChan*2+1:nChan*3
        in_focus(i).cam   =  3;
        in_focus(i).ch    =  i-nChan*2;
        % get image sequence for 2 channels to compare
        tmp=squeeze(chData3c(10:end-10,10:end-10,i-(nChan*2),:));
        
        % see when each channel is in focus
        %focus_met(:,i) = squeeze(mean(max(tmp)));
        focus_met(:,i) = squeeze(mean(max(imgradient3(tmp))));
        
        [ zFocus, fitTmp ] = mpSetup.cali.getSubResPlanePosition(focus_met(:,i),Z1);
        fit(1:length(fitTmp),2*i-1:2*i,:) = fitTmp;
        in_focus(i).zpos = zFocus;
        [~,in_focus(i).frame] = min(abs(Z1-zFocus));
    end

    for i = nChan*3+1:nChan*4
        in_focus(i).cam   =  4;
        in_focus(i).ch    =  i-nChan*3;
        % get image sequence for 2 channels to compare
        tmp=squeeze(chData4c(10:end-10,10:end-10,i-(nChan*3),:));
        
        % see when each channel is in focus
        %focus_met(:,i) = squeeze(mean(max(tmp)));
        focus_met(:,i) = squeeze(mean(max(imgradient3(tmp))));
        
        [ zFocus, fitTmp ] = mpSetup.cali.getSubResPlanePosition(focus_met(:,i),Z2);
        fit(1:length(fitTmp),2*i-1:2*i,:) = fitTmp;
        in_focus(i).zpos = zFocus;
        [~,in_focus(i).frame] = min(abs(Z2-zFocus));
    end
end