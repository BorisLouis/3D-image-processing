function plotCal(ch1,ch2,ROI)

if size(ROI,1) <=4
    it1 = size(ROI,1);
    dim = 1;

    figure
    for i= 1:it1
        subplot(dim,it1,i)
        imagesc(ch1(ROI(i,2):ROI(i,2)+ROI(i,4)-1,ROI(i,1):ROI(i,1)+ROI(i,3)-1));
    end

    
else
    it1 = 4;
    dim = 2;
    
    figure
    Channel1 = ch1(ROI(1,2):ROI(1,2)+ROI(1,4)-1,ROI(1,1):ROI(1,1)+ROI(1,3)-1);
    for i= 1:it1
        subplot(dim,it1,i)
        Channelx = ch1(ROI(i,2):ROI(i,2)+ROI(i,4),ROI(i,1):ROI(i,1)+ROI(i,3));
        Channelx = imresize(Channelx, size(Channel1));
        imshowpair(Channel1, Channelx)
        % imagesc(ch1(ROI(i,2):ROI(i,2)+ROI(i,4)-1,ROI(i,1):ROI(i,1)+ROI(i,3)-1));
        title(append('Plane 1 x Plane ', num2str(i)))
    end

    for i= 1:it1
        subplot(dim,it1,it1+i)
        Channelx = ch2(ROI(4+i,2):ROI(4+i,2)+ROI(4+i,4),ROI(4+i,1):ROI(4+i,1)+ROI(4+i,3));
        Channelx = imresize(Channelx, size(Channel1));
        imshowpair(Channel1, Channelx)
        title(append('Plane 1 x Plane ', num2str(i+4)))
    end

    
end

end