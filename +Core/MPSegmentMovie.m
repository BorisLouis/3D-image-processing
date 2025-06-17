classdef MPSegmentMovie < Core.MPMovie

    properties
        SegmentMap
    end
    
    methods
        function obj = MPSegmentMovie(raw,cal,info)
            
            obj  = obj@Core.MPMovie(raw,cal,info);
            
        end
        
        function getSegmentMovie(obj, q)
            f = waitbar(0, 'Initializing segmentation algorithm...');
            if strcmp(obj.info.frame2Load, 'all')
                nFrames = obj.calibrated{1, 1}.nFrames; 
            elseif isa(obj.info.frame2Load, 'double')
                nFrames = max(obj.info.frame2Load);
            end

            for i = 1:nFrames
                waitbar(i./nFrames, f, append('Running segmentation: frame ', num2str(i), ' out of ', num2str(nFrames)));
                CurrFrame = obj.getFrame(i,q);
                imgaussfilt(CurrFrame, 20);
                CurrFrameBg = CurrFrame - imgaussfilt(CurrFrame, 25);
                CurrFrameBg(CurrFrameBg < 2000) = 0;

                CurrFrameBg = imbinarize(CurrFrameBg);
                BurrFrameBg = bwareaopen(CurrFrameBg, 15, 8);

                se = strel('disk', 2);
                Mask(:,:,i) = imopen(CurrFrameBg, se);

                if obj.info.ShowSegment == 1
                    if i == 10
                        figure()
                        subplot(1,2,1)
                        imagesc(CurrFrame)
                        axis image
                        subplot(1,2,2)
                        imshowpair(CurrFrame, CurrFrameBg)
                        axis image
                    end
                end
            end
            close(f)
            obj.SegmentMap = Mask;
        end
    end
end

