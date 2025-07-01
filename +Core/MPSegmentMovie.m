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
                waitbar(i./nFrames, f, append('Segmenting frame ', num2str(i), ' out of ', num2str(nFrames), '...'))
                CurrFrame = obj.getFrame(i, q);
                CurrFrameBg = CurrFrame - imgaussfilt(CurrFrame, obj.info.GlobalBgCorr);
                CurrFrameBg(CurrFrameBg < 0) = 0;
                CurrFrameBg = mat2gray(CurrFrameBg);

                CurrFrameEnhanced = imadjust(CurrFrameBg);
                [~,mask(:,:,i)] = imSegmentation.segmentStack(CurrFrameEnhanced, 'method', 'adaptive');

                if strcmp(obj.info.ShowSegmentation, 'on')
                    if i == obj.info.TestFrame
                        Fig = figure();
                        subplot(1,3,1)
                        imagesc(CurrFrame)
                        title('Raw data')
                        axis image
                        subplot(1,3,2)
                        imagesc(mask(:,:,i))
                        title('Segment mask')
                        axis image
                        subplot(1,3,3)
                        imshowpair(CurrFrame, mask(:,:,i))
                        title('Overlay')
                        axis image
                        sgtitle(append('TestFrame ', num2str(i)));

                        Filename = append(obj.raw.movInfo.Path, filesep, 'Segmentation_Testframe_', num2str(i), '.png');
                        saveas(Fig, Filename);
                    end
                end       
            end
            close(f)
            obj.SegmentMap = mask;
        end

        function SaveMask(obj, q)
            FileName = append(obj.calibrated{1, 1}.mainPath, filesep, 'SegmentMask');
            Mask = obj.SegmentMap;
            save(FileName, "Mask");
            disp(append("== Segmentmask saved - Channel ", num2str(q), " =="))
        end
    end
end

