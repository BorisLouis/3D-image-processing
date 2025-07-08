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

            mkdir(append(obj.raw.movInfo.Path, filesep, 'SegmentMovie'));

            n = 1;
            Step = 0;
            ChunkSize = 100;
            for k = 1:ChunkSize:nFrames
                Step = Step + 1;
                idx = k:min(k+ChunkSize-1, nFrames);
                Startidx = idx(1)-1;
                for i = idx
                    waitbar(n./nFrames, f, append('Segmenting frame ', num2str(n), ' out of ', num2str(nFrames), '...'))
                    CurrFrameAll = obj.getFrame(n, q);
                    for j = 1:size(CurrFrameAll, 3)
                        CurrFrame = CurrFrameAll(:,:,j);
                        CurrFrameBg = CurrFrame - imgaussfilt(CurrFrame, obj.info.GlobalBgCorr);
                        CurrFrameBg(CurrFrameBg < 0) = 0;
                        CurrFrameBg = mat2gray(CurrFrameBg);
                        CurrFrameEnhanced = imadjust(CurrFrameBg);
                        [~,mask(:,:,i-Startidx, j)] = imSegmentation.segmentStack(CurrFrameEnhanced, 'method', 'adaptive');
                    end
        
                    if strcmp(obj.info.ShowSegmentation, 'on')
                        if i == obj.info.TestFrame
                            if strcmp(obj.info.Dimension, '3D')
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
                                imshowpair(CurrFrame, mask(:,:,i, 4))
                                title('Overlay')
                                axis image
                                sgtitle(append('TestFrame ', num2str(i), ' - plane 4'));
                            elseif strcmp(obj.info.Dimension, '2D')
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
                            end
    
                            Filename = append(obj.raw.movInfo.Path, filesep, 'Segmentation_Testframe_', num2str(i), '.png');
                            saveas(Fig, Filename);
                        end
                    end    
                    n = n+1;
                end
                Filename = append(obj.raw.movInfo.Path, filesep, 'SegmentMovie', filesep, 'SegmentMovie', num2str(Step), '.mat');
                save(Filename, 'mask');
                mask = [];
            end
            close(f)
        end

        function SaveMask(obj, q)
            FileName = append(obj.calibrated{1, 1}.mainPath, filesep, 'SegmentMask');
            Mask = obj.SegmentMap;
            save(FileName, "Mask");
            disp(append("== Segmentmask saved - Channel ", num2str(q), " =="))
        end
    end
end

