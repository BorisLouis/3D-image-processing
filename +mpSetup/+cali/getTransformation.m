function [transformation, ZDiff] = getTransformation(obj)
    for i = 1:obj.cal.nFiles %% loop through files
        fmov = append('MPCal', num2str(i));
        [frameInfo, movInfo, ~ ]= Load.Movie.ome.getInfo(obj.MPCalibrations.(fmov).raw.fullPath);
        [ movC1, movC2, idx ] = Load.Movie.ome.load( frameInfo, movInfo  , 1:movInfo.maxFrame );
        if obj.MPCalibrations.(fmov).cal.file.flipCam2
            movC2 = flip(movC2,2);
        end

        [ chData1c, chData2c ] = mpSetup.cali.getChData( movC1, movC2, obj.cal.file.ROI1 );
        [ chData3c, chData4c ] = mpSetup.cali.getChData( movC1, movC2, obj.cal.file.ROI2FullCam);
        f = waitbar(0, append('Calculating transformation: Sample ', num2str(i), '/', num2str(obj.cal.nFiles)));
        for j = 1:size(obj.cal.file.inFocus1 ,2)
            waitbar(j/8, f, append('Calculating transformation: Sample ', num2str(i), '/', num2str(obj.cal.nFiles)));
            %%Get images of plane in focus 
            Movie = append('MPCal', num2str(i));
            Ch1Frame = obj.MPCalibrations.(Movie).cal.file.inFocus1(j).frame;
            Ch2Frame = obj.MPCalibrations.(Movie).cal.file.inFocus2(j).frame;
            Ch1ZPos = obj.MPCalibrations.(Movie).cal.file.inFocus1(j).zpos;
            Ch2ZPos = obj.MPCalibrations.(Movie).cal.file.inFocus2(j).zpos;
            ZDiffFocus = Ch2ZPos - Ch1ZPos;

            if j < 5
                Channel1 = chData1c(:,:,j, Ch1Frame);
                Channel2 = chData3c(:,:,j, Ch1Frame);
            else
                Channel1 = chData2c(:,:,j-4, Ch1Frame);
                Channel2 = chData4c(:,:,j-4, Ch1Frame);
            end
            
            se = strel('disk', 10);
            bgChannel2 = imopen(Channel2, se);
            bgChannel1 = imopen(Channel1, se);
            Channel2 = Channel2 - bgChannel2;
            Channel1 = Channel1 - bgChannel1;         
            similarity = 0;
            n = 0;

            % while similarity < 0.99985
                Channel2 = Channel2(1+n:end-1, 1+n:end-n);
                figure()
                subplot(1,2,1)
                imshowpair(Channel1,Channel2)
                title("before correction");
                hold on
    
                config = "monomodal";
                transf = "similarity";
                    
                [optimizer,metric] = imregconfig(config);
    
                tform = imregcorr(Channel2,Channel1,transf);
                movingRegistered = imwarp(Channel2,tform,"OutputView",imref2d(size(Channel1)));
                similarity = multissim(movingRegistered, Channel1); 
    
                if similarity > 0.99985
                    if n == 0
                        Scale(i,j) = tform.Scale;
                    else
                        Scale(i,j) = 1;
                    end
                    RotationAngle(i,j) = tform.RotationAngle;
                    Translation1(i,j) = tform.Translation(1);
                    Translation2(i,j) = tform.Translation(2);
                    ZDiff(i,j) = ZDiffFocus;
                else
                    Scale(i,j) = NaN;
                    RotationAngle(i,j) = NaN;
                    Translation1(i,j) = NaN;
                    Translation2(i,j) = NaN;
                    ZDiff(i,j) = NaN;
                end
    
                subplot(1,2,2)
                imshowpair(Channel1,movingRegistered);
                title("after correction");
                sgtitle(append("Plane ", num2str(j), " x Plane ", num2str(j+8)));
                n = n+1;
            %     if n > 3
            %         break
            %     end
            % end
        end
        close all
        close(f)
    end
    for i = 1:size(obj.cal.file.inFocus1 ,2)
        AngleList = deg2rad(RotationAngle(:,i));
        mean_x = mean(cos(AngleList));
        mean_y = mean(sin(AngleList));
        mean_angle = rad2deg(atan2(mean_y, mean_x));
        tform = simtform2d(nanmedian(Scale(:,i)), mean_angle, [nanmedian(Translation1(:,i)), nanmedian(Translation2(:,i))]);
        transformation{i,1} = tform;
    end
end

