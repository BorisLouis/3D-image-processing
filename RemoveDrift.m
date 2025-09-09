% Drift correction script
% Root folder
rootDir = 'D:\Multimodal tracking\20250724\alldata';

% Get list of subfolders
subFolders = dir(rootDir);
subFolders = subFolders([subFolders.isdir]); % keep only directories
subFolders = subFolders(~ismember({subFolders.name},{'.','..'})); % remove . and ..

for i = 1:numel(subFolders)
    try
        folderPath = fullfile(rootDir, subFolders(i).name);
        matFile = fullfile(folderPath, 'TraceswPhase.mat');
    
        if ~isfile(matFile)
            fprintf('No TraceswPhase.mat in %s, skipping...\n', folderPath);
            continue;
        end
        
        % Load data
        S = load(matFile, 'Traces3D');
        Traces3D = S.Traces3D;
    
        % Collect all displacements for drift estimation
        allDisp = []; % store frame-to-frame displacements
    
        for t = 1:numel(Traces3D)
            traj = Traces3D{t};
            if size(traj,1) > 1
                dx = diff(traj.row);
                dy = diff(traj.col);
                dz = diff(traj.z);
                allDisp = [allDisp; [dx, dy, dz]];
            end
        end
    
        % Estimate drift per frame as the mean displacement across all trajectories
        driftPerFrame = mean(allDisp,1,'omitnan');  % [drift_x drift_y drift_z]
    
        fprintf('Folder: %s, Estimated drift/frame = (%.4f, %.4f, %.4f)\n', ...
                subFolders(i).name, driftPerFrame(1), driftPerFrame(2), driftPerFrame(3));
    
        % Apply correction
        for t = 1:numel(Traces3D)
            traj = Traces3D{t};
            nFrames = size(traj,1);
            correction = (0:nFrames-1)' * driftPerFrame; % linear correction
            traj.row = traj.row - correction(:,1);
            traj.col = traj.col - correction(:,2);
            traj.z   = traj.z   - correction(:,3);
            Traces3D{t} = traj;
        end
    
        % Save back into the same file
        save(fullfile(folderPath, 'TraceswPhaseCorr.mat'), 'Traces3D');
        fprintf('Corrected drift and saved to %s\n', matFile);
    catch
    end
end

disp('Drift correction completed for all subfolders.');