MainFolder = 'F:\Data Hannah\202512117\Darkfield';
SubFolders = {'500nmPS', '1000nmPS', '2000nmPS'};
for k = 1:3
    for i = 1:5
        try
            if strcmp(SubFolders{k}, '500nmPS')
                SubsubFolder = '500nmPS_Darkfield_sample__';
            elseif strcmp(SubFolders{k}, '1000nmPS')
                SubsubFolder = '1000nmPS_Darkfield_sample__';
            elseif strcmp(SubFolders{k}, '2000nmPS')
                SubsubFolder = '2000nmPS_Darkfield_sample__';
            end
        
            Name = append(MainFolder, filesep, SubFolders{k}, filesep, SubsubFolder, num2str(i));
            for j = 1:8
                try
                    Name2 = append(Name, filesep, 'PhaseIm_minProjection_plane', num2str(j), '.mat');
                    load(Name2);
                    MinPhase(i,j) = max(MeanImage, [], 'all');
                catch
                    disp('idk');
                end
            end
        catch
            disp('idk');
        end
    end
    Results.(append('Phase_', SubFolders{k})) = MinPhase;
    MinPhase = [];
end