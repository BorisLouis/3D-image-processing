clear; clc; close all;

%% ------------------ Parameters ------------------
nBlack = 500;       
nGray  = 80;        
nGreen = 40;        
nFrames = 120;      
boxSize = 7;        
dt = 0.9;           
collisionThresh = 0.15; 
spacing = 0.15;      
maxClusterSize = 4;  

%% ------------------ Initialize positions and velocities ------------------
blackPos = boxSize*rand(nBlack,2);
grayPos  = boxSize*rand(nGray,2);
greenPos = boxSize*rand(nGreen,2);

blackVel = 0.1*randn(nBlack,2);
grayVel  = 0.05*randn(nGray,2);
greenVel = 0.08*randn(nGreen,2);

%% ------------------ Initialize clusters ------------------
blackClusters = struct('members', num2cell((1:nBlack)'), ...
                       'relPos', num2cell(zeros(1,2),2), ...
                       'center', num2cell(blackPos,2), ...
                       'vel', num2cell(blackVel,2), ...
                       'angle', num2cell(zeros(nBlack,1)), ...
                       'rotVel', num2cell(0.2*randn(nBlack,1)));

grayClusters = struct('members', num2cell(nBlack+(1:nGray)'), ...
                      'relPos', num2cell(zeros(1,2),2), ...
                      'center', num2cell(grayPos,2), ...
                      'vel', num2cell(grayVel,2), ...
                      'angle', num2cell(zeros(nGray,1)), ...
                      'rotVel', num2cell(0.1*randn(nGray,1)));

%% ------------------ Video 1: Initial cluster formation ------------------
figure('Color','w'); axis equal off;
filename1 = 'initial_clusters.gif';

for frame = 1:nFrames
    cla; hold on;
    allClusters = [blackClusters; grayClusters];

    %% Move clusters
    for c = 1:length(allClusters)
        clusterSize = length(allClusters(c).members);
        speedFactor = 1/sqrt(clusterSize); 

        allClusters(c).center = reshape(allClusters(c).center,1,2);
        allClusters(c).vel    = reshape(allClusters(c).vel,1,2);

        allClusters(c).center = allClusters(c).center + allClusters(c).vel*dt*speedFactor;
        allClusters(c).angle  = allClusters(c).angle  + allClusters(c).rotVel*dt*speedFactor;

        pos = allClusters(c).center; vel = allClusters(c).vel;
        if pos(1)<=0 || pos(1)>=boxSize; vel(1)=-vel(1); end
        if pos(2)<=0 || pos(2)>=boxSize; vel(2)=-vel(2); end
        allClusters(c).vel = vel;
    end

    %% Move green dots
    for g = 1:nGreen
        greenPos(g,:) = greenPos(g,:) + greenVel(g,:)*dt;
        if greenPos(g,1)<=0 || greenPos(g,1)>=boxSize; greenVel(g,1)=-greenVel(g,1); end
        if greenPos(g,2)<=0 || greenPos(g,2)>=boxSize; greenVel(g,2)=-greenVel(g,2); end
    end

    %% Cluster collisions
    i = 1;
    while i < length(allClusters)
        j = i+1;
        while j <= length(allClusters)
            ci = reshape(allClusters(i).center,1,2);
            cj = reshape(allClusters(j).center,1,2);
            distVec = ci - cj;

            if norm(distVec) < collisionThresh
                totalMembers = numel(allClusters(i).members) + numel(allClusters(j).members);
                if totalMembers <= maxClusterSize
                    allMembers = [allClusters(i).members, allClusters(j).members];
                    newRelPos = generateRelPosPolymer(allMembers, spacing);
                    allClusters(i).members = allMembers;
                    allClusters(i).relPos  = newRelPos;
                    allClusters(i).vel     = mean([reshape(allClusters(i).vel,1,2); reshape(allClusters(j).vel,1,2)],1);
                    allClusters(i).center  = mean([ci; cj],1);
                    allClusters(j) = [];
                else
                    v1 = reshape(allClusters(i).vel,1,2);
                    v2 = reshape(allClusters(j).vel,1,2);
                    n = distVec/norm(distVec);
                    allClusters(i).vel = v1 - 2*dot(v1,n)*n;
                    allClusters(j).vel = v2 - 2*dot(v2,n)*n;
                    j = j + 1;
                end
            else
                j = j + 1;
            end
        end
        i = i + 1;
    end

    %% Green dots bounce off clusters
    for g = 1:nGreen
        for c = 1:length(allClusters)
            clusterRel = sanitizeRelPos(allClusters(c).relPos);
            ctr = reshape(allClusters(c).center,1,2);
            clusterPos = clusterRel + repmat(ctr, size(clusterRel,1), 1);

            for m = 1:size(clusterPos,1)
                diffVec = greenPos(g,:) - clusterPos(m,:);
                dist = norm(diffVec);
                if dist < collisionThresh && dist>0
                    normal = diffVec / dist;
                    greenVel(g,:) = greenVel(g,:) - 2*dot(greenVel(g,:), normal)*normal;
                end
            end
        end
    end

    %% Split back
    blackClusters = allClusters(cellfun(@(m) all(m<=nBlack), {allClusters.members}));
    grayClusters  = allClusters(~cellfun(@(m) all(m<=nBlack), {allClusters.members}));

    %% Plot clusters
    for c = 1:length(allClusters)
        theta = allClusters(c).angle;
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        relPos = sanitizeRelPos(allClusters(c).relPos);
        ctr = reshape(allClusters(c).center,1,2);
        clusterPos = (R*relPos')' + repmat(ctr, size(relPos,1), 1);

        blackIdx = allClusters(c).members <= nBlack;
        grayIdx  = allClusters(c).members > nBlack;
        if any(blackIdx)
            plot(clusterPos(blackIdx,1), clusterPos(blackIdx,2), 'ko','MarkerFaceColor','k','MarkerSize',6);
        end
        if any(grayIdx)
            plot(clusterPos(grayIdx,1), clusterPos(grayIdx,2), 'ko','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',8);
        end
    end

    plot(greenPos(:,1), greenPos(:,2), 'go','MarkerFaceColor','g','MarkerSize',6);
    xlim([0 boxSize]); ylim([0 boxSize]);
    drawnow;

    frameImg = getframe(gcf);
    [A,map] = rgb2ind(frame2im(frameImg),256);
    if frame == 1
        imwrite(A,map,filename1,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename1,'gif','WriteMode','append','DelayTime',0.1);
    end
end

%% ------------------ VIDEO 2: Monomer-level polymerization ------------------
figure('Color','w'); axis equal off;
filename2 = 'full_network_monomer_level.gif';

collisionThreshNetwork = 0.2;  
blackClusters = [blackClusters; grayClusters];  

for frame = 1:nFrames
    cla; hold on;

    %% Move clusters
    for c = 1:length(blackClusters)
        clusterSize = length(blackClusters(c).members);
        speedFactor = 1/sqrt(clusterSize); 
        blackClusters(c).center = reshape(blackClusters(c).center,1,2) + reshape(blackClusters(c).vel,1,2)*dt*speedFactor;
        blackClusters(c).angle = blackClusters(c).angle + blackClusters(c).rotVel*dt*speedFactor;

        pos = blackClusters(c).center; vel = blackClusters(c).vel;
        if pos(1)<=0 || pos(1)>=boxSize; vel(1)=-vel(1); end
        if pos(2)<=0 || pos(2)>=boxSize; vel(2)=-vel(2); end
        blackClusters(c).vel = vel;
    end

    %% Green bead motion
    for g = 1:nGreen
        greenPos(g,:) = greenPos(g,:) + greenVel(g,:)*dt;
        if greenPos(g,1)<=0 || greenPos(g,1)>=boxSize; greenVel(g,1)=-greenVel(g,1); end
        if greenPos(g,2)<=0 || greenPos(g,2)>=boxSize; greenVel(g,2)=-greenVel(g,2); end
    end

    %% Monomer-level sticking
    i = 1;
    while i < length(blackClusters)
        j = i + 1;
        while j <= length(blackClusters)
            rel_i = sanitizeRelPos(blackClusters(i).relPos);
            rel_j = sanitizeRelPos(blackClusters(j).relPos);

            pos_i = rel_i + repmat(reshape(blackClusters(i).center,1,2), size(rel_i,1), 1);
            pos_j = rel_j + repmat(reshape(blackClusters(j).center,1,2), size(rel_j,1), 1);

            if any(pdist2(pos_i, pos_j) < collisionThreshNetwork, 'all')
                allMembers = [blackClusters(i).members, blackClusters(j).members];
                newRelPos = generateRelPosPolymer(allMembers, spacing);

                blackClusters(i).members = allMembers;
                blackClusters(i).relPos = newRelPos;
                blackClusters(i).vel = mean([blackClusters(i).vel; blackClusters(j).vel],1);
                blackClusters(i).center = mean([blackClusters(i).center; blackClusters(j).center],1);
                blackClusters(j) = [];
            else
                j = j + 1;
            end
        end
        i = i + 1;
    end

    %% Plot clusters
    for c = 1:length(blackClusters)
        theta = blackClusters(c).angle;
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        relPos = sanitizeRelPos(blackClusters(c).relPos);
        clusterPos = (R*relPos')' + repmat(reshape(blackClusters(c).center,1,2), size(relPos,1), 1);
        grayIdx = blackClusters(c).members > nBlack;
        blackIdx = ~grayIdx;
        if any(blackIdx)
            plot(clusterPos(blackIdx,1), clusterPos(blackIdx,2), 'ko','MarkerFaceColor','k','MarkerSize',6);
        end
        if any(grayIdx)
            plot(clusterPos(grayIdx,1), clusterPos(grayIdx,2), 'ko','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',8);
        end
    end

    plot(greenPos(:,1), greenPos(:,2), 'go','MarkerFaceColor','g','MarkerSize',6);
    xlim([0 boxSize]); ylim([0 boxSize]);
    drawnow;

    frameImg = getframe(gcf);
    [A,map] = rgb2ind(frame2im(frameImg),256);
    if frame == 1
        imwrite(A,map,filename2,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename2,'gif','WriteMode','append','DelayTime',0.1);
    end
end

%% ------------------ Helper functions ------------------
function relPos = generateRelPosPolymer(members, spacing)
    % Modified: generate branched, web-like clusters
    N = length(members);
    relPos = zeros(N,2);
    relPos(1,:) = [0 0];
    if N == 2
        relPos(2,:) = [spacing 0];
        return;
    end

    for k = 2:N
        % Randomly pick an existing node to branch from
        parentIdx = randi([1, k-1]);

        % Random branching angle (0–360°)
        ang = rand()*360;

        % Random radial length (with slight noise)
        r = spacing*(0.9 + 0.2*rand());

        % Position new monomer near chosen parent
        relPos(k,:) = relPos(parentIdx,:) + r*[cosd(ang) sind(ang)];
    end

    % Slight jitter to prevent perfect symmetry
    relPos = relPos + 0.02*randn(size(relPos));
end

function rel = sanitizeRelPos(rel)
    if isempty(rel)
        rel = [0 0];
    elseif iscell(rel)
        rel = cell2mat(rel);
    end
    rel = double(rel);
    rel = reshape(rel, [], 2);
end
