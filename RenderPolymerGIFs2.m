clear; clc; close all;

%% ------------------ PARAMETERS ------------------
nChains = 60;              % number of initial polymer chains
chainLengthRange = [5 25]; % each chain has 5–25 monomers
nGray = 100;               % gray branch points
nGreen = 40;               % green moving probes
nFrames = 150;             % total frames per video
boxSize = 10;              % simulation box size
dt = 0.9;                  % timestep
collisionThresh = 0.15;    % distance threshold for merging
spacing = 0.15;            % monomer spacing
branchConnectRadius = 0.25;% how close a chain end must be to a branch point to attach
maxBranchLinks = 3;        % each gray node can connect up to 3 chains

%% ------------------ INITIALIZE CHAINS ------------------
chains = [];
monomerCounter = 0;
for i = 1:nChains
    L = randi(chainLengthRange);
    monomerIdx = monomerCounter + (1:L);
    monomerCounter = monomerCounter + L;
    relPos = zeros(L,2);
    for k = 2:L
        relPos(k,:) = relPos(k-1,:) + spacing*[cosd(rand*360) sind(rand*360)];
    end
    chains(i).members = monomerIdx;
    chains(i).relPos = relPos - mean(relPos); % center chain
    chains(i).center = boxSize*rand(1,2);
    chains(i).vel = 0.05*randn(1,2);
    chains(i).angle = 0;
    chains(i).rotVel = 0.2*randn;
    chains(i).isAttached = false;
end
nBlack = monomerCounter;

%% ------------------ INITIALIZE GRAY BRANCH POINTS ------------------
grayPos = boxSize*rand(nGray,2);
grayVel = 0.02*randn(nGray,2);
grayLinks = cell(nGray,1); % store which chains are linked to each gray node

%% ------------------ INITIALIZE GREEN PROBES ------------------
greenPos = boxSize*rand(nGreen,2);
greenVel = 0.08*randn(nGreen,2);

%% ------------------ VIDEO 1: Initial chain motion ------------------
figure('Color','w'); axis equal off;
filename1 = 'chains_free_motion.gif';

for frame = 1:nFrames
    cla; hold on;

    % Move chains
    for c = 1:length(chains)
        speedFactor = 1/sqrt(length(chains(c).members));
        chains(c).center = reshape(chains(c).center,1,2) + reshape(chains(c).vel,1,2)*dt*speedFactor;
        chains(c).angle = chains(c).angle + chains(c).rotVel*dt*speedFactor;

        pos = chains(c).center; vel = chains(c).vel;
        if pos(1)<=0 || pos(1)>=boxSize; vel(1)=-vel(1); end
        if pos(2)<=0 || pos(2)>=boxSize; vel(2)=-vel(2); end
        chains(c).vel = vel;
    end

    % Move gray and green
    grayPos = grayPos + grayVel*dt;
    greenPos = greenPos + greenVel*dt;
    for g = 1:nGray
        if grayPos(g,1)<=0 || grayPos(g,1)>=boxSize; grayVel(g,1)=-grayVel(g,1); end
        if grayPos(g,2)<=0 || grayPos(g,2)>=boxSize; grayVel(g,2)=-grayVel(g,2); end
    end
    for g = 1:nGreen
        if greenPos(g,1)<=0 || greenPos(g,1)>=boxSize; greenVel(g,1)=-greenVel(g,1); end
        if greenPos(g,2)<=0 || greenPos(g,2)>=boxSize; greenVel(g,2)=-greenVel(g,2); end
    end

    % Plot
    for c = 1:length(chains)
        theta = chains(c).angle;
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        relPos = (R*chains(c).relPos')';
        absPos = relPos + repmat(chains(c).center, size(relPos,1),1);
        plot(absPos(:,1), absPos(:,2), 'k.-','MarkerSize',10,'LineWidth',1);
    end
    plot(grayPos(:,1), grayPos(:,2), 'ko','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',8);
    plot(greenPos(:,1), greenPos(:,2), 'go','MarkerFaceColor','g','MarkerSize',6);
    xlim([0 boxSize]); ylim([0 boxSize]); axis equal off;
    drawnow;

    frameImg = getframe(gcf);
    [A,map] = rgb2ind(frame2im(frameImg),256);
    if frame == 1
        imwrite(A,map,filename1,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename1,'gif','WriteMode','append','DelayTime',0.1);
    end
end

%% ------------------ VIDEO 2: NETWORK FORMATION ------------------
filename2 = 'network_formation.gif';
figure('Color','w'); axis equal off;

for frame = 1:nFrames
    cla; hold on;

    % Move unattached chains
    for c = 1:length(chains)
        if ~chains(c).isAttached
            speedFactor = 1/sqrt(length(chains(c).members));
            chains(c).center = reshape(chains(c).center,1,2) + reshape(chains(c).vel,1,2)*dt*speedFactor;
            chains(c).angle = chains(c).angle + chains(c).rotVel*dt*speedFactor;
            pos = chains(c).center; vel = chains(c).vel;
            if pos(1)<=0 || pos(1)>=boxSize; vel(1)=-vel(1); end
            if pos(2)<=0 || pos(2)>=boxSize; vel(2)=-vel(2); end
            chains(c).vel = vel;
        end
    end

    % Move green beads
    for g = 1:nGreen
        greenPos(g,:) = greenPos(g,:) + greenVel(g,:)*dt;
        if greenPos(g,1)<=0 || greenPos(g,1)>=boxSize; greenVel(g,1)=-greenVel(g,1); end
        if greenPos(g,2)<=0 || greenPos(g,2)>=boxSize; greenVel(g,2)=-greenVel(g,2); end
    end

    % Try attaching chain ends to gray nodes
    for c = 1:length(chains)
        if chains(c).isAttached, continue; end

        theta = chains(c).angle;
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        relPos = (R*chains(c).relPos')';
        absPos = relPos + repmat(chains(c).center, size(relPos,1),1);

        chainEnds = [absPos(1,:); absPos(end,:)];

        for e = 1:2
            distances = sqrt(sum((grayPos - chainEnds(e,:)).^2,2));
            [minDist, minIdx] = min(distances);
            if minDist < branchConnectRadius && length(grayLinks{minIdx}) < maxBranchLinks
                % Attach chain to this gray node
                chains(c).isAttached = true;
                grayLinks{minIdx}{end+1} = c;
                % Move chain center to align the end
                shiftVec = grayPos(minIdx,:) - chainEnds(e,:);
                chains(c).center = chains(c).center + shiftVec;
                break;
            end
        end
    end

    % Plot chains and branches
    for c = 1:length(chains)
        theta = chains(c).angle;
        R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
        relPos = (R*chains(c).relPos')';
        absPos = relPos + repmat(chains(c).center, size(relPos,1),1);
        if chains(c).isAttached
            plot(absPos(:,1), absPos(:,2), 'r.-','MarkerSize',10,'LineWidth',1.5);
        else
            plot(absPos(:,1), absPos(:,2), 'k.-','MarkerSize',8,'LineWidth',1);
        end
    end

    % Draw gray nodes and connecting lines
    for g = 1:nGray
        plot(grayPos(g,1), grayPos(g,2), 'ko','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',8);
        if ~isempty(grayLinks{g})
            for linkIdx = 1:length(grayLinks{g})
                c = grayLinks{g}{linkIdx};
                theta = chains(c).angle;
                R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
                relPos = (R*chains(c).relPos')';
                absPos = relPos + repmat(chains(c).center, size(relPos,1),1);
                endPos = absPos(1,:);
                if norm(grayPos(g,:) - absPos(end,:)) < norm(grayPos(g,:) - absPos(1,:))
                    endPos = absPos(end,:);
                end
                plot([grayPos(g,1) endPos(1)], [grayPos(g,2) endPos(2)], 'Color',[0.3 0.3 0.3]);
            end
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
