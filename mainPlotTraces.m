%% 
% This code assumes that some trackRes has been loaded onto matlab
% no clearing
clc
close all

%% User input

path2Save = 'D:\Polymer Dynamics\PAA_3x_bAA\time12\20251008_3xbAA_sample1_min12_1';
ext = '.gif';
filename=sprintf('%s%sdatamovie%s', path2Save,filesep,ext);

minSize = 10;%number of frame the traces needs to last to be plotted.
expTime = 0.010; %sec
sizeParticles = 200; % diameter in nm
frameRate = 10;
trailing = 20; %frame the traces stays in the movie

%% Plot all Traces with time color-coding (4D plot)
CM = zeros(size(traces,2),3);
maxFr = zeros(size(traces,2),1);
figure
hold on
for i = 1:size(traces,2)
    
    currTrace = traces{1,i};
    
    if height(currTrace) > minSize
        colPlot = currTrace.col;
        rowPlot = currTrace.row;
        zPlot   = currTrace.z;
        tPlot   = currTrace.t*expTime;
        plot3(colPlot,rowPlot,zPlot)
        %plot with time color coding
        patch([colPlot(:)' nan],[rowPlot(:)' nan],[zPlot(:)' nan],[tPlot(:)' nan],'EdgeColor','interp','FaceColor','none')

    end        
    CM(i,:) = [mean(currTrace.row),mean(currTrace.col),mean(currTrace.z)];
    maxFr(i,:) = max(currTrace.t);
end
CM = mean(CM,1);
maxFr = max(maxFr);
axis image
%% Make Awesome movie

radius = 1000;
yLimit = [0 35000];
xLimit = [0 25000];
zLimit = [-2000 2000];


Fig = figure;
view(3);
xlim(xLimit);
ylim(yLimit);
zlim(zLimit);
xlim manual;
ylim manual;
gcf;
hold on

[x,y,z] = sphere(32);
x = x*sizeParticles/2;
y = y*sizeParticles/2;
z = z*sizeParticles/2;

camlight
lighting('gouraud');
for i = 1 :maxFr
    for j = 1:size(traces,2)
        currTrace = traces{1,j};
        if height(currTrace) > minSize
            idx2Frame = currTrace.t==i;
            idx = i-trailing:i;
            idx = ismember(currTrace.t,idx);
            if and(~all(idx==0), ~all(idx2Frame ==0))

                data2Plot = currTrace(idx,:);
                axis image
                xlim(xLimit);
                ylim(yLimit);
                zlim(zLimit)
                xlim manual;
                ylim manual;
                zlim manual;

                view(3);

                gcf;
                hold on

                plot3(data2Plot.col,data2Plot.row,data2Plot.z,'color',[0 0 1])

                X = x+currTrace.col(idx2Frame);
                Y = y+currTrace.row(idx2Frame);
                Z = z+currTrace.z(idx2Frame);

                surf(X,Y,Z,'LineStyle','none','Facecolor',[0.4,0.4,0.4])
                

            end
        end
    end
    drawnow;
    frame = getframe(Fig);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    if i == 1

        imwrite(imind,cm,filename,'gif','DelayTime',1/frameRate, 'loopcount',inf);

    else

        imwrite(imind,cm,filename,'gif','DelayTime',1/frameRate, 'writemode','append');

    end
    clf;

end

%% plot only a single traced
figure
idx = 7;

currTrace = trackRes.traces{idx,1};
    
colPlot = currTrace.col;
rowPlot = currTrace.row;
zPlot   = currTrace.z;
tPlot   = currTrace.t*expTime;

patch([colPlot(:)' nan],[rowPlot(:)' nan],[zPlot(:)' nan],[tPlot(:)' nan],'EdgeColor','interp','FaceColor','none')

view(3)
colorbar
xlabel('Position (nm)')
ylabel('Position (nm)')
zlabel('Position (nm)')
axis image