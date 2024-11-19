clc
close all
expTime = 0.010;

for ii = 1:length(trackRes.traces)
    tPlot = [];
    if size(trackRes.traces{ii,1}, 1) > 19
        currTrace = trackRes.traces{ii,1};
        colPlot = currTrace.col;
        rowPlot = currTrace.row;
        zPlot   = currTrace.z;
        tPlot   = currTrace.t*expTime;
        plot3(colPlot,rowPlot,zPlot)
        %plot with time color coding
        patch([colPlot(:)' nan],[rowPlot(:)' nan],[zPlot(:)' nan],[tPlot(:)' nan],'EdgeColor','interp','FaceColor','none', 'LineWidth',1)
        grid on
        hold on
    end
end

view(50, 35)