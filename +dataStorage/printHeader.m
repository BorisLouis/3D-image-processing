fileID = fopen('celldata.dat','w');


test = strsplit(frameHead,'<TiffData');

for i =1:length(test)
    fprintf(fileID,['<TiffData' test{i}]);
    fprintf(fileID,'\n');
end

test2 = test{end};
test2 = strsplit(test2,'<Plane DeltaT');
for i =2:length(test2)
    fprintf(fileID,['<Plane DeltaT' test2{i}]);
    fprintf(fileID,'\n');
end

fclose(fileID);