close all
Path = append('D:\Multimodal tracking\20250724\alldata\50mgmL_sample4_min0__1', filesep, '50mgmL_sample4_min0__1_MMStack_Pos0.ome.tif');
[FrameInfo, MovInfo, ~] =  Load.Movie.ome.getInfo(Path);
[ movC1, movC2, idx ] = Load.Movie.ome.load( FrameInfo, MovInfo, 50);
StackToPlot = movC1(42:42+530, 1610:1610+360);
figure()
imagesc(StackToPlot)
colormap('gray');
clim([100 180])
axis image

Path = append('D:\Multimodal tracking\20250724\alldata\50mgmL_sample4_min4__1', filesep, '50mgmL_sample4_min4__1_MMStack_Pos0.ome.tif');
[FrameInfo, MovInfo, ~] =  Load.Movie.ome.getInfo(Path);
[ movC1, movC2, idx ] = Load.Movie.ome.load( FrameInfo, MovInfo, 50);
StackToPlot = movC1(42:42+530, 1610:1610+360);
figure()
imagesc(StackToPlot)
colormap('gray');
clim([100 180])
axis image

Path = append('D:\Multimodal tracking\20250724\alldata\50mgmL_sample4_min8__1', filesep, '50mgmL_sample4_min8__1_MMStack_Pos0.ome.tif');
[FrameInfo, MovInfo, ~] =  Load.Movie.ome.getInfo(Path);
[ movC1, movC2, idx ] = Load.Movie.ome.load( FrameInfo, MovInfo, 50);
StackToPlot = movC1(42:42+530, 1610:1610+360);
figure()
imagesc(StackToPlot)
colormap('gray');
clim([100 180])
axis image

Path = append('D:\Multimodal tracking\20250724\alldata\50mgmL_sample4_min10__1', filesep, '50mgmL_sample4_min8__2_MMStack_Pos0.ome.tif');
[FrameInfo, MovInfo, ~] =  Load.Movie.ome.getInfo(Path);
[ movC1, movC2, idx ] = Load.Movie.ome.load( FrameInfo, MovInfo, 50);
StackToPlot = movC1(42:42+530, 1610:1610+360);
figure()
imagesc(StackToPlot)
colormap('gray');
clim([100 180])
axis image

Path = append('D:\Multimodal tracking\20250724\alldata\50mgmL_sample4_min13__1', filesep, '50mgmL_sample4_min13__1_MMStack_Pos0.ome.tif');
[FrameInfo, MovInfo, ~] =  Load.Movie.ome.getInfo(Path);
[ movC1, movC2, idx ] = Load.Movie.ome.load( FrameInfo, MovInfo, 50);
StackToPlot = movC1(42:42+530, 1610:1610+360);
figure()
imagesc(StackToPlot)
colormap('gray');
clim([100 180])
axis image