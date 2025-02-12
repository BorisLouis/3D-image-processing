# RawLIFreader
Reading the raw data from a Leica lif-file in Matlab for further analysis.
Adapted from [bioformats](https://github.com/ome/bioformats/blob/develop/components/formats-gpl/src/loci/formats/in/LIFReader.java). 

Bioformats will be more user friendly than the RawLIFreader. This small script is only to make it easier to access the bytes in the .lif file directly. If you don't need that, use bioformats. Or the [load Leica lif file](https://nl.mathworks.com/matlabcentral/fileexchange/48774-load-leica-lif-file) by [Ron Hoebe](https://nl.mathworks.com/matlabcentral/profile/authors/91123)

```
LF = LIFfile(filename); %read the content of the .lif file
dat=LF.getLIFBinaryBlockData('MemBlock_6122'); %contains raw data. Find the name of the memoryblock in LF.binaryBlockHeaders
