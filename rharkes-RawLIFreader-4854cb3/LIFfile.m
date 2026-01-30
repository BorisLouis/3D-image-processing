classdef LIFfile
    % Abbriviations
    % - Leica Image File (lif) <-- what we read here
    % - Extended Leica File (xlef)
    % - Extended Leica Image Files (xlif)
    % - Leica Object File (lof)
    % - Common Block Header (CBH)
    % - Metadata Block Header (MBH)
    % - Binary Block Header (BBH)
    
    properties
        binaryBlockHeaders
        xmlRaw
        xmlData
    end
    properties (Hidden)
        fid
        closed
        filesize
        
    end
    properties (Constant, Hidden)
        CBHi = 0x70;  %Common Block Header identifier (112)
        MBHi = 0x2A;  %Metadata Block Header identifier (42)
        LBBHi = 0x2A; %LIF Binary Block Header identifier (42)
        extension = '.lif'
    end
    
    methods
        function obj = LIFfile(filename)            
            % LIF file contains a Metadata Block followed by LIF Binary Blocks.
            if isstruct(filename) %file as structure
                filename = fullfile(filename.folder,filename.name);
            end
            if ~strcmp(filename(end-3:end),obj.extension),filename=[filename,obj.extension];end
            if ~exist(filename,'file')
                error('cannot find file')
            else
                obj.fid = fopen(filename,'r','l');
            end
            d = dir(filename);
            obj.filesize = d.bytes;
            
            % Read the MetadataBlock data
            obj.xmlRaw = obj.readMetadataBlock('First Metadata Block');
            xmlStream = java.io.ByteArrayInputStream(uint8(obj.xmlRaw));
            xmlSource = org.xml.sax.InputSource(xmlStream);
            xmlData = obj.parseLIFXML(xmlread(xmlSource));
            
            % Get information about the binary blocks
            blockCounter=1;
            while ftell(obj.fid)<obj.filesize
                header = getLIFBinaryBlockHeader(obj,sprintf('Binary Block %.0f', blockCounter));
                binaryBlockHeaders(blockCounter) = header;
                blockCounter=blockCounter+1;
            end
            
            %couple binaryBlockHeaders to XML blocks
            for ct = 1:length(binaryBlockHeaders)
                element = xmlData(ismember(xmlData(:,2),binaryBlockHeaders(ct).identifier),1);
                binaryBlockHeaders(ct).metadata = element{1};
                %extract the name of the element
                theAttributes = binaryBlockHeaders(ct).metadata.getAttributes;
                numAttributes = theAttributes.getLength;
                for count = 1:numAttributes
                    attrib = theAttributes.item(count-1);
                    if strcmp(char(attrib.getName),'Name')
                        binaryBlockHeaders(ct).Name=char(attrib.getValue);
                    end
                end  
            end
            obj.binaryBlockHeaders=binaryBlockHeaders;
        end
        function data = getLIFBinaryBlockData(obj,BBid) %get the raw data from a lif file
            if ischar(BBid) %find the name in the description
                BBidName = BBid;
                BBid = find(ismember({obj.binaryBlockHeaders.identifier},BBidName));
                assert(~isempty(BBid),'Did not find binary block with the name %s',BBidName);
            end
            assert(BBid<=length(obj.binaryBlockHeaders),'Could not find binary block %.0f, only %.0f present',BBid,length(obj.binaryBlockHeaders))
            header = obj.binaryBlockHeaders(BBid);
            fseek(obj.fid,header.offset,'bof');
            data = fread(obj.fid,header.datasize,'*uint8');
        end
        function header = getLIFBinaryBlockHeader(obj,name)
            % consists of CBH LifBBH and binary data
            CBHsize = readCommonBlockHeader(obj,name);
            ID = fread(obj.fid,1,'*uint8');
            assert(ID==obj.LBBHi,'Found %x instead of %x in %s. Not a valid block.',ID, obj.LBBHi,name)
            BDsize = fread(obj.fid,1,'*uint64'); %binary data size
            ID = fread(obj.fid,1,'*uint8');
            assert(ID==obj.LBBHi,'Found %x instead of %x in %s. Not a valid block.',ID, obj.LBBHi,name)
            BDIsize = fread(obj.fid,1,'*uint32'); %binary data identifier size
            BDI = char(fread(obj.fid,BDIsize,'*uint16')); %binary data indentifier
            offset = ftell(obj.fid);
            fseek(obj.fid,BDsize,'cof'); %forward over BDsize
            header.datasize = BDsize;
            header.offset = offset;
            header.identifier = BDI';
        end
        function mdb = readMetadataBlock(obj,name)
            %consist of CBH MDB and metadata
            CBHsize = obj.readCommonBlockHeader(name);
            MBHsize = obj.readMetadataBlockHeader(name);
            mdb = char(fread(obj.fid,MBHsize,'*uint16'));
        end
        function CBHsize = readCommonBlockHeader(obj,name)
            ID = fread(obj.fid,1,'*uint32');
            assert(ID==obj.CBHi,'Found %x instead of %x in %s. Not a valid block.',ID, obj.CBHi,name)
            CBHsize = fread(obj.fid,1,'*uint32');
        end
        function MBHsize = readMetadataBlockHeader(obj,name)
            ID = fread(obj.fid,1,'*uint8');
            assert(ID==obj.MBHi,'Found %x instead of %x in %s. Not a valid block.',ID, obj.MBHi,name)
            MBHsize = fread(obj.fid,1,'*uint32');
        end
        function close(obj)
            fclose(obj.fid);
        end
    end
    
    methods(Static) %taken from Matlab https://nl.mathworks.com/help/matlab/ref/xmlread.html
        function elements = parseLIFXML(theNode)
            elements = findElements(theNode,{});
            function elements = findElements(theNode,elements)
                % Adapted from the ParseXML function
                if strcmp(char(theNode.getNodeName),'Element')
                    temp = theNode.cloneNode(true);
                    %remove the child-node called Children
                    childNodes = temp.getChildNodes;
                    numChildNodes = childNodes.getLength;
                    for count = 1:numChildNodes
                        theChild = childNodes.item(count-1);
                        if strcmp(char(theChild.getNodeName),'Children')
                            elements = findElements(theChild,elements); %Search the Children
                            temp.removeChild(theChild);
                        end
                        if strcmp(char(theChild.getNodeName),'Memory')
                            theAttributes = theChild.getAttributes;
                            numAttributes = theAttributes.getLength;
                            for count = 1:numAttributes
                                attrib = theAttributes.item(count-1);
                                if strcmp(char(attrib.getName),'MemoryBlockID')
                                    MemoryBlockID=char(attrib.getValue);
                                end
                            end
                        end
                    end
                    elements{end+1,1}=temp;
                    elements{end,2}=MemoryBlockID;
                elseif theNode.hasChildNodes %Look at the children
                    childNodes = theNode.getChildNodes;
                    numChildNodes = childNodes.getLength;
                    for count = 1:numChildNodes
                        theChild = childNodes.item(count-1);
                        elements = findElements(theChild,elements);
                    end
                end
            end
        end
    end
end
