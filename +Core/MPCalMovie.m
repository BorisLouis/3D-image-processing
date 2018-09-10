classdef MPCalMovie < Core.MPParticleMovie
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        traces
    end
    
    methods
        function obj = MPCalMovie(raw,cal)
            
            obj  = obj@Core.MPParticleMovie(raw,cal);
        end
        
        function [traces, counter] = trackInZ(obj,trackParam)
            %track the particle in the Z direction (3rd dimension here)
            %Here we do not expect any big movement from one frame to the
            %other so we give a warning if the tracking parameter seems to
            %soft.
            assert(isstruct(trackParam),'Tracking parameter is expected to be a struct with two field "euDistPx" and "ellip"')
            assert(and(isfield(trackParam,'euDistPx'),isfield(trackParam,'ellip')),...
                'Tracking parameter is expected to be a struct with two field "euDistPx" and "ellip"')
            
            if trackParam.euDistPx > 1
                warning('Current euclidian distance thresholds is probably too high, we do not expect much movement from one frame to the next here')
            end
            
            if or(trackParam.ellip > 6, trackParam.ellip<=3)
                warning('Requested ellipticity thresold is better to be close to 5 which means that at least 2 of the best focus plane should be consistent ([1 2 3 2 1])');
            end
            
            [traces,counter] = obj.zTracking(trackParam);
            
            obj.traces.trace = traces;
            obj.traces.nTrace = counter;
            
        end
       
    end
        
    methods (Access = private)
        function [traces,counter] = zTracking(obj, trackParam)
            %track the particle in the Z direction (3rd dimension here)
            assert(~isempty(obj.calibrated),'Data should be calibrated to do ZzCalibrationration');
            assert(~isempty(obj.candidatePos), 'No candidate found, please run findCandidatePos before zzCalibrationration');
            assert(~isempty(obj.particles), 'No particles found, please run superResConsolidate method before doing ZzCalibrationration');
            assert(isstruct(trackParam),'Tracking parameter is expected to be a struct with two field "euDistPx" and "ellip"')
            assert(and(isfield(trackParam,'euDistPx'),isfield(trackParam,'ellip')),...
                'Tracking parameter is expected to be a struct with two field "euDistPx" and "ellip"')
            %We copy the List as boolean to keep track of where there are
            %still particles left
            [listBool] = Core.trackingMethod.copyList(obj.particles.List,1);
            %We copy as NaN for storage of the traces;
            [traces]   = Core.trackingMethod.copyList(obj.particles.List,NaN);
            %We pick the first particle available
            [idx] = Core.trackingMethod.pickParticle(listBool);
            counter = 1;
            errCount =1;
            while (idx)
                %loop until there is no particle (pickParticle return false)
                if errCount>1000
                    warning('While loop ran for unexpectedly longer time');
                    break;
                    
                end
                %Connect particles (cf consolidation but across frames
                [listIdx] = Core.trackingMethod.connectParticles(obj.particles.List,listBool,idx, trackParam);
                %if the particle was connected in less than 5 frames we remove
                % its appearance from the list bool
                if length(listIdx) < 5
                    
                    [listBool] = Core.trackingMethod.removeParticles(listBool,listIdx);
                    
                else
                    %Otherwise we store traces, increment counter and remove.
                    [traces]  = Core.trackingMethod.storeTraces(traces,listIdx,counter);
                    counter = counter +1;
                    [listBool] = Core.trackingMethod.removeParticles(listBool,listIdx);
                    
                end
                % We pick a new particle and start all over again
                [idx] = Core.trackingMethod.pickParticle(listBool);
                errCount = errCount +1;
            end
            counter = counter -1;
        end
    end
end
