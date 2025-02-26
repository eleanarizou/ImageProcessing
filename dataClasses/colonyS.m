classdef colonyS
    %% stores colony path, radius and radialprofile (in microns) for each colony
    
    %AW note: SC function. Duplicates IH "Colony" data structure. Need to
    %merge
    
    properties
        id
        fileFolder
        fileName
        radius % microns
        radialProfile
    end
    
%     properties(SetAccess = private)
%         radialProfile
%         % structure - Non-normalized, dapiNormalized (values(mean,std, stdError), pixels)
%     end
    
    methods
        % contstructor function
        function object = colonyS(id, fileFolder, fileName,  radius)
            if nargin > 0
                object.id = id;
                object.fileFolder = fileFolder;
                object.fileName = fileName;
                object.radius = radius;
            end
        end
        
        function mask = saveNuclearMembraneMask(object, meta, channelName)
            %% --------------- save nuclear/membrane mask from ilastik segmentation
            % channelName = DAPI/betacatenin;
            
            filePrefix = strtok(object.fileName, '.');
            dapiChannel = find(cell2mat(cellfun(@(c)strcmp(c,channelName),upper(meta.channelNames),'UniformOutput',false)));
            maskName = [filePrefix '_ch' int2str(dapiChannel) '_Simple Segmentation.h5'];
            
            mask = readIlastikFile([object.fileFolder filesep maskName]);
            if strcmpi('DAPI', channelName)
                location = 'nuclear';
                sizeThreshold = floor(0.1*pi*(10^2)/meta.xres);
                mask = bwareaopen(mask, sizeThreshold);
            else
                location = 'membrane';
                mask = ~mask;
            end
            
            colonyMask = imread([object.fileFolder filesep 'colony' int2str(object.id) '_colonyMask.tif']);
            mask = mask & colonyMask;
            imwrite(mask, [object.fileFolder filesep 'colony' int2str(object.id) '_' location 'Mask.tif']);
        end
        
        
        %% =============================================================================================================
        %% =============================================================================================================
        
        function object = calculateRadialProfile(object, meta, bins, tooHighIntensity, colonyMask, nuclearMask)
            %% ---------------- calculate radial profile for a static colony
            % 1) read images, masks.
            %images = uint16(zeros(meta.xSize, meta.ySize, meta.nChannels));
            for ii = 1:meta.nChannels
                images(:,:,ii) = imread([object.fileFolder filesep 'colony' int2str(object.id) '.tif'], ii);
            end
            
            if ~exist('colonyMask', 'var')
                colonyMask = imread([object.fileFolder filesep 'colony' int2str(object.id) '_colonyMask.tif']);
            end
            
            if ~isempty(find(colonyMask, 1)) % for regular dishes some images are empty
                
                if ~exist('nuclearMask', 'var')
                    %nuclearMask = imread([object.fileFolder filesep 'colony' int2str(object.id) '_nuclearMask.tif']);
                    nuclearMask = saveNuclearMembraneMask(object, meta, 'DAPI');
                end
                nonNuclearMask = colonyMask &~ nuclearMask;
                
                % membrane mask, radial average
                membraneChannel = find(ismember(upper(meta.channelLabel), 'NONMEMBRANE'), 1);
                if ~isempty(membraneChannel)
                    if ~exist('membraneMask', 'var')
                        membraneMask = saveNuclearMembraneMask(object, meta, 'BETACATENIN');
                    end
                    membraneMask = colonyMask & membraneMask;
                end
                
                % 2) get nuclear/non-nuclear/non-membrane/spots intensity
                % profiles (for spots, use the entire colony).
                % select approprite mask for each channel
                
                intensityProfile = zeros(size(images,1), size(images,2), meta.nChannels);
                for ii = 1:meta.nChannels
                    masktype = meta.channelLabel{ii};
                    switch upper(masktype)
                        case 'NUCLEAR'
                            mask1 = nuclearMask;
                        case 'NONNUCLEAR'
                            mask1 = nonNuclearMask;
                        case 'NONMEMBRANE'
                            mask1 = colonyMask & ~membraneMask;
                        case 'SPOTS' % fish spots, use entire colony
                            mask1 = colonyMask;
                    end
                    intensityProfile(:,:,ii) = getColonyNuclearIntensityProfile(images(:,:,ii), mask1, meta, tooHighIntensity);
                end
                
                nuclearChannels = ismember(upper(meta.channelLabel), 'NUCLEAR');
                dapiChannel = find(ismember(upper(meta.channelNames), 'DAPI'), 1);
                
                intensityProfile(:,:,dapiChannel) = intensityProfile(:,:,dapiChannel)./max(max(intensityProfile(:,:,dapiChannel)));
                intensityProfile_dapiNormalized = intensityProfile;
                intensityProfile_dapiNormalized(:,:,nuclearChannels) = intensityProfile(:,:,nuclearChannels)./intensityProfile(:,:,dapiChannel);
                intensityProfile_dapiNormalized(isnan(intensityProfile_dapiNormalized)) = 0;
                
                
            end
            
            % 3) radial average, aranged from edge to center
            object.radialProfile.bins = bins;
            dists = bwdist(~colonyMask);
            dists = dists*meta.xres; % convert to microns
            
            for ii = 1:length(bins)-1
                bin1 = [bins(ii) bins(ii+1)];
                idx1 = find(dists>bin1(1) & dists<=bin1(2));
                object.radialProfile.pixels(1,ii) = numel(idx1);
                
                if~isempty(idx1)
                    for jj = 1:meta.nChannels
                        image1 = intensityProfile(:,:,jj);
                        object.radialProfile.notNormalized.mean(jj,ii) = mean(image1([idx1]));
                        object.radialProfile.notNormalized.std(jj,ii) = std(image1([idx1]));
                        
                        image1 = intensityProfile_dapiNormalized(:,:,jj);
                        image1(isinf(image1)) = 0;
                        object.radialProfile.dapiNormalized.mean(jj,ii) = mean(image1([idx1]));
                        object.radialProfile.dapiNormalized.std(jj,ii) = std(image1([idx1]));
                    end
                end
            end
            %             % membraneChannel
            %             if ~isempty(membraneChannel)
            %                 umToPixel = 1/meta.xres;
            %                 [rA_mean, rA_std, nPixels] =   radialAverageOneColonyOnetimePoint_nonMembrane(colonyMask, membraneMask, images(:,:,membraneChannel), bins, umToPixel);
            %                 object.radialProfile.notNormalized.mean(membraneChannel,:) = rA_mean;
            %                 object.radialProfile.notNormalized.std(membraneChannel,:) = rA_std;
            %                 object.radialProfile.membranePixels(1,:) = nPixels;
            %
            %                 object.radialProfile.dapiNormalized.mean(membraneChannel,:) = rA_mean;
            %                 object.radialProfile.dapiNormalized.std(membraneChannel,:) = rA_std;
            %             end
        end
        
        
        %% =============================================================================================================
        %% =============================================================================================================
        
        function object = calculateRadialProfileD(object, meta, bins, cMask_timepoints)
            %%  ----------------------- calculate radial profiles for a dynamic colony
            % cMask_timepoints: timepoints for which colony masks are saved
            % --- 1) read images, masks
            ch = find(cell2mat(cellfun(@(c)strcmp(c,'BETACATENIN'),upper(meta.channelNames),'UniformOutput',false)));
            prefix = strtok(object.fileName, '.');
            for ii = 1:meta.nTime
                images(:,:,ii) = imread([object.fileFolder filesep prefix '_ch' int2str(ch) '.tif'], ii);
            end
            
            membraneMasks = readIlastikFile([object.fileFolder filesep prefix '_ch' int2str(ch) '_Simple Segmentation.h5']);
            for ii = 1:numel(cMask_timepoints)
                colonyMasks(:,:,ii) = imread([object.fileFolder filesep prefix '_colonyMask' int2str(cMask_timepoints(ii)) '.tif']);
            end
            %%
            % --- 2) make colonyMaksId array - to know when to use which
            % mask
            timeIds = [cMask_timepoints meta.nTime+1];
            counter = 1;
            for ii = 1:numel(timeIds)-1
                colonyMasksId(:, counter:timeIds(ii+1)-1) = ii;
                counter = timeIds(ii+1);
            end
            %%
            % -- 3) calculate radialProfile
            umToPixel = 1/meta.xres;
            for ii = 1:meta.nTime
                colonyImage1 = images(:,:,ii);
                colonyMask1 = colonyMasks(:,:,colonyMasksId(ii));
                membraneMask1 = membraneMasks(:,:,ii).*colonyMask1;
                [rA_mean, rA_std, nPixels] =   radialAverageOneColonyOnetimePoint_nonMembrane(colonyMask1, membraneMask1, colonyImage1, bins, umToPixel);
                object.radialProfile.notNormalized.mean(1,:,ii) = rA_mean;
                object.radialProfile.notNormalized.std(1,:,ii) = rA_std;
                object.radialProfile.pixels(1,:,ii) = nPixels;
            end
            object.radialProfile.bins = bins;
            
        end
        
        function makeNonMembraneMovies(object, meta, ch, colonyMask_timepoints, colonyMasksId)
            %% -------------------------- saves non membrane movies for channel ch of all colonies
            % --- 1) read membrane, colonymasks
            prefix = strtok(object.fileName, '.');
            membraneMasks = readIlastikFile([object.fileFolder filesep prefix '_ch' int2str(ch) '_Simple Segmentation.h5']);
            for ii = 1:numel(colonyMask_timepoints)
                colonyMasks(:,:,ii) = imread([object.fileFolder filesep prefix '_colonyMask' int2str(colonyMask_timepoints(ii)) '.tif']);
            end
            %%
            % --- 2) declare save path
            newFilePath = [object.fileFolder filesep prefix '_ch' int2str(ch) '_nonMembrane.tif'];
            
            % --- 3) read image, make non-membrane cast, save
            for ii = 1:meta.nTime
                colonyImage1 = imread([object.fileFolder filesep prefix '_ch' int2str(ch) '.tif'], ii);
                colonyMask1 = colonyMasks(:,:, colonyMasksId(ii));
                nonMembraneMask1 = colonyMask1 & ~membraneMasks(:,:,ii);
                colonyImage_nonMembrane = bsxfun(@times, colonyImage1, cast(nonMembraneMask1, class(colonyImage1)));
                
                if ii == 1
                    imwrite(colonyImage_nonMembrane, newFilePath);
                else
                    imwrite(colonyImage_nonMembrane, newFilePath, 'WriteMode', 'append');
                end
                
            end
            
            % --
            
            
        end
    end
end