function colony = processOneColonyImage(filename, dataDir, varargin)
    % colony = processOib(filename, dataDir, varargin)
    % 
    % process micropatterns from Oib files
    % this is like processVsi, with the main difference that it assume one
    % colony per oib file, and we pass the colID
    % whereas processVsi can deal with multiple colonies per file

    % colMIPDir = fullfile(colDir,'MIP');
    % if ~exist(colMIPDir,'dir')  mkdir(colMIPDir);   end

    % output dirs
    %---------------------
    colDir = fullfile(dataDir,'colonies');
    if ~exist(colDir,'dir')     mkdir(colDir);      end

    previewDir = fullfile(colDir,'previews');
    if ~exist(previewDir,'dir')  mkdir(previewDir); end

    % input parameters
    %---------------------
    in_struct = varargin2parameter(varargin);
    load(fullfile(dataDir,'metaData.mat'),'meta');
    
    type = 'MIP';
    s = round(20/meta.xres);
    colID = '';
    adjustmentFactor = [];
    
    if isfield(in_struct,'type')
        type = in_struct.type;
    end
    if isfield(in_struct,'colID')
        colID = in_struct.colID;
    end
    if isfield(in_struct,'adjustmentFactor')
        adjustmentFactor = in_struct.adjustmentFactor;
    end
    if isfield(in_struct,'cleanScale')
        s = in_struct.cleanScale;
    end
    
    if isfield(in_struct,'clparameters')
        clparameters = in_struct.clparameters;
    else
        clparameters =  struct('sclose', 6, 'sopen', 8, 'checkcontained', false,...
                            'minArea', [],'convhull', true);
    end
    
    if isfield(in_struct,'DAPIChannel')
        DAPIChannel = in_struct.DAPIChannel;
    else
        error('specify DAPI channel');
    end
    
    if isfield(in_struct,'thresh')
        thresh = in_struct.thresh;
    end
    
    if isfield(in_struct,'prenormImage')
        prenormImage = in_struct.prenormImage;
    else
        prenormImage = false;
    end

    %---------------------
    
    disp(['processing ' fullfile(dataDir,filename)]);
    img = readStack2(fullfile(dataDir,filename));

    disp('determine threshold');
    if strcmp(type,'MIP')
        IP = max(img,[],4);
    elseif strcmp(type,'SIP')
        IP = sum(img,4);
        % adjust default for SIP
        if isempty(adjustmentFactor) adjustmentFactor = 0.5; end
    elseif strcmp(type,'plane')
        if isfield(in_struct,'planeIdx')
            planeIdx = in_struct.planeIdx;
        else
            planeIdx = ceil(size(img,4)/2);
        end
        IP = img(:,:,:,planeIdx);
    else
        error('specify SIP or MIP for type');
    end
    nuclei = IP(:,:,DAPIChannel);
    t = thresholdMP(nuclei, adjustmentFactor);
    mask = nuclei > t;
    
    disp('find colonies');
    [colony, cleanmask, welllabel] = findColonies(mask, [], meta, clparameters);
    
    disp('process individual colonies')
    b = colony.boundingBox;
    
    b(b < 1) = 1; %AW added 
    
    si = size(squeeze(img(:,:,1)));
    
    if b(4) > si(1)
        b(4) = si(1);
    end
    
    if b(2) > si(2)
        b(2) = si(2);
    end
    
    colnucmask = mask(b(3):b(4),b(1):b(2));
    colimg = IP(b(3):b(4),b(1):b(2),:);
    
    if(prenormImage)
        dapi = colimg(:,:,DAPIChannel);
        for ii = 1:size(colimg,3)
            colimg(:,:,ii) = normalizeImageToDAPI(colimg(:,:,ii),dapi,3*clparameters.sopen);
        end
    end
    
    if exist('s','var')
        colimg = backGroundSubImOpen(colimg,4*s);
        %colimg = imfilter(colimg,fspecial('gaussian',s,floor(s/3)));
    end
    
    if exist('thresh','var')
        colony.makeRadialAvgNoSeg(colimg, colnucmask,[], meta.colMargin,1,thresh)
    else
        colony.makeRadialAvgNoSeg(colimg, colnucmask,[], meta.colMargin)
    end
    
    disp('save mask');
    imshow(cat(3,mask,cleanmask,0*mask))
    bbox = b;
    rec = [bbox(1), bbox(3), bbox(2)-bbox(1), bbox(4)-bbox(3)];
    rectangle('Position',rec,'LineWidth',2,'EdgeColor','b')
    saveas(gcf, fullfile(previewDir,['col_id' num2str(colID) '_mask.tif']));
    
%     disp('save preview');
%     preview = double(colimg(:,:,setdiff(1:meta.nChannels,DAPIChannel)));
%     for i = 1:size(preview,3)
%         preview(:,:,i) = imadjust(mat2gray(preview(:,:,i)));
%     end
    %imwrite(preview, fullfile(previewDir,['col_id' num2str(colID) '_' type '_preview.tif']));
    
%     % write DAPI z-stack separately for Ilastik    
%     MIP = false;
%     colonies(coli).saveImage(img, colDir, DAPIChannel, MIP);
%     MIP = true;
%     colonies(coli).saveImage(img, colMIPDir, DAPIChannel, MIP);

end