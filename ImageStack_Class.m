classdef  ImageStack_Class < dynamicprops
% author: Oleg Krichevsky okrichev@bgu.ac.il, okrichev@gmail.com
    
    properties
        gate_ns = [];
        PixSizeXYZ = [];
        PixSizeHVP = [];
        Pic1 = []; 
        Pic2 = [];
        PhotonInd = struct('phInd1', [], 'phInd2', []);
        %ROI = struct('x', [], 'y', []);
        ScanType = 'XYscan';
        
        fpathToBioformatFile = '';
        ImageElementCount
        ImElementNo
        NoTimePoints 
        NoColors
        NoZplanes
        NoPlanesPerTimePoint
        TotalPlanes
        stackSizeX
        stackSizeY
        ColorInd
        TimeInd
        ZInd
        X0
        Y0
        Z0
        TimeIncrement
        ChannelNames = {}
        CLim
        
        PCAcoeff
        CompMatrix
        hfig % figure handle
    end
    
    properties (Hidden = true)
        dynprop_meta = struct();
    end
    
    methods
        function obj = ImageStack_Class(varargin)
            if (~isempty(varargin)),
                if isa(varargin{1}, 'Image_Class'),
                    [roi, IsInList] = ParseInputs('From ROI', struct('x', [], 'y', []), varargin); %if interactive pass [] for ROI
                    prnt = varargin{1};
                    if IsInList, % create from ROI
                        obj.gate_ns = prnt.gate_ns;
                        obj.PixSizeXYZ = prnt.PixSizeXYZ;
                        obj.PixSizeHVP = prnt.PixSizeHVP;
                        figure;
                        if isempty(roi.x),
                            CenterPlane = round((size(prnt.Pic1, 3)+1)/2);
                            imagesc(prnt.Pic1(:, :, CenterPlane));
                            set(gca, 'YDir', 'normal');
                            [BW, xi, yi] = roipoly;
                            obj.ROI.x = xi;
                            obj.ROI.y = yi;
                        else % apply given ROI
                            obj.ROI = roi;
                            BW = roipoly(prnt.Pic1(:, :, 1), roi.x, roi.y);
                        end

                        BW3 = repmat(BW, [1, 1, size(prnt.Pic1, 3)]);
                        obj.Pic1 = prnt.Pic1.*BW3;
                        obj.Pic2 = prnt.Pic2.*BW3;
                        obj.Pic1(isnan(obj.Pic1)&(~BW)) = 0; %otherwise product of NaN and 0 gives Nan
                        obj.Pic2(isnan(obj.Pic2)&(~BW)) = 0;
                    end %if IsInList
                elseif isa(varargin{1}, 'char'), % load from file
                    fpath = varargin{1};
                    [ROIrect, isROIrect] = ParseInputs('ROIrect', [], varargin);
                    s = ParseInputs('SeriesNo', 1, varargin);
                    if isROIrect,
                        obj.X0 = ROIrect(1);
                        obj.Y0 = ROIrect(3);
%                         w = ROIrect(2) - ROIrect(1) + 1;
%                         h = ROIrect(4) - ROIrect(3) + 1;
%                         data = bfopen(fpath, obj.X0, obj.Y0, w, h); % assume this is bioformat data
                    else
%                         data = bfopen_v1(fpath); % assume this is bioformat data
                        obj.X0 = 0;
                        obj.Y0 = 0;
                    end
                    seriesInFile = ImageStack_Class.DoGetSeriesData(fpath);
                    PixSizeX = seriesInFile(s).seriesMetadata.get('Pixel width (in microns)');
                    PixSizeY = seriesInFile(s).seriesMetadata.get('Pixel height (in microns)');
                    PixSizeZ = seriesInFile(s).seriesMetadata.get('Z step (in microns)');
                    obj.PixSizeXYZ = [PixSizeX PixSizeY PixSizeZ];
                    
                    data = bfopen_v1(fpath, varargin{:});
                    obj.fpathToBioformatFile = fpath;
                    NoSeries = size(data, 1);
%                     if (NoSeries > 1),
%                         errordlg('Number of experiments larger than 1! Check what this means');
%                     end;
 
                    omeMeta = data{s, 4};
                    obj.ImageElementCount = omeMeta.getImageCount();
%                     if (obj.ImageElementCount > 1),
%                         errordlg('Image element count larger than 1! Check what this means');
%                     end;
                    obj.ImElementNo = s - 1;
                    ImElementNo = obj.ImElementNo;
                    obj.NoTimePoints = omeMeta.getPixelsSizeT(ImElementNo).getValue();
                    obj.NoColors = omeMeta.getPixelsSizeC(ImElementNo).getValue();
                    obj.NoZplanes = omeMeta.getPixelsSizeZ(ImElementNo).getValue();
                    obj.TotalPlanes = omeMeta.getPlaneCount(ImElementNo);
                    obj.stackSizeX = omeMeta.getPixelsSizeX(ImElementNo).getValue(); % image width, pixels
                    obj.stackSizeY = omeMeta.getPixelsSizeY(ImElementNo).getValue();  
                    obj.NoPlanesPerTimePoint = obj.NoColors*obj.NoZplanes;
                    obj.TimeInd = 1:obj.NoTimePoints;
                    obj.ZInd = 1:obj.NoZplanes;
                    obj.ColorInd = 1:obj.NoColors;
                    disp(['No of colors: ', num2str(obj.NoColors)])
                    
                    %get default channel Names
                    if ~isprop(obj, 'ChannelNames'),
                        obj.addprop('ChannelNames');
                    end
                    
                    for chanNo = 1:obj.NoColors,
                        ChanNames{chanNo} = omeMeta.getChannelName(ImElementNo, chanNo - 1);
                    end
                    
                    if seriesInFile(s).seriesMetadata.containsKey('Channel #4'),
                        DefaultChannelNames = {'SHG', 'GFP', 'Tomato', 'Auto'};
                    elseif seriesInFile(s).seriesMetadata.containsKey('Channel #3'),
                        DefaultChannelNames = {'SHG', 'GFP', 'Tomato'};
                    end
                    
                    obj.ChannelNames = ParseInputs('ChannelNames', DefaultChannelNames, varargin);
                                        
                    %get plane stack
                    
                    if ~strcmp(omeMeta.getPixelsDimensionOrder(ImElementNo), 'XYCZT'),
                        error('Nonstandard Dimension order: adjust the program!');
                    end
                    
%                     obj.PixSizeXYZ = double([omeMeta.getPixelsPhysicalSizeX(ImElementNo).value, ...
%                         omeMeta.getPixelsPhysicalSizeY(ImElementNo).value, ...
%                         omeMeta.getPixelsPhysicalSizeZ(ImElementNo).value]);
                    obj.PixSizeHVP = obj.PixSizeXYZ;
                    timeIncr = omeMeta.getPixelsTimeIncrement(ImElementNo);
                    if ~isempty(timeIncr),
                        obj.TimeIncrement = double(timeIncr.value);
                    end
                    
                    for planeNo = 1:obj.TotalPlanes,
                        planeColor = omeMeta.getPlaneTheC(ImElementNo, planeNo-1).getValue() + 1;
                        planeTime = omeMeta.getPlaneTheT(ImElementNo, planeNo-1).getValue() + 1;
                        planeZ = omeMeta.getPlaneTheZ(ImElementNo, planeNo-1).getValue() + 1;
                        [Color, Z, Time] = obj.NoInStackToCZT(planeNo);
                        if (planeColor == Color) & (planeTime == Time) & (planeZ == Z),
%                            Img = data{s, 1}{planeNo, 1};
%                            obj.Pic1(:, :, planeNo) = Img;
                        else
                            error('Mismatch in Color, Time or Z values!')
                        end
                    end
                    obj.Pic1 = reshape([data{s, 1}{:, 1}], [size(data{s, 1}{1, 1}), obj.TotalPlanes]);
                    
                elseif isa(varargin{1}, 'ImageStack_Class'), % create SubStack
                    wholeStack = varargin{1}; 
                    [planeNo, Color, Z, Time, ChanNames] = wholeStack.CZTtoNoInStack(varargin{2:end});
                    ROIrect = ParseInputs('ROIrect', [], varargin); 
                    if isempty(ROIrect),
                        ROIrect = [1 size(wholeStack.Pic1, 2) 1 size(wholeStack.Pic1, 1)];
                    end
                    
                    obj = ImageStack_Class;
                    %copy properties from wholeStack
                    props = properties(wholeStack);
                    for ii = 1:length(props),
                        if ~isprop(obj, props{ii}),
                            obj.addprop(props{ii});
                        end
                        obj.(props{ii}) = wholeStack.(props{ii});
                    end
                    % correct properties
                    obj.ColorInd = Color;
                    obj.ZInd = Z;
                    obj.TimeInd = Time;
                    obj.NoTimePoints = length(Time);
                    obj.NoZplanes = length(Z);
                    obj.NoColors = length(Color);
                    obj.TotalPlanes = obj.NoTimePoints*obj.NoZplanes*obj.NoColors;
                    obj.stackSizeX = ROIrect(2) - ROIrect(1) + 1;
                    obj.stackSizeY = ROIrect(4) - ROIrect(3) + 1;
                    obj.X0 = ROIrect(1) - 1;
                    obj.Y0 = ROIrect(3) - 1;
                    obj.Pic1 = wholeStack.Pic1(ROIrect(3):ROIrect(4), ROIrect(1):ROIrect(2), planeNo);
                    if ~isempty(wholeStack.Pic2),
                        obj.Pic2 = wholeStack.Pic2(ROIrect(3):ROIrect(4), ROIrect(1):ROIrect(2), planeNo);
                    end
                    obj.NoPlanesPerTimePoint = obj.NoColors*obj.NoZplanes;
                    if isprop(obj, 'ChannelNames')
                        obj.ChannelNames = ChanNames;
                    end;
                end
            end
        end
        
        function [Color, Z, Time] = NoInStackToCZT(obj, planeNo)
            if ~isempty(obj.fpathToBioformatFile), %bioformat image stack
                Time = floor((planeNo - 1)./obj.NoPlanesPerTimePoint) + 1;
                CZ = mod(planeNo - 1, obj.NoPlanesPerTimePoint);
                Z = floor(CZ./obj.NoColors) + 1;
                Color = mod(CZ, obj.NoColors) + 1;
            else 
                error('This function works only for Bioformat data')
            end
        end
                 
        
        function [planeNo, Color, Z, Time, ChanNames] = CZTtoNoInStack(obj, varargin)
            if ~isempty(obj.fpathToBioformatFile), %bioformat image stack
                Color = ParseInputs('Color', [], varargin);
                Time = ParseInputs('Time', [], varargin);
                Z = ParseInputs('Z', [], varargin);
                
                if isempty(Color),
                    Color = 1:obj.NoColors; % all colors
                elseif ischar(Color), % defined by color string
                    Color = find(strcmp(obj.ChannelNames, Color));
                elseif iscell(Color), % cell array of color strings
                    tempColor = [];
                    for i = 1:length(Color),
                        tempColor = [tempColor find(strcmp(obj.ChannelNames, Color(i)))];
                    end
                    Color = tempColor;
                end
                
                if isempty(Time),
                    Time = 1:obj.NoTimePoints; % all time points
                end
                
                if isempty(Z),
                    Z = 1:obj.NoZplanes; % all Z
                end
                
                Zadder = (Z - 1)*obj.NoColors;
                Tadder = (Time - 1)*obj.NoPlanesPerTimePoint;
                
                [COLOR, ZADDER] = meshgrid(Color, Zadder);
                cz = COLOR + ZADDER;
                cz = cz(:);
                
                [CZ, TADDER] = meshgrid(cz, Tadder);
                
                planeNo = CZ + TADDER;
                planeNo = planeNo(:);
                planeNo = sort(planeNo);
                ChanNames = obj.ChannelNames(Color);
            else 
                error('This function works only for Bioformat data')
            end
        end
        
        function SubStack = DoGetSubStack(obj, varargin)
            SubStack = ImageStack_Class(obj, varargin{:});           
        end
        
        function Cln = DoCloneObject(obj)
            Cln = obj.DoGetSubStack;
        end
        
        function RatioStack = DoCreateRatioStack(obj, Color1, Color2, varargin)
            SubStack1 = ImageStack_Class(obj, 'Color', Color1, varargin{:}); 
            SubStack2 = ImageStack_Class(obj, 'Color', Color2, varargin{:});
            SubStack1.Pic1 = double(SubStack1.Pic1)./(double(SubStack2.Pic1) + eps);
            %prepare for type casting
            origType = class(obj.Pic1);
            if ~strcmp(origType, 'double'),
                maxPic = max(SubStack1.Pic1(:));
                rescaleFactor = double(intmax(origType))/maxPic;
                SubStack1.Pic1 = cast(SubStack1.Pic1*rescaleFactor, origType);
            end
            SubStack1.ChannelNames = {[SubStack1.ChannelNames{1} '/' SubStack2.ChannelNames{1}]};
            RatioStack = SubStack1;
        end
        
        function Chan = DoSeparateChannels(obj, varargin)
            Color = ParseInputs('Color', 1:obj.NoColors, varargin);
            Time = ParseInputs('Time', 1:obj.NoTimePoints, varargin);
            for i = 1:length(Color),
                Chan(i) = obj.DoGetSubStack('Color', Color(i), 'Time', Time);
            end
        end
        
        function obj = DoGetStats(obj, varargin)
            for i = 1:size(obj.Pic1, 3),
                StatsStruc.XYplanes.Xcm(i) = sum(sum(obj.Pic1(:, :, i), 1).*(1:obj.stackSizeX))./sum(sum(obj.Pic1(:, :, i)));
                StatsStruc.XYplanes.Ycm(i) = sum(sum(obj.Pic1(:, :, i), 2)'.*(1:obj.stackSizeY))./sum(sum(obj.Pic1(:, :, i)));
                StatsStruc.XYplanes.Z = (1:size(obj.Pic1, 3))*obj.PixSizeXYZ(3);                
            end
            
            for i = 1:size(obj.Pic1, 1),
                StatsStruc.XZplanes.Xcm(i) = sum(sum(obj.Pic1(i, :, :), 3).*(1:obj.stackSizeX))./sum(sum(obj.Pic1(i, :, :)));
                StatsStruc.XZplanes.Zcm(i) = sum(squeeze(sum(obj.Pic1(i, :, :), 2))'.*(1:obj.TotalPlanes))./squeeze(sum(sum(obj.Pic1(i, :, :))));
                StatsStruc.XZplanes.Y = (1:size(obj.Pic1, 1))*obj.PixSizeXYZ(2);
            end
            
            for i = 1:size(obj.Pic1, 2),
                StatsStruc.YZplanes.Ycm(i) = sum(squeeze(sum(obj.Pic1(:, i, :), 3))'.*(1:obj.stackSizeY))./squeeze(sum(sum(obj.Pic1(:, i, :))));
                StatsStruc.YZplanes.Zcm(i) = sum(squeeze(sum(obj.Pic1(:, i, :), 1))'.*(1:obj.TotalPlanes))./squeeze(sum(sum(obj.Pic1(:, i, :))));
                StatsStruc.YZplanes.X = (1:size(obj.Pic1, 2))*obj.PixSizeXYZ(1);
            end


            if ~isprop(obj, 'Stats'),
                obj.addprop('Stats');                
            end
            obj.Stats = StatsStruc;                
        end
        
        
        function obj = DoBinning(obj, varargin)
            BinSize = ParseInputs('BinSize', 2, varargin);
            if size(BinSize) == [1 1],
                BinSize = BinSize*ones(1, 3);
            elseif size(BinSize) == [1 2],
                BinSize = [BinSize 1];
            elseif size(BinSize) ~= [1 3],
                error('BinSize should be a row of 1 to 3 values')
            end
            
            obj.DoBinning1dim(1, BinSize(1));
            obj.DoBinning1dim(2, BinSize(2));
            obj.DoBinning1dim(3, BinSize(3));
            
        end
        
        function obj = DoPCA(obj, varargin)
            T = ParseInputs('Time', 1, varargin);
            Color = ParseInputs('Color', 1:obj.NoColors, varargin)
            
            for i = 1:length(Color),
                P = obj.DoGetSubStack('Color', Color(i), 'Time', T);               
                col(:, i) = P.Pic1(:);
                % sqrt on the first SHG channel
                if i == 1,
                    col(:, i) = sqrt(col(:, i) - min(col(:, i)));
                end
            end
            col = col./repmat(std(col), size(col, 1), 1);
            [COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(col);
            obj.PCAcoeff = COEFF;
        end
        
        function PCAcomp = DoGetPCAcomp(obj, varargin)
            if isempty(obj.PCAcoeff),
                error('run obj.DoPCA first!');
            end

            for C = 1:obj.NoColors,
                P = obj.DoGetSubStack('Color', C, varargin{:});
                col(:, C) = P.Pic1(:);
                if C == 1, % sqrt on SHG channel,
                    col(:, C) = sqrt(col(:, C) - min(col(:, C)));
                end
             end

            col = (col - repmat(mean(col), size(col, 1), 1))./repmat(std(col), size(col, 1), 1);
            SCORE = col*obj.PCAcoeff; 

            for compNo = 1:obj.NoColors;
                comp = SCORE(:, compNo);
                temp = P.DoCloneObject;
                temp.Pic1 = reshape(comp, size(P.Pic1));
                PCAcomp(compNo) = temp;
            end
        end
        
        function obj = DoCompensateChannels(obj, varargin)
            ChannelList = ParseInputs('ChannelList', {'Bone', 'Macrophages', 'HSC'}, varargin);
            for i = 1:length(ChannelList),
            end
        end
        
        function obj = DoShowChannelStacks(obj, varargin)
            OutlierFracCLim = ParseInputs('OutlierFracCLim', 0.02, varargin);
            obj.hfig = figure('Position', get(0, 'Screensize'));
            Z = 1;
            setappdata(obj.hfig,'Z', 1);
            gcf = obj.hfig;
            
            T = ParseInputs('Time', 1, varargin);
            Color = ParseInputs('Color', 1:obj.NoColors, varargin);
            [PCAcomp, IsPCA] = ParseInputs('PCA', 1:obj.NoColors, varargin);
            [ChanRatio, IsChannelRatio] = ParseInputs('Channel Ratio', [3 2], varargin);
            [ColorTitlesIn, IsColorTitlesIn] = ParseInputs('ColorTitles', {}, varargin);
            
            extP = []
            for i = 1:length(varargin),
                if strcmp(class(varargin{i}), 'ImageStack_Class'),
%                     if ~exist('extP'),
%                         extP = varargin{i};
%                     else
                        extP = [extP varargin{i}];
%                    end
                end
            end
             
            
            for i = 1:length(Color),
                P(i) = obj.DoGetSubStack('Color', Color(i), 'Time', T);
                ColorTitle{i} = ['CH' num2str(Color(i))];
            end

            if IsPCA,
                if isempty(obj.PCAcoeff),
                    error('run obj.DoPCA first!');
                end

                for i = 1:length(Color),
                    col(:, i) = P(i).Pic1(:);
                    ColorTitle{i} = ['PCA comp ' num2str(Color(i))];
                end

                col = (col - repmat(mean(col), size(col, 1), 1))./repmat(std(col), size(col, 1), 1);
                SCORE = col*obj.PCAcoeff; 

                for compNo = 1:length(Color);
                    comp = SCORE(:, compNo);
                    P(compNo).Pic1 = reshape(comp, size(P(compNo).Pic1));
                end
            end

            if IsChannelRatio,
               for i = 1:size(ChanRatio, 1),
                   P(length(Color)+i).Pic1 = P(ChanRatio(i, 1)).Pic1./P(ChanRatio(i, 2)).Pic1;
                   ColorTitle{length(Color)+i} = ['CH' num2str(ChanRatio(i, 1)) ' / ' 'CH' num2str(ChanRatio(i, 2))];
               end
            end
            
            P = [P extP];
            if IsColorTitlesIn,
                ColorTitle = ColorTitlesIn;
            end
            
            for i = 1:length(P),
                
                %CLim(i, :) = [min(P(i).Pic1(:)) max(P(i).Pic1(:))];
                if ~isprop(P(i), 'CLim')
                    P(i).addprop('CLim');
                end
                P(i).CLim = prctile(double(P(i).Pic1(:)), [OutlierFracCLim,  1-OutlierFracCLim]*100);
            end
   
            setappdata(obj.hfig,'ChannelData', P);
                        
%             hManager = uigetmodemanager(obj.hfig);
%             try
%                  % this should work for versions of MATLAB <= R2014a
%                set(hManager.WindowListenerHandles, 'Enable', 'off');
%             catch
%                % this works in R2014b, and maybe beyond; your mileage may vary
%                [hManager.WindowListenerHandles.Enabled] = deal(false);
%             end
              % these lines are common to all versions up to R2014b (and maybe beyond)
           % set(obj.hfig, 'WindowKeyPressFcn', @obj.keypress_callback);
            %set(figure_handle, 'KeyPressFcn', @(obj, evt) keypress_cb(obj, evt));

            
            set(obj.hfig, 'KeyPressFcn', {@obj.keypress_callback, @obj.DoShowZplaneInAStack, 'Z', obj.NoZplanes}, ...                
                'NumberTitle', 'off', 'Name', 'Press left or right arrow to navigate or Esc to leave');
            %'WindowKeyPressFcn', [], ...
            
            % set up compensation toolbar
            tBar = uitoolbar(obj.hfig);
            [img,map] = imread('/Users/oleg/Documents/MATLAB/ImageStack OO/round-black-left-arrow.jpg');
            p = uipushtool(tBar,'TooltipString','Previous',...
                             'ClickedCallback',...
                             @obj.CompensationClickedCallback);
            p.CData = img;

            [img,map] = imread('/Users/oleg/Documents/MATLAB/ImageStack OO/round-black-right-arrow.jpg');
            p = uipushtool(tBar,'TooltipString','Next',...
                             'ClickedCallback',...
                             @obj.CompensationClickedCallback);
            p.CData = img;
            
            % Read an image
            [img,map] = imread('/Users/oleg/Documents/MATLAB/ImageStack OO/comp2web.png');

            % Convert image from indexed to truecolor
            %icon = ind2rgb(img,map);

            icon = img;
            % Create a uipushtool in the toolbar
            p = uipushtool(tBar,'TooltipString','Compensation',...
                             'ClickedCallback',...
                             @obj.CompensationClickedCallback);

            % Set the button icon
            p.CData = icon;
            
            obj.DoShowZplaneInAStack(Z);
        end
        
        function RGBpic = DoShowRGB(obj, varargin)
            OutlierFracCLim = ParseInputs('OutlierFracCLim', 0.02, varargin);
            Z = ParseInputs('Z', 1, varargin);
            T = ParseInputs('Time', 1, varargin);
            Color = ParseInputs('Color', 1:obj.NoColors, varargin);
%            [ChanRatio, IsChannelRatio] = ParseInputs('Channel Ratio', [3 2], varargin);
            
            for i = 1:length(Color),
                P = obj.DoGetSubStack('Color', Color(i), 'Time', T);
                RGBpic(:, :, i) = double(P.Pic1(:, :, Z))/max(max(double(P.Pic1(:, :, Z))));
            end

            
   %         extP = []
            for i = 1:length(varargin),
                if strcmp(class(varargin{i}), 'ImageStack_Class'),
                    if i < length(varargin),
                        if iscell(varargin{i+1}),
                            Z = ParseInputs('Z', Z, varargin{i+1});
                            T = ParseInputs('Time', T, varargin{i+1});
                            Color = ParseInputs('Color', 1, varargin{i+1}); 
                        else
                            Color = 1;                         
                        end
                    else
                        Color = 1;
                    end
                    
                    for j = 1:length(Color),
                        P = varargin{i}.DoGetSubStack('Color', Color(j), 'Time', T);
                        RGBpic(:, :, end+1) = double(P.Pic1(:, :, Z))/max(max(double(P.Pic1(:, :, Z))));
                    end  
                end               
            end
            
            if size(RGBpic, 3) > 3,
                warning('There are more than three colors!')
            elseif size(RGBpic, 3) == 2,
                RGBpic(:, :, 3) = 0;
            end
            
                        
%             for i = 1:length(P),
%                 
%                 %CLim(i, :) = [min(P(i).Pic1(:)) max(P(i).Pic1(:))];
%                 if ~isprop(P(i), 'CLim')
%                     P(i).addprop('CLim');
%                 end
%                 P(i).CLim = prctile(double(P(i).Pic1(:)), [OutlierFracCLim,  1-OutlierFracCLim]*100);
%            end

            imagesc(RGBpic);
            axis image;
            figure(gcf)
        end
        
        function keypress_callback(obj, hFig, e, WhatToShowFunc, PlaneParamName, UpperLim, varargin)
            Z = getappdata(hFig, PlaneParamName);
            %getappdata(obj, 'ImageStack_Class_object')
            switch e.Key,
                case 'leftarrow',
                    if Z > 1,
                        Z = Z - 1;
                    end
                case 'rightarrow'
                    if Z < UpperLim,
                        Z = Z + 1;
                    end
            end  
            setappdata(hFig, PlaneParamName, Z);
            %obj.DoShowZplaneInAStack(Z);
            WhatToShowFunc(Z, varargin{:});
        end
        
        function CompensationClickedCallback(obj,src, e)
            hFig = obj.hfig;
            allAxes = findobj(hFig,'type','axes','Tag','');
            % place checkboxes near axes
            if isappdata(src, 'CompChkBx'), % remove checkboxes
                hBx = getappdata(src, 'CompChkBx');
                delete(hBx);
                rmappdata(src, 'CompChkBx');
                set(gcf, 'Name', 'Press left or right arrow to navigate or Esc to leave')
            else
                for i = 1:length(allAxes),
                    pos =  get(allAxes(i), 'Position');
                    hBx(i) = uicontrol('Parent', hFig, 'Style','checkbox', ...
                        'String', 'Compensate','FontSize', 14, ...
                        'Units','normalized', ...
                        'Position', [pos(1), pos(2) - 0.05, 0.2, 0.05]);
                end
                setappdata(src, 'CompChkBx', hBx);
                set(gcf, 'Name', 'Select channels to compensate')
            end
        end
        
        function obj = DoShowZplaneInAStack(obj, Z, varargin)       
            P = getappdata(obj.hfig,'ChannelData');
%             for i = 1:length(P),
%                 CLim(i, :) = [min(P(i).Pic1(:)) max(P(i).Pic1(:))];
%             end

            
            for i = 1:length(P),
                ax(i) = subplot(2, ceil(length(P)/2), i);
                imagesc(P(i).Pic1(:, :, Z), P(i).CLim);
                colorbar
                axis image;
                title([P(i).ChannelNames{1} ', Z = ', num2str(Z)]);
            end
            linkaxes(ax);
           
        end
        
        function obj = DoGoOverZplanesNearCell(obj, varargin)
            Zhsc = ParseInputs('Zhsc', [], varargin);
            [hFig, ishFig] = ParseInputs('hFig', [], varargin);
            
            T = ParseInputs('Tstart', 1, varargin);
            if isempty(Zhsc),
                NoOfTimePoints = obj.NoTimePoints;
            else
                NoOfTimePoints = length(Zhsc);
            end
            if ishFig,
                obj.hfig = hFig;
                set(hFig, 'Position', get(0, 'Screensize'));
            else
                obj.hfig = figure('Position', get(0, 'Screensize'));
            end
            %T = 1;
            setappdata(obj.hfig,'T', T);
            gcf = obj.hfig;
            set(obj.hfig, 'KeyPressFcn', {@obj.keypress_callback, @obj.ShowZplanesNearCell, 'T', NoOfTimePoints, varargin{:}}, ...                
                'NumberTitle', 'off', 'Name', 'Press left or right arrow to navigate or Esc to leave');
            obj.ShowZplanesNearCell(T, varargin{:});
        end
        
        function obj = ShowZplanesNearCell(obj, T, varargin)       
            %P = getappdata(obj.hfig,'ChannelData');
            Color = ParseInputs('Color', 1, varargin);
            Zhsc = ParseInputs('Zhsc', [], varargin);
            Tstart = ParseInputs('Tstart', 1, varargin);
            [CellShapes, isShowCell] = ParseInputs('CellShapes', {}, varargin);
            CellShapeColor =  ParseInputs('CellShapeColor', 1, varargin);
            [CursorLine, isCursorLine] =  ParseInputs('CursorLine', {}, varargin);
            [NoSubplots, isNoSubplots] = ParseInputs('NumberOfSubplots', [], varargin);
            %patching for SCF - constant array
            if obj.NoTimePoints == 1,
                Ch = obj.DoGetSubStack('Color', Color,'Time', 1);
            else
                Ch = obj.DoGetSubStack('Color', Color,'Time', T);
            end
%             for i = 1:length(P),
%                 CLim(i, :) = [min(P(i).Pic1(:)) max(P(i).Pic1(:))];
%             end
            CellSize = ParseInputs('CellSize_um', 15, varargin);
            PixSizeZ = obj.PixSizeXYZ(3);
            % use correct property name 
            if isprop(obj, 'CellCoordImport'),
                CoordFieldName = 'CellCoordImport';
            elseif isprop(obj, 'PythonDataImport'),
                CoordFieldName = 'PythonDataImport';
            end
            
            if isempty(Zhsc),
                Zhsc = [obj.(CoordFieldName).Zhsc];
            end
            cellZ = Zhsc(T);
            ValidZ = find(abs((1:obj.NoZplanes)-cellZ)< (CellSize/(2*PixSizeZ)));
            if isShowCell,
                CellShape = CellShapes{T};
                tempImg = Ch.Pic1(:, :, ValidZ);
                maxPix = max(tempImg(:));
                Ch.Pic1(CellShape) = maxPix*CellShapeColor;
            end
            if T ~= Tstart,
                axLim = axis;
            end
 %           subplot(1,1,1)
            if isNoSubplots,
                NoSubplotRows = NoSubplots(1);
                NoSubplotCols = NoSubplots(2);
            else
                NoSubplotRows = 2;
                NoSubplotCols = ceil(length(ValidZ)/2);               
            end
            
            for i = 1:length(ValidZ),
                ax(i) = subplot(NoSubplotRows, NoSubplotCols, i);
                imagesc(Ch.Pic1(:, :, ValidZ(i)));%, P(i).CLim);
                colorbar
                axis image;
                title(['T = ', num2str(T) ', Z = ', num2str(ValidZ(i)), ...
                    ' Cell Z = ', num2str(cellZ)]);
                if T ~= Tstart, 
                    axis(axLim);
                end
            end
                           
            linkaxes(ax);
            if isCursorLine,
                if isa(CursorLine, 'matlab.graphics.chart.primitive.Line'),
                    set(CursorLine, 'XData', [T T])
                else
                    warning('This is not a valid Line handle!')
                end

            end
           
        end
        
        function CellParam = DoFindHSCAtTimeT(obj, T, varargin)
            CellSize = ParseInputs('CellSize_um', 10, varargin);
            Color = ParseInputs('Color', 1, varargin);
            rectSize = ParseInputs('Rect Size', 40, varargin);
            NoThresholds = ParseInputs('NoThresholds', 1, varargin);
            GradFilter = ParseInputs('GradientFilter', false, varargin);
            [VolumeToAdapt, isAdaptVolume] = ParseInputs('AdaptVolume_um3', NaN, varargin);
            VolumeToAdapt = VolumeToAdapt/prod(obj.PixSizeXYZ); % to voxels
            [ProjectedAreaToAdapt, isAdaptProjectedArea] = ParseInputs('AdaptProjectedArea_um2', NaN, varargin);
            ProjectedAreaToAdapt = ProjectedAreaToAdapt/prod(obj.PixSizeXYZ(1:2)); % to pixels
            [AreaToAdapt, isAdaptArea] = ParseInputs('AdaptArea_um2', NaN, varargin);
            AreaToAdapt = AreaToAdapt/prod(obj.PixSizeXYZ(1:2)); % to pixels
            % analyze background from the last SubtractBackgroundFramesNo
            % frames without a cell. If NaN do not subtract, if (-1)
            % subtract all
            SubtractBackgroundFramesNo = ParseInputs('Subtract Background', NaN, varargin);
            Run3D = ParseInputs('Run3D', true, varargin);
            
            Ch = obj.DoGetSubStack('Color', Color,'Time', T);
            Ch.Pic1 = double(Ch.Pic1);
            
            PixSizeZ = obj.PixSizeXYZ(3);
            if isprop(obj, 'CellCoordImport'),
                CoordFieldName = 'CellCoordImport';
            elseif isprop(obj, 'PythonDataImport'),
                CoordFieldName = 'PythonDataImport';
            end
            
            if iscolumn(obj.(CoordFieldName)(1).Xhsc)
                Xhsc = [];
                Yhsc = [];
                Zhsc = [];

                for i = 1:length(obj.(CoordFieldName)),
                    Xhsc = [Xhsc(:); obj.(CoordFieldName)(i).Xhsc(:)];
                    Yhsc = [Yhsc(:); obj.(CoordFieldName)(i).Yhsc];
                    Zhsc = [Zhsc(:); obj.(CoordFieldName)(i).Zhsc];
                end
            else    
                Xhsc = [obj.(CoordFieldName).Xhsc];
                Yhsc = [obj.(CoordFieldName).Yhsc];
                Zhsc = [obj.(CoordFieldName).Zhsc];
            end
            
            if T <= length(Xhsc),
                cellX = Xhsc(T);
                cellY = Yhsc(T);
                cellZ = Zhsc(T);
            else
                cellX = Xhsc(end);
                cellY = Yhsc(end);
                cellZ = Zhsc(end);
            end
            ValidZ = (abs((1:obj.NoZplanes)-cellZ)< (CellSize/(2*PixSizeZ)));
            img3D = Ch.Pic1(:, :, ValidZ);
            img3Dfin = img3D;
            img3Dfin(~isfinite(img3D)) = 0;
            img = sum(img3Dfin, 3);
%            img0 = img(round(cellY + ((-rectSize/2):(rectSize/2))), ...
%                            round(cellX + ((-rectSize/2):(rectSize/2))));
            H = fspecial('gaussian', [10 10], 5);
            se = strel('disk',2);
            
            img = imfilter(img, H);
            
            ROIx1 = max(round(cellX-rectSize/2), 1);
            ROIx2 = min(round(cellX+rectSize/2), size(img, 2));
            ROIy1 = max(round(cellY-rectSize/2), 1);
            ROIy2 = min(round(cellY+rectSize/2), size(img, 1));
%             img = img(round(cellY + ((-rectSize/2):(rectSize/2))), ...
%                             round(cellX + ((-rectSize/2):(rectSize/2))));
            img = img(ROIy1:ROIy2, ROIx1:ROIx2);
           
            if ~isnan(SubtractBackgroundFramesNo),
                % find frames with no cell
                bkgFramesXY = sqrt((Xhsc - cellX).^2 + (Yhsc - cellY).^2) > ...
                    rectSize/sqrt(2) + CellSize/(2*PixSizeZ);
                bkgFramesZ = abs(Zhsc - cellZ) > (CellSize/PixSizeZ); % could be half, but want to be sure
                bkgFrames = find(bkgFramesXY | bkgFramesZ);
                if SubtractBackgroundFramesNo > 0, % > -1 in fact
                    [~, II] = sort(abs(bkgFrames - T));
                    bkgFrames = bkgFrames(II(1:SubtractBackgroundFramesNo)); % closest SubtractBackgroundFramesNo frames
                end
                bkgCh = obj.DoGetSubStack('Color', Color,'Time', bkgFrames, 'Z', find(ValidZ));
                bkgPic = mean(double(bkgCh.Pic1), 3)*sum(ValidZ);
                bkgPic = imfilter(bkgPic, H);
%                 bkgPic = bkgPic(round(cellY + ((-rectSize/2):(rectSize/2))), ...
%                             round(cellX + ((-rectSize/2):(rectSize/2))));
                bkgPic = bkgPic(ROIy1:ROIy2, ROIx1:ROIx2);
            end
            
            if GradFilter,
                [Gmag, Gdir] = imgradient(img,'prewitt');
                thresh = multithresh(Gmag, NoThresholds);
                quantImg = imquantize(Gmag, thresh);
                boundSoma = bwboundaries(imopen((quantImg == NoThresholds+1), se), 'noholes');
            else
            %NoThresholds = 1;
                thresh = multithresh(img, NoThresholds);
                quantImg = imquantize(img, thresh);
                [boundSoma, boundInd] = bwboundaries(imopen((quantImg == NoThresholds+1), se), 'noholes');
            end
            if NoThresholds == 2,
                boundLam = bwboundaries(imopen(((quantImg == NoThresholds)| (quantImg == NoThresholds+1)), se), 'noholes');
            end
                         
            if isempty(boundSoma),
                CellParam.Zc = NaN;
                CellParam.widthZ_um = NaN;
                CellParam.brightnessZ = NaN;
                CellParam.backgroundZ = NaN;
                CellParam.contour.X = NaN;
                CellParam.contour.Y = NaN;
                CellParam.area = NaN;
                CellParam.Xc = NaN; %centroid
                CellParam.Yc = NaN; %centroid
                CellParam.perimeter = NaN;
                CellParam.D_um = NaN;
                CellParam.shapeFactor = NaN; 
                CellParam.Xcm = NaN;
                CellParam.Ycm = NaN;
                CellParam.polygeom.iner = [];
                CellParam.polygeom.cpmo = [];
                CellParam.Pix3D = [];
                return
            end
            [max_size, max_index] = max(cellfun('size', boundSoma, 1)); % find largest
            distFromCenter = cell2mat(cellfun(@mean, boundSoma, 'UniformOutput', false)) - repmat(size(img), length(boundSoma), 1)/2; % find closest to the center
            [~, closest_index] = min(sum(distFromCenter.^2, 2));
            brightest_ind = 0;
            brightest_mean =  mean(img(boundInd == 0));
            for bi = 1:length(boundSoma), % find brightest
                 newmean = mean(img(boundInd == bi));
                 if newmean > brightest_mean,
                     brightest_ind  = bi;
                     brightest_mean = newmean;
                 end
            end
            
            if ~((max_index == brightest_ind)&(max_index == closest_index)),
                if brightest_ind == closest_index,
                    max_index = closest_index;
                    warningMes = 'Warning: not the largest selection';
                elseif (max_index == brightest_ind),
                    warningMes = 'Warning: not the closest to center';
                elseif (max_index == closest_index),
                    warningMes = 'Warning: not the brightest selection';
                else 
                    max_index = closest_index;
                    warningMes = 'Warning: not the brightest or largest selection';
                end
            else
                warningMes = '';
            end
            
            boundSoma = boundSoma{max_index};
            % find boundLam containing boundSoma
            if NoThresholds == 2,
                for bL = 1:length(boundLam),
                    bLgood(bL) = any(inpolygon(boundSoma(:, 1), boundSoma(:, 2), boundLam{bL}(:, 1), boundLam{bL}(:, 2)));
                end
                bLgood = find(bLgood);
                if length(bLgood)~=1,
                    warning(['problem identifying soma at T = ' num2str(T)]);
                end
                boundSoma = boundLam{bLgood};
            end
                
            BW = poly2mask(boundSoma(:, 2), boundSoma(:, 1), size(img, 1), size(img, 2));
%             Xbw = sum(BW, 1)*(1:size(BW, 2))'/sum(BW(:))+round(cellX - rectSize/2)-1;
%             Ybw = (1:size(BW, 1))*sum(BW, 2)/sum(BW(:))+round(cellY - rectSize/2)-1;
            Xbw = sum(BW, 1)*(1:size(BW, 2))'/sum(BW(:))+ROIx1-1;
            Ybw = (1:size(BW, 1))*sum(BW, 2)/sum(BW(:))+ROIy1-1;


            CellArea = sum(BW(:));
            for z = 1:Ch.NoZplanes,
%                 planePic = double(Ch.Pic1(round(cellY + ((-rectSize/2):(rectSize/2))), ...
%                             round(cellX + ((-rectSize/2):(rectSize/2))), z));
                planePic = double(Ch.Pic1(ROIy1:ROIy2, ROIx1:ROIx2, z));

                BG = mean(mean(planePic.*(~BW)));
                CellBrightness(z) = sum(sum(planePic.*BW)) - BG*CellArea;
                if all(planePic(:) == 0),
                    CellBrightness(z) = NaN;
                end
                    
            end
            %prepare fit
            z = max(1, round(cellZ - CellSize/PixSizeZ)):min(Ch.NoZplanes, round(cellZ + CellSize/PixSizeZ));
            CB = CellBrightness(z);
            z = z(isfinite(CB));
            CB = CB(isfinite(CB));
            [mCB, mPos] = max(CB);
            A = mCB - min(CB);
            wSq = 3/10*(CellSize/PixSizeZ)^2;
            Zc = z(mPos);
            BL = min(CB);
            FitParam = nlinfitWeight2(z, CB, @GaussianBL, [A wSq Zc BL], [], []);
%             if (FitParam.beta(3) >= 1) & (FitParam.beta(3) <= obj.NoZplanes) & ... 
%                     (FitParam.beta(1) > FitParam.beta(4)),
            if (FitParam.beta(3) >= 1) & (FitParam.beta(3) <= obj.NoZplanes)

                CellParam.Zc = FitParam.beta(3);
            else
               CellParam.Zc = NaN;
            end
            
            CellParam.widthZ_um = sqrt(10*FitParam.beta(2)/3)*PixSizeZ;
            CellParam.brightnessZ = FitParam.beta(1);
            CellParam.backgroundZ = FitParam.beta(4);
            if Run3D,
                hSubplot = subplotplus({{[]},{{[]};{[]}}});
                view(3);
                camlight;
                lighting phong; 
                set(gcf,'CurrentAxes',hSubplot(1));
            else
                subplot(1, 2, 1)
            end
            
            plot(CellBrightness, '-o');
            hold all
            zz = (0:0.01:1)*(max(z) - min(z))+min(z);
            plot(z, CB,'o', zz, GaussianBL(FitParam.beta, zz), '-');
%             hold on
            plot([cellZ cellZ], get(gca, 'Ylim'), '--');
            plot([CellParam.Zc CellParam.Zc], get(gca, 'Ylim'), '-');
            hold off
            title(['T = ', num2str(T)]);
            
            
            % now recalculate XY shape having a better estimate of
            % parameters
            if ~isnan(CellParam.Zc)
                ValidZ = (abs((1:obj.NoZplanes) - CellParam.Zc) < (CellParam.widthZ_um/(2*PixSizeZ)));
                img3D = Ch.Pic1(:, :, ValidZ);
                img3Dfin = img3D;
                img3Dfin(~isfinite(img3D)) = 0;
                img = sum(img3Dfin, 3);

                %img = sum(Ch.Pic1(:, :, ValidZ), 3); 
                CellParam.image = img;
%                 img0 = img(round(Ybw + ((-rectSize/2):(rectSize/2))), ...
%                                 round(Xbw + ((-rectSize/2):(rectSize/2))));
                
                ROIx1bw = max(round(Xbw-rectSize/2), 1);
                ROIx2bw = min(round(Xbw+rectSize/2), size(img, 2));
                ROIy1bw = max(round(Ybw-rectSize/2), 1);
                ROIy2bw = min(round(Ybw+rectSize/2), size(img, 1));

                img0 = img(ROIy1bw:ROIy2bw, ROIx1bw:ROIx2bw);

                img = imfilter(img, H);
%                 img = img(round(Ybw + ((-rectSize/2):(rectSize/2))), ...
%                                 round(Xbw + ((-rectSize/2):(rectSize/2))));
                img = img(ROIy1bw:ROIy2bw, ROIx1bw:ROIx2bw);

                 % NoThresholds = 1;
                if Run3D,
                    for zz = 1:size(img3D, 3),
                        img3D(:, :, zz) = imfilter(img3D(:, :, zz), H);
                    end
%                     img3D = img3D(round(Ybw + ((-rectSize/2):(rectSize/2))), ...
%                                 round(Xbw + ((-rectSize/2):(rectSize/2))), :);
                    img3D = img3D(ROIy1bw:ROIy2bw, ROIx1bw:ROIx2bw, :);

                end

                if ~isnan(SubtractBackgroundFramesNo),
                    % find frames with no cell
                    bkgFramesXY = sqrt((Xhsc - Xbw).^2 + (Yhsc - Ybw).^2) > ...
                        rectSize/sqrt(2) + CellSize/(2*PixSizeZ);
                    bkgFramesZ = abs(Zhsc - CellParam.Zc) > (CellParam.widthZ_um/PixSizeZ); % could be half, but want to be sure
                    bkgFrames = find(bkgFramesXY | bkgFramesZ);
                    if SubtractBackgroundFramesNo > 0, % > -1 in fact
                        [~, II] = sort(abs(bkgFrames - T));
                        bkgFrames = bkgFrames(II(1:SubtractBackgroundFramesNo)); % closest SubtractBackgroundFramesNo frames
                    end
                    bkgCh = obj.DoGetSubStack('Color', Color,'Time', bkgFrames, 'Z', find(ValidZ));
                    bkgPic = mean(double(bkgCh.Pic1), 3)*sum(ValidZ);
                    CellParam.background = bkgPic;
                    bkgPic = imfilter(bkgPic, H);
%                     bkgPic = bkgPic(round(Ybw + ((-rectSize/2):(rectSize/2))), ...
%                                 round(Xbw + ((-rectSize/2):(rectSize/2))));
                    bkgPic = bkgPic(ROIy1bw:ROIy2bw, ROIx1bw:ROIx2bw);

                    img = img - bkgPic;
                    if Run3D,
                        z1 = 1;
                        if sum(ValidZ),
                            for zz = find(ValidZ),
                                bkgCh = obj.DoGetSubStack('Color', Color,'Time', bkgFrames, 'Z', zz);
                                tempPic = mean(double(bkgCh.Pic1), 3);
                                tempPic = imfilter(tempPic, H);
%                                 bkgPic3D(:, :, z1) = tempPic(round(Ybw + ((-rectSize/2):(rectSize/2))), ...
%                                     round(Xbw + ((-rectSize/2):(rectSize/2))));
                                bkgPic3D(:, :, z1) = tempPic(ROIy1bw:ROIy2bw, ROIx1bw:ROIx2bw);

                                z1 = z1 + 1;
                            end
                            img3D = img3D - bkgPic3D;
                        end
                    end
                    
                end
                
                if GradFilter,
                    [Gmag, Gdir] = imgradient(img,'prewitt');
                    imgForThreshold = Gmag;
%                     thresh = multithresh(Gmag, NoThresholds);
%                     quantImg = imquantize(Gmag, thresh);
                else
                    imgForThreshold = img;
%                     thresh = multithresh(img, NoThresholds);
%                     quantImg = imquantize(img, thresh);                   
                end
                thresh = multithresh(imgForThreshold, NoThresholds);
                

%                 [CellParam, boundSoma, warningMes] = obj.Apply2Dthreshold(CellParam, imgForThreshold, thresh, se, rectSize, Xbw, Ybw);
                [CellParam, boundSoma, warningMes] = obj.Apply2Dthreshold(CellParam, imgForThreshold, thresh, se, ROIx1bw, ROIy1bw);
                if isAdaptArea,
                    increment = 1.02;
                    currentValue = CellParam.area;
                    if (currentValue < AreaToAdapt),
                        while (currentValue < AreaToAdapt),
                            lowBound = currentValue;
                            lowThresh = thresh;

                            thresh = thresh/increment;
%                             CellParam = obj.Apply2Dthreshold(CellParam, imgForThreshold, thresh, se, rectSize, Xbw, Ybw);  
                            CellParam = obj.Apply2Dthreshold(CellParam, imgForThreshold, thresh, se, ROIx1bw, ROIy1bw);                     

                            currentValue = CellParam.area;
                        end
                        highBound = currentValue;
                        highThresh = thresh;
                        thresh = interp1([lowBound highBound], [lowThresh highThresh], ...
                            AreaToAdapt, 'linear','extrap');  

                    elseif (currentValue > AreaToAdapt)
                        while (currentValue > AreaToAdapt),
                            highBound = currentValue;
                            highThresh = thresh;

                            thresh = thresh*increment;
%                             CellParam = obj.Apply2Dthreshold(CellParam, imgForThreshold, thresh, se, rectSize, Xbw, Ybw);
                            CellParam = obj.Apply2Dthreshold(CellParam, imgForThreshold, thresh, se, ROIx1bw, ROIy1bw);                     

                            currentValue = CellParam.area;
                        end
                        lowBound = currentValue;
                        lowThresh = thresh;
                        thresh = interp1([lowBound highBound], [lowThresh highThresh], ...
                            AreaToAdapt, 'linear','extrap');  

                    end
%                    [CellParam, boundSoma, warningMes] = obj.Apply2Dthreshold(CellParam, imgForThreshold, thresh, se, rectSize, Xbw, Ybw);
                    [CellParam, boundSoma, warningMes] = obj.Apply2Dthreshold(CellParam, imgForThreshold, thresh, se, ROIx1bw, ROIy1bw);
                end
                %CellParam.image = sum(Ch.Pic1(:, :, ValidZ), 3);
                if ~isempty(warningMes),
                    disp(warningMes);
                end
                
                if Run3D,
                    thr3D = multithresh(img3D, 1);
                    [PixIdx, ii3, jj3, kk3, Vol, bound3Dprojection, Aproj, iShell, jShell, kShell] = ...
                                obj.Apply3Dthreshold(img3D, thr3D);
                    if isAdaptVolume | isAdaptProjectedArea,
                        increment = 1.02;
                        if isAdaptVolume,
                            valueToAdapt = VolumeToAdapt;
                            currentValue = Vol;
                        else
                            valueToAdapt = ProjectedAreaToAdapt;
                            currentValue = Aproj;
                        end
                        
                        if (currentValue < valueToAdapt),
                            while (currentValue < valueToAdapt),
                                lowBound = currentValue;
                                lowThresh = thr3D;

                                thr3D = thr3D/increment;
                                [~, ~, ~, ~, Vol, ~, Aproj] = ...
                                    obj.Apply3Dthreshold(img3D, thr3D);
                                if isAdaptVolume,
                                    currentValue = Vol;
                                else
                                    currentValue = Aproj;
                                end
                            end
                            highBound = currentValue;
                            highThresh = thr3D;
                            thr3D = interp1([lowBound highBound], [lowThresh highThresh], ...
                                valueToAdapt, 'linear','extrap');  

                        elseif (currentValue > valueToAdapt)
                            while (currentValue > valueToAdapt),
                                highBound = currentValue;
                                highThresh = thr3D;

                                thr3D = thr3D*increment;
                                [~, ~, ~, ~, Vol, ~, Aproj] = ...
                                    obj.Apply3Dthreshold(img3D, thr3D);
                                if isAdaptVolume,
                                    currentValue = Vol;
                                else
                                    currentValue = Aproj;
                                end
                            end
                            lowBound = currentValue;
                            lowThresh = thr3D;
                            thr3D = interp1([lowBound highBound], [lowThresh highThresh], ...
                                valueToAdapt, 'linear','extrap');  
                        
                        end
                        [PixIdx, ii3, jj3, kk3, Vol, bound3Dprojection, Aproj, iShell, jShell, kShell] = ...
                                obj.Apply3Dthreshold(img3D, thr3D);
                    end
                        
%                     ii3 = ii3 + round(Ybw - rectSize/2)-1;
%                     jj3 = jj3 + round(Xbw - rectSize/2)-1;
                    ii3 = ii3 + ROIy1bw-1;
                    jj3 = jj3 + ROIx1bw-1;

                    kk3 = kk3 + find(ValidZ, 1) - 1;
                    Pix3D = sub2ind(size(Ch.Pic1), ii3, jj3, kk3);
                    CellParam.Pix3D = Pix3D;
%                     iShell = iShell + round(Ybw - rectSize/2)-1;
%                     jShell = jShell + round(Xbw - rectSize/2)-1;
                    iShell = iShell + ROIy1bw -1;
                    jShell = jShell + ROIx1bw -1;

                    kShell = kShell + find(ValidZ, 1) - 1;
                    Pix3Dshell = sub2ind(size(Ch.Pic1), iShell, jShell, kShell);
                    CellParam.Pix3Dshell = Pix3Dshell;
                    
%                     CellParam.projection(:, 1) = bound3Dprojection(:, 1)+round(Ybw - rectSize/2)-1;
%                     CellParam.projection(:, 2) = bound3Dprojection(:, 2)+round(Xbw - rectSize/2)-1;
                    CellParam.projection(:, 1) = bound3Dprojection(:, 1)+ROIy1bw-1;
                    CellParam.projection(:, 2) = bound3Dprojection(:, 2)+ROIx1bw-1;

                    CellParam.Volume_um3 = Vol*prod(obj.PixSizeXYZ);
                    CellParam.projArea_um2 = Aproj*prod(obj.PixSizeXYZ(1:2));
                end
                
                if Run3D,
                    set(gcf,'CurrentAxes',hSubplot(2));
                else
                    subplot(1, 2, 2)
                end
                %imagesc(img0);
                imagesc(img);
                axis image;
                hold on
                plot(boundSoma(:,2), boundSoma(:,1), '-r', 'LineWidth', 2);
                if Run3D,
                    plot(bound3Dprojection(:, 2), bound3Dprojection(:, 1), '-w','LineWidth', 2); 
                end
%                 arrow([CellParam.Xcm - round(Xbw - rectSize/2)+1, ...
%                     CellParam.Ycm - round(Ybw - rectSize/2)+1], ...
%                     [CellParam.Xc - round(Xbw - rectSize/2)+1, ...
%                     CellParam.Yc - round(Ybw - rectSize/2)+1], 'Color', 'k')
                arrow([CellParam.Xcm - ROIx1bw+1, ...
                    CellParam.Ycm - ROIy1bw+1], ...
                    [CellParam.Xc - ROIx1bw+1, ...
                    CellParam.Yc - ROIy1bw+1], 'Color', 'k')

                hold off
            
                if Run3D,
                    set(gcf,'CurrentAxes',hSubplot(3));
                    img3D = zeros(size(Ch.Pic1));
                    img3D(Pix3D) = 1;
                    isosurface(img3D, 0.5)
                end


            end
            figure(gcf)
            
        end
        
        function obj = DoFindHSCparam(obj, varargin)
            [~, isNoPause] = ParseInputs('NoPause', NaN, varargin);
            Trange = ParseInputs('Trange',  1:obj.NoTimePoints, varargin);

            for T = Trange,
                T
                CP = obj.DoFindHSCAtTimeT(T, varargin{:});
                if ~isNoPause,
                    pause;
                end
                fldnames = fieldnames(CP);
                for i = 1:length(fldnames),
                    if strcmp(fldnames{i}, 'image') | strcmp(fldnames{i}, 'background'),
                        CellParam.(fldnames{i})(:, :, T) = CP.(fldnames{i});
                    elseif strcmp(fldnames{i}, 'Pix3D'),
                        CellParam.Pix3D{T} = CP.Pix3D; 
                    elseif strcmp(fldnames{i}, 'Pix3Dshell'),
                        CellParam.Pix3Dshell{T} = CP.Pix3Dshell; 
                    elseif strcmp(fldnames{i}, 'projection'),
                        CellParam.projection{T} = CP.projection; 
                    elseif ~strcmp(fldnames{i}, 'imageRef'),
                        CellParam.(fldnames{i})(T) = CP.(fldnames{i});
                    end
                end
            end
            
            if ~isprop(obj, 'CellParam'),
                obj.addprop('CellParam');
            end
            obj.CellParam = CellParam;
        end

        function obj = DoShowCellXYmotion(obj, varargin)
            ShowContour = ParseInputs('ShowContour', true, varargin);
            ImageField = ParseInputs('ImageField', 'image', varargin);
            ShowWatershed = ParseInputs('ShowWatershed', false, varargin);
            ShowDirectStretches = ParseInputs('ShowDirectStretches', false, varargin);
            ShowPolarity = ParseInputs('ShowPolarity', false, varargin);
            MovieName = ParseInputs('CreateMovie', '', varargin);
            if isprop(obj, 'CellParam'),
                if ~isempty(MovieName),
                    vidObj = VideoWriter(MovieName);
                    set(vidObj, 'FrameRate', 5);
                    open(vidObj);
                end
                CP = obj.CellParam;
                subplot(1, 1, 1);
                for T = 1:size(CP.image, 3),
                    if iscell(ImageField),
                        for imgFieldNo = 1:length(ImageField)
                            img(:, :, imgFieldNo) = CP.(ImageField{imgFieldNo})(:, :, T);
                            img(:, :, imgFieldNo) = img(:, :, imgFieldNo)/max(max(img(:, :, imgFieldNo)));
                        end
                        if size(img, 3) < 3,
                            img(:, :, end:3) = 0;
                        end
                        
                    else
                        img = CP.(ImageField)(:, :, T);
                    end
                    if ShowWatershed,
                        L = CP.imageRefWatershed(:, :, T);
                        img(L == 0) = min(img(:));
                    end
                   
                    if T > 1,
                        v = axis;
                        imagesc(img);
                        axis(v);
                        set(gca, 'DataAspectRatio', [1 1 1])
                    else
                        imagesc(img);
                        axis image;
                    end
                    title(['T = ',  num2str(T), ' Z = ', num2str(CP.Zc(T))])
                    if ShowContour,
                        hold on
                        plot(CP.contour(T).X, CP.contour(T).Y, '-r', 'LineWidth', 2);
                        if isfield(CP, 'projection'),
                            if ~isempty(CP.projection{T}),
                                plot(CP.projection{T}(:, 2), CP.projection{T}(:, 1), '-w', 'LineWidth', 2);
                            end
                        end
                        hold off
                    end
                    
                    if ShowDirectStretches,
                        hold on
                        if isprop(obj, 'DirectedStretch'),
                            for i = 1:length(obj.DirectedStretch),
                                plot(obj.DirectedStretch(i).Xhsc, obj.DirectedStretch(i).Yhsc, '-w');
                            end
                        elseif isprop(obj, 'PythonDataImport'),
                            DS = obj.PythonDataImport.DirectStretchList;
                            for i = 1:length(DS),
                                plot(DS{i}.Xhsc, DS{i}.Yhsc, '-w');
                            end
                        end
                        hold off
                    end
                    
                    if ShowPolarity,
%                         if exist('arrowH'),
%                             delete(arrowH);
%                         end
                        hold on
%                         set(gcf, 'Units', 'pixels');
%                         pos = get(gcf, 'position');
                                          
                        arrowH = arrow([CP.Xcm(T), CP.Ycm(T)], ...
                            [CP.Xcm(T), CP.Ycm(T)] + 30*([CP.Xc(T)-CP.Xcm(T), CP.Yc(T)-CP.Ycm(T)]), ...
                            'Color', 'r')
%                        set(gcf, 'Units', 'normalized');
                        hold off
                    end
                    
                    figure(gcf)
                    if ~isempty(MovieName)
                        axis off
                        %[Xfig Yfig] = AxesToFigureCoordinates([160 180], [20 40])
                        if exist('hAn'),
                            delete(hAn);
                        end
                        hAn = annotation('textbox', [0.68, 0.8, 0.1, 0.1], 'String', num2str(T), 'FontSize', 20, 'Color', [1 1 1], 'LineStyle', 'none');  % time of drug application

                        currFrame = getframe;
                        writeVideo(vidObj,currFrame);
                    end
                    pause
                end
                
                if ~isempty(MovieName),
                    close(vidObj);
                end
            else
                error('Run DoFindHSCparam first!');
            end
        end
        
        function TimeAveragedStack = DoGetTimeAveragedStack(obj, varargin)
            Color = ParseInputs('Color', 1, varargin);
            T = ParseInputs('Time', 1:obj.NoTimePoints, varargin);
            GFPave = [];
            for Z = 1:obj.NoZplanes,
        %        disp([T, Z])
                subStack = obj.DoGetSubStack('Z', Z, 'Color', Color, 'Time', T);
                PicAve(:, :, Z) = mean(double(subStack.Pic1), 3);
            end
            TimeAveragedStack = obj.DoGetSubStack('Color', Color, 'Time', T(1));
            TimeAveragedStack.Pic1 = PicAve;
        end
        
        function obj = DoGetHSCplanesFromOtherStack(obj, gfp, varargin)
            AddWatershed = ParseInputs('AddWatershed', true, varargin);
            WeightType = ParseInputs('WeightType', 'sphere', varargin);
            CellSize_um = ParseInputs('CellSize_um', 10, varargin);
            if isprop(obj, 'CellParam'),
                CP = obj.CellParam;
                Zplanes = 1:obj.NoZplanes;
                for T = 1:size(CP.image, 3),
                    dd = Zplanes - CP.Zc(T);
                    %CellDia_plns = CP.widthZ_um(T)/obj.PixSizeXYZ(3);
                    %CellDia_plns = median(CP.D_um)/obj.PixSizeXYZ(3);
                    %CellDia_plns = median(CP.widthZ_um(~isnan(CP.widthZ_um)))/obj.PixSizeXYZ(3);
                    CellDia_plns = CellSize_um/obj.PixSizeXYZ(3);
                    if strcmp(WeightType, 'sphere'),
                        weight = pi*sqrt(CellDia_plns^2 - 4*dd.^2).*(abs(dd) < (CellDia_plns/2));
                        
                    elseif strcmp(WeightType, 'step'),
                        weight  = abs(dd) < (CellDia_plns/2);
                    end
                    weight = weight/sum(weight);
                    weightedImage = reshape(reshape(double(gfp.Pic1), [], gfp.NoZplanes)*weight(:), size(gfp.Pic1, 1), size(gfp.Pic1, 2));
                    CP.imageRef(:, :, T) = weightedImage;
                    if AddWatershed,
                        if ~any(isnan(weightedImage)),
                            ggfp = imfilter(weightedImage, fspecial('disk', 5));
                            CP.imageRefWatershed(:, :, T) = watershed(ggfp);
                        end
                    end
                end
                obj.CellParam = CP;
            else
                error('Run DoFindHSCparam first!');
            end
        end
        
        function WSstruc = DoGetOnWaterShedFractionsAtTimeT(obj, T, varargin)
            Ntrials = ParseInputs('No of Trials', 100, varargin);
            OnFracDistance = ParseInputs('OnFractDistance', 1, varargin);
            [WSdil, isWSdil] = ParseInputs('WSdilated', [], varargin);
            Dia = obj.CellParam.D_um(T)/obj.PixSizeXYZ(1);
            if ~isWSdil,
                WS = double(obj.CellParam.imageRefWatershed(:, :, T) == 0);
            end
            contr = obj.CellParam.contour(T);
            if isempty(contr.X),
                WSstruc.dist = NaN;
                WSstruc.onFrac = NaN;
                return
            end
            
            for i = 1 : (Ntrials+1),
                if i == 1,
                    Xshift = 0;
                    Yshift = 0;
                else
                    Xshift  = round(Dia*(2*rand - 1));
                    Yshift  = round(Dia*(2*rand - 1)); 
                end
                clear dist
                if ~isWSdil, 
                    ind = sub2ind(size(WS), contr.Y + Yshift, contr.X + Xshift);
                    dist(1) = sum(WS(ind));
                    R = 1;
                    while dist(end) < length(contr.Y),
                        WSdil = imdilate(WS, strel('disk', R, 0));
                        dist(R+1) = sum(WSdil(ind));
                        R = R + 1;
                        WSstruc(i).dist = dist; % diff([0 dist]);
                        %[~, ii] = find(Psim(i).dist > (Psim(i).dist(end)/3), 1);
                        %Psim(i).medianDist = ii-1;
                    end
                    WSstruc(i).onFrac = WSstruc(i).dist(OnFracDistance+1)/WSstruc(i).dist(end);
                else % do 3D analysis
                    Pix3Dshell = obj.CellParam.Pix3Dshell{T};
                    [ii, jj, kk] = ind2sub(size(WSdil{1}), Pix3Dshell);
                    newSub = [ii, jj, kk]  + repmat([Yshift, Xshift, 0], size(Pix3Dshell));
                    ind = sub2ind(size(WSdil{1}), newSub(:, 1), newSub(:, 2), newSub(:, 3));
                    for R = 1:length(WSdil),
                        dist(R) = sum(WSdil{R}(ind));
                        WSstruc(i).dist = dist; % diff([0 dist]);
                    end
                    WSstruc(i).onFrac = WSstruc(i).dist(OnFracDistance+1)/length(Pix3Dshell);
                end
            end
            %disp([Psim.medianDist])
        end
        
        function WSobj = DoApplyWaterShedTransform(obj, varargin)
            FilterSize_um = ParseInputs('FilterSize_um', 2.5, varargin);
            FilterSize = round(FilterSize_um/obj.PixSizeXYZ(1));
            ggfp = imfilter(obj.Pic1, fspecial('disk', FilterSize));
            WSimage = zeros(size(obj.Pic1));
            for k = 1:size(obj.Pic1, 3),
                WSimage(:, :, k) = (watershed(ggfp(:, :, k)) == 0);
            end
            WSobj = obj.DoCloneObject;
            WSobj.Pic1 = WSimage;
        end
            
        function obj = DoGetOnWaterShedFractions(obj, varargin)
            [WSImage, isWSImage] = ParseInputs('WaterShedImage', [], varargin);
            [WSobj, isWSobj] = ParseInputs('WaterShedObject', [], varargin);
            if isWSobj,
                WSdil{1} = WSobj.Pic1;
                for R = 1:6,
                    WSdil{R+1} = imdilate(WSdil{1}, strel('disk', R, 0));
                end
            end
 
            for T = 1:length(obj.CellParam.contour),
              %  disp(T)
                if isWSobj,
                    WSstruc = obj.DoGetOnWaterShedFractionsAtTimeT(T, 'WSdilated', WSdil, varargin{:});
                else
                    WSstruc = obj.DoGetOnWaterShedFractionsAtTimeT(T, varargin{:});
                end
                if isnan(WSstruc(1).onFrac),
                    WSstruc(1).percentile = NaN;
                else
                    WSstruc(1).percentile = ...
                        (sum([WSstruc.onFrac] >= WSstruc(1).onFrac)-1)/(length(WSstruc)-1);
                end

                obj.CellParam.WaterShedDistances(T) = WSstruc(1);
                obj.CellParam.WaterShedDistSimulated{T} = WSstruc(2:end);
            end
        end
        
        function [NeighborStruc, Pic1quant] = DoGetCellNeighborhoodAtTimeT(obj, MapObj, T, varargin)
            shellNos = ParseInputs('Shell No', [-2, 4], varargin);
            DoZshell = ParseInputs('DoZshell', true, varargin);
            Ntrials = ParseInputs('No of Trials', 100, varargin);
            bgMethod = ParseInputs('Background Method', 'local', varargin);
            percentiles =  ParseInputs('Percentiles', 10:5:95, varargin);
            Dia = obj.CellParam.D_um(T)/obj.PixSizeXYZ(1);
            MaxDispDia = ParseInputs('MaxDisplacement in dia', 1, varargin);
            MaxDisp = MaxDispDia*Dia;
            pix3D = obj.CellParam.Pix3D{T};
            if isempty(pix3D),
                NeighborStruc.shell = NaN;
                NeighborStruc.pix = NaN;
                NeighborStruc.Xshift = [];
                NeighborStruc.Yshift = [];
                NeighborStruc.pixVal = NaN;              
                return
            end
          
            Pic1quant = zeros(size(MapObj.Pic1));
            for planeNo = 1:size(MapObj.Pic1, 3),
                tempPlane = double(MapObj.Pic1(:,:, planeNo));
                if all(tempPlane(:) == mean(tempPlane(:))),
                    Pic1quant(:, :, planeNo) = NaN;
                else
                    prc(planeNo, :) = prctile(tempPlane(:), percentiles);
                    Pic1quant(:, :, planeNo) = imquantize(tempPlane,  prc(planeNo, :));
                end
            end

            img3D = zeros([obj.stackSizeY obj.stackSizeX obj.NoZplanes]);
            img3D(pix3D) = 1;
            %startPix = MapObj.Pic1(obj.CellParam.Pix3D{T});
            %check average outside the cell
            if strcmp(bgMethod, 'local'),
                H = strel('disk', round(MaxDisp), 0);
                bkgImg = zeros(size(img3D));
                outImg = ((imdilate(img3D, H) - img3D)>0);
                bkgImg(outImg) = MapObj.Pic1(outImg);
                meanBkgZ = double(sum(sum(bkgImg, 1), 2))./...
                    double(sum(sum(outImg, 1), 2));
                meanBkgZ = squeeze(squeeze(meanBkgZ));
                stdBkgZ = [];
                for k = 1:size(bkgImg, 3),
                    bbkg = bkgImg(:, :, k);
                    bbkg = bbkg(outImg(:, :, k));
                    stdBkgZ(k) = std(bbkg);
                end
            else % background global
                meanBkgZ = median(reshape(MapObj.Pic1, prod(size(MapObj.Pic1(:, :, 1))), size(MapObj.Pic1, 3)))';
                stdBkgZ = mad(reshape(MapObj.Pic1, prod(size(MapObj.Pic1(:, :, 1))), size(MapObj.Pic1, 3)))';
            end

            Xshift  = [0, round(Dia*(2*rand(1, Ntrials) - 1))];
            Yshift  = [0, round(Dia*(2*rand(1, Ntrials) - 1))]; 

            NeighborStruc(1).shell = 0;
            NeighborStruc(1).pix = obj.CellParam.Pix3D{T};
            NeighborStruc(1).Npix = length(NeighborStruc(1).pix);
            NeighborStruc(1).pixVal = MapObj.Pic1(NeighborStruc(1).pix);
            NeighborStruc(1).meanBkgZ = meanBkgZ;
            NeighborStruc(1).stdBkgZ = stdBkgZ;
            [ii3, jj3, kk3] = ind2sub(size(img3D), NeighborStruc(1).pix);
            NeighborStruc(1).relPixVal = MapObj.Pic1(NeighborStruc(1).pix)./...
                meanBkgZ(kk3);
            NeighborStruc(1).Xshift = Xshift;
            NeighborStruc(1).Yshift = Yshift;
            NeighborStruc(1).percentiles = percentiles;
            NeighborStruc(1).prcntlsValues = prc;
           % NeighborStruc(1).MapImgQuantized = Pic1quant;
           
            for iShift = 1:length(Xshift),
                newSub = [ii3, jj3, kk3]  + repmat([Yshift(iShift), Xshift(iShift), 0], size(ii3));
                try 
                    ind = sub2ind(size(img3D), newSub(:, 1), newSub(:, 2), newSub(:, 3));
                    HH = hist(Pic1quant(ind), 1:(length(percentiles)+1));
                    NeighborStruc(1).quantizedHistograms(:, iShift) = HH(:);
                    NeighborStruc(1).meanShell(iShift) = mean(MapObj.Pic1(ind));
                catch
                    NeighborStruc(1).meanShell(iShift) = NaN;
                end
            end



%             NeighborStruc(1).relTostdPixVal = MapObj.Pic1(NeighborStruc(1).pix)./...
%                 stdBkgZ(kk3)';

            
            H = strel('disk', 1, 0);
            %split shells into those to erode and those to dilate
            shellNeg = (-1):(-1):shellNos(1);
            shellPos = 1:shellNos(end);
            prevImg = img3D;
            for shell = shellNeg,
                nextImg = imerode(prevImg, H);
                NeighborStruc(end+1).shell = shell;
                NeighborStruc(end).pix = find(prevImg - nextImg);
                NeighborStruc(end).Npix = length(NeighborStruc(end).pix);
                NeighborStruc(end).pixVal = MapObj.Pic1(NeighborStruc(end).pix);
                NeighborStruc(end).meanBkgZ = meanBkgZ;
                NeighborStruc(end).stdBkgZ = stdBkgZ;
                [ii3, jj3, kk3] = ind2sub(size(img3D), NeighborStruc(end).pix);
                NeighborStruc(end).relPixVal = MapObj.Pic1(NeighborStruc(end).pix)./...
                    meanBkgZ(kk3);                        
                for iShift = 1:length(Xshift),
                    newSub = [ii3, jj3, kk3]  + repmat([Yshift(iShift), Xshift(iShift), 0], size(ii3));
                    try 
                        ind = sub2ind(size(img3D), newSub(:, 1), newSub(:, 2), newSub(:, 3));
                        HH = hist(Pic1quant(ind), 1:(length(percentiles)+1));
                        NeighborStruc(end).quantizedHistograms(:, iShift) = HH(:);
                        NeighborStruc(end).meanShell(iShift) = mean(MapObj.Pic1(ind));
                    catch
                        NeighborStruc(end).meanShell(iShift) = NaN;
                    end
                end
    %                 NeighborStruc(end).relTostdPixVal = MapObj.Pic1(NeighborStruc(end).pix)./...
%                 stdBkgZ(kk3)';
                prevImg = nextImg;
            end
            
            prevImg = img3D;
            for shell = shellPos,
                nextImg = imdilate(prevImg, H);
                NeighborStruc(end+1).shell = shell;
                NeighborStruc(end).pix = find(nextImg - prevImg);
                NeighborStruc(end).Npix = length(NeighborStruc(end).pix);
                NeighborStruc(end).pixVal = MapObj.Pic1(NeighborStruc(end).pix);
                NeighborStruc(end).meanBkgZ = meanBkgZ;
                NeighborStruc(end).stdBkgZ = stdBkgZ;
                [ii3, jj3, kk3] = ind2sub(size(img3D), NeighborStruc(end).pix);
                NeighborStruc(end).relPixVal = MapObj.Pic1(NeighborStruc(end).pix)./...
                    meanBkgZ(kk3);
                for iShift = 1:length(Xshift),
                    newSub = [ii3, jj3, kk3]  + repmat([Yshift(iShift), Xshift(iShift), 0], size(ii3));
                    try
                        ind = sub2ind(size(img3D), newSub(:, 1), newSub(:, 2), newSub(:, 3));
                        HH = hist(Pic1quant(ind), 1:(length(percentiles)+1));
                        NeighborStruc(end).quantizedHistograms(:, iShift) = HH(:);
                        NeighborStruc(end).meanShell(iShift) = mean(MapObj.Pic1(ind));
                    catch
                        NeighborStruc(end).meanShell(iShift) = NaN;
                    end
                end
 
%                 NeighborStruc(end).relTostdPixVal = (MapObj.Pic1(NeighborStruc(end).pix))./...
%                 stdBkgZ(kk3)';
                prevImg = nextImg;
            end
            
            if DoZshell,
                imgTemp = zeros(size(img3D));
                imgTemp(:, :, 2:end) = img3D(:, :, 1:(end-1));
                Zshell1 =  find((imgTemp - img3D)> 0.5);
                
                imgTemp = zeros(size(img3D));
                imgTemp(:, :, 1:(end -1)) = img3D(:, :, 2:end);
                Zshell2 = find((imgTemp - img3D)> 0.5);
                
                Zshell = unique([Zshell1; Zshell2]);
                
                NeighborStruc(1).ZshellPix = Zshell;
                NeighborStruc(1).NoZshellPix = length(Zshell);
                [ii3, jj3, kk3] = ind2sub(size(img3D), Zshell);
                
                for iShift = 1:length(Xshift),
                    newSub = [ii3, jj3, kk3]  + repmat([Yshift(iShift), Xshift(iShift), 0], size(ii3));
                    try
                        ind = sub2ind(size(img3D), newSub(:, 1), newSub(:, 2), newSub(:, 3));
                        NeighborStruc(1).meanZShell(iShift) = mean(MapObj.Pic1(ind));
                    catch
                        NeighborStruc(1).meanZShell(iShift) = NaN;
                    end
                end

                
            end
            
%             Xshift  = [0, round(MaxDisp*(2*rand(1, Ntrials) - 1))];
%             Yshift  = [0, round(MaxDisp*(2*rand(1, Ntrials) - 1))]; 
% 
%             for i = 1:length(NeighborStruc),
%                 [ii3, jj3, kk3] = ind2sub(size(img3D), NeighborStruc(i).pix);
%                 for j = 1 : length(Xshift),
%                     newPix = sub2ind(size(img3D), ...
%                         ii3+Yshift(j), jj3+Xshift(j), kk3);
%                     NeighborStruc(i).pixVal{j} = MapObj.Pic1(newPix);
%                 end
%                 NeighborStruc(i).Xshift = Xshift;
%                 NeighborStruc(i).Yshift = Yshift;
%             end
        end
        
        function obj = DoGetCellNeighborhood(obj, MapObj, varargin)
            TimeRange = ParseInputs('Time', 1:length(obj.CellParam.contour), varargin);
            for T = TimeRange,
                disp(T)
                NS = obj.DoGetCellNeighborhoodAtTimeT(MapObj, T, varargin{:});
                obj.CellParam.Neighborhood{T} = NS;
            end
        end
        
        function obj = DoApplyFilterToStack(obj, varargin)
            FilterType = ParseInputs('FilterType', 'Gaussian3D', varargin);
            Rg = 3; % um: corresponds roughly to cell diameter of 8 um
            FilterWidth = ParseInputs('FilterWidth um', Rg, varargin);
            
            if strcmp(FilterType, 'Gaussian3D'),
                % take filter half size to be at least 2Rg
                GaussFiltSizeXYZ = ceil(2*[Rg; Rg; Rg]./obj.PixSizeXYZ);
                [xx,yy,zz] = ndgrid(-GaussFiltSizeXYZ(1):GaussFiltSizeXYZ(1), ...
                    -GaussFiltSizeXYZ(2):GaussFiltSizeXYZ(2), -GaussFiltSizeXYZ(3):GaussFiltSizeXYZ(3));
                GaussFilter = exp(-3/2*(xx.^2*obj.PixSizeXYZ(1)^2/Rg^2 + yy.^2*obj.PixSizeXYZ(2)^2/Rg^2 + ...
                    zz.^2*obj.PixSizeXYZ(3)^2/Rg^2));  
                FLT = GaussFilter/sum(GaussFilter(:));
            end
            % apply filter - remove mean before the application of Gaussian
            % filter
            
            % Test how many channels and time points are involved
            for T = 1:obj.NoTimePoints,
                for C = 1:obj.NoColors,
                    disp(['running T = ' num2str(T) ' point out of ' num2str(obj.NoTimePoints) ...
                        ' Color = ' num2str(C) ' point out of ' num2str(obj.NoColors)]);
                    P = obj.DoGetSubStack('Color', C, 'Time', T);           
                    PicFLT = imfilter(P.Pic1 - mean(P.Pic1(:)), FLT);
                    planeNo = obj.CZTtoNoInStack('Color', C, 'Time', T);
                    obj.Pic1(:, :, planeNo) = PicFLT;
                end
            end            
        end
        
        function obj = DoFindCells(obj, varargin)
             if isprop(obj, 'Cells'),
                 delete(obj.dynprop_meta.Cells);
                 obj.dynprop_meta = rmfield(obj.dynprop_meta, 'Cells');
             end
             CellDiameter = ParseInputs('Cell size (um)', 8, varargin);
             CellMinDia = ParseInputs('Min Cell size (um)', 2, varargin);
             MaxCellNo = ParseInputs('Max Cell No', 6, varargin);
             SelectField = ParseInputs('SelectField', 'maxPixValue', varargin);
             [~, IsFilter] = ParseInputs('ApplyFilter', [], varargin);
             
                % recalculate into Gaussian width through gyration radius
             Rg = CellDiameter/2*sqrt(3/5);
                % take filter half size to be at least 2Rg
             if IsFilter,
                 GaussFiltSizeXYZ = ceil(2*[Rg; Rg; Rg]./obj.PixSizeXYZ);
                 [xx,yy,zz] = ndgrid(-GaussFiltSizeXYZ(1):GaussFiltSizeXYZ(1), ...
                        -GaussFiltSizeXYZ(2):GaussFiltSizeXYZ(2), -GaussFiltSizeXYZ(3):GaussFiltSizeXYZ(3));
                 GaussFilter = exp(-3/2*(xx.^2*obj.PixSizeXYZ(1)^2/Rg^2 + yy.^2*obj.PixSizeXYZ(2)^2/Rg^2 + ...
                        zz.^2*obj.PixSizeXYZ(3)^2/Rg^2));
                % apply filter - remove mean before the application of Gaussian
                % filter
                 Pic1gauss = imfilter(obj.Pic1 - mean(obj.Pic1(:)), GaussFilter);
             else
                 Pic1gauss = obj.Pic1;
             end
             
            % auto threshold at this point
            % Thr = graythresh(Pic1gauss(:))*max(Pic1gauss(:));
             Thr = mymultithresh(Pic1gauss(:), 10);
             CellNo = 0;
             for i = 1:length(Thr),
                 Pic1gaussBW = (Pic1gauss > Thr(end-i+1));
                 CC = bwconncomp(Pic1gaussBW);
                 if CC.NumObjects > MaxCellNo, % get somewhat more cells to discard the smallest
                     break;
                 end
             end
                
            
            % find cells
            for i = 1:CC.NumObjects,
                [II, JJ, KK] = ind2sub(size(Pic1gauss), CC.PixelIdxList{i});
                R = Pic1gauss(CC.PixelIdxList{i})'*[JJ II KK]/sum(Pic1gauss(CC.PixelIdxList{i}));
                XYZgyr = Pic1gauss(CC.PixelIdxList{i})'*[JJ II KK].^2/sum(Pic1gauss(CC.PixelIdxList{i}));
                XYZgyr = sqrt(XYZgyr - R.^2).*obj.PixSizeXYZ';
                if all(XYZgyr(1:2) > (CellMinDia/2)*sqrt(3/5)),
                    Cell1.maxPixValue = max(Pic1gauss(CC.PixelIdxList{i}));
                    Cell1.numPix = length(CC.PixelIdxList{i});
                    Cell1.sumPix = sum(Pic1gauss(CC.PixelIdxList{i}));
                    Cell1.XYZ_pix = R;
                    Cell1.XYZ_um = R.*obj.PixSizeXYZ';
                    Cell1.XYZgyr_um = XYZgyr; 
                    Cell1.PixelIdxList = CC.PixelIdxList{i};
                    if ~isprop(obj, 'Cells'),
                        o = addprop(obj, 'Cells');
                        obj.dynprop_meta.('Cells') = o;
                        obj.Cells = Cell1;
                    else 
                        obj.Cells(end+1) = Cell1;
                    end                
                end                
            end
            
            % select the brightest objects
            [~, JJ] = sort([obj.Cells.(SelectField)], 'descend');
            obj.Cells = obj.Cells(JJ(1:min(MaxCellNo, end)));
        end
        
        function DoShowCells(obj, varargin)
            % get contours
            % prepare Binary Stack
            BW = zeros(size(obj.Pic1));
            for i = 1:length(obj.Cells)
                BW(obj.Cells(i).PixelIdxList) = 1;
            end
            %B = bwboundaries(BW);
            B = BW - imerode(BW, true(3));
            B = B > 0.5;
            for i = 1:obj.NoZplanes,
                subplot(1, 2, 1)
                imagesc(obj.Pic1(:, :, i))
                axis image
                
                subplot(1, 2, 2)
                Pic1withContours = obj.Pic1(:, :, i);
                Pic1withContours(B(:, :, i)) = min(Pic1withContours(:));
                imagesc(Pic1withContours)
                axis image
                
                title(num2str(i));
                figure(gcf);
                pause
            end
            
        end
        
        function [Pic1out, Pic2out, Xscale, Yscale] = DoShowImage(obj, varargin)
              if isempty(obj.Pic1),
                  Pic1out = [];
                  Pic2out = [];
                  warning('No data in Pic1');
                  return
              end
              
              
              CenterPlane = round((size(obj.Pic1, 3)+1)/2);
              PlaneNo = ParseInputs('PlaneNo', CenterPlane, varargin); 
              OutliersFrac = ParseInputs('OutliersFrac', 0.001, varargin); 
              PlotPic2 = ParseInputs('Plot Pic2?', 1, varargin); 
              ImType = ParseInputs('ImageType', ImageType.Image, varargin); 
              PhotonData = ParseInputs('PhotonData', [], varargin);
              gate = ParseInputs('Gate_ns', [0 inf], varargin);
              binSize = ParseInputs('Bins', 1, varargin);
              binSize = round(binSize);
              if (length(binSize(:)) > 2),
                  error('Bins should be a scalar or a vector of two values1');
              elseif (length(binSize(:)) == 1)
                  binSize(2) = binSize(1);
              end
                  
              
              Pic2out = [];
              
              if isempty(obj.Pic2), 
                  PlotPic2 = 0;
              end
              
              % to determine limits kick out the outliers
              if (ImType == ImageType.Image),
                Pic1out = obj.Pic1(:, :, PlaneNo);
                if binSize(1) > 1,
                    tempPic = cumsum([zeros(size(Pic1out(1, :))); Pic1out]);
                    Pic1out = diff(tempPic(1:binSize(1):end, :));
                end
                
                if binSize(2) > 1,
                    tempPic = cumsum([zeros(size(Pic1out(:, 1))), Pic1out], 2);
                    Pic1out = diff(tempPic(:, 1:binSize(2):end), 1, 2);
                end

                
                if PlotPic2,
                    Pic2out = obj.Pic1(:, :, PlaneNo);
                    if binSize(1) > 1,
                        tempPic = cumsum([zeros(size(Pic2out(1, :))); Pic2out]);
                        Pic2out = diff(tempPic(1:binSize(1):end, :));
                    end

                    if binSize(2) > 1,
                        tempPic = cumsum([zeros(size(Pic2out(:, 1))), Pic2out], 2);
                        Pic2out = diff(tempPic(:, 1:binSize(2):end), 1, 2);
                    end
                end
              else % (if ImType ~= ImageType.Image) 
                  if isempty(PhotonData),
                      error('No photon data supplied!');
                  end
                  phInd1 = obj.PhotonInd(PlaneNo).phInd1;
                  delayTimes = [PhotonData(PlaneNo).delayTime];
                  
                  Pic1out = zeros(floor(size(obj.Pic1(:, :, PlaneNo))./binSize));
                  if PlotPic2
                      Pic2out = zeros(floor(size(obj.Pic2(:, :, PlaneNo))./binSize));
                      phInd2 = obj.PhotonInd(PlaneNo).phInd2;
                  end
                  
                  for ii = 1:binSize(1):(size(obj.Pic1, 1)-binSize(1)+1),
                      for jj = 1:binSize(2):(size(obj.Pic1, 2)-binSize(2)+1),
                          dT = [];
                          dT2 = [];
                          for k = 0:(binSize(1)-1),
                              dT = [dT delayTimes(phInd1(ii+k, jj, 1):phInd1(ii+k, jj+binSize(2)-1, 2))];
                              if PlotPic2,
                                  dT2 = [dT2 delayTimes(phInd2(ii+k, jj, 1):phInd2(ii+k, jj+binSize(2)-1, 2))];
                              end
                          end
                          dT = dT(dT >= gate(1) & dT <= gate(2));
                          if PlotPic2,
                              dT2 = dT2(dT2 >= gate(1) & dT2 <= gate(2));
                          end

                          if (ImType == ImageType.LifeTime)
                              Pic1out((ii-1)/binSize(1)+1, (jj-1)/binSize(2)+1) = mean(dT);
                              if PlotPic2,
                                  Pic2out((ii-1)/binSize(1)+1, (jj-1)/binSize(2)+1) = mean(dT2);
                              end
                          elseif (ImType == ImageType.LifeTimeMedian)
                              Pic1out((ii-1)/binSize(1)+1, (jj-1)/binSize(2)+1) = median(dT);
                              if PlotPic2,
                                  Pic2out((ii-1)/binSize(1)+1, (jj-1)/binSize(2)+1) = median(dT2);
                              end

                          elseif (ImType == ImageType.Gated)
                              Pic1out((ii-1)/binSize(1)+1, (jj-1)/binSize(2)+1) = length(dT);
                              if PlotPic2,
                                  Pic2out((ii-1)/binSize(1)+1, (jj-1)/binSize(2)+1) = length(dT2);
                              end
                          end
                      end
                  end                  
              end
              
              Pic11 = Pic1out;
              Pic11(isnan(Pic11)) = [];
              pic1sorted = sort(Pic11(:));
              Cmin = pic1sorted(round(length(pic1sorted)*OutliersFrac)+1);
              Cmax = pic1sorted(round(length(pic1sorted)*(1-OutliersFrac)));
              Xscale = (0:(size(Pic1out, 2)-1))*binSize(2)*obj.PixSizeHVP(1);
              Yscale = (0:(size(Pic1out, 1)-1))*binSize(1)*obj.PixSizeHVP(2);
              %figure;
              
              if PlotPic2,
                  subplot(1, 2, 1);
              else 
                  subplot(1, 1, 1);
              end;
              
              %imagesc(Xscale, Yscale, Pic1out, [Cmin Cmax])
              imagesc(Pic1out, [Cmin Cmax])
              set(gca, 'YDir', 'normal');
              axis image
              
              if PlotPic2,
                  subplot(1, 2, 2)
                  %imagesc(Xscale, Yscale, Pic2out, [Cmin Cmax]);
                  imagesc(Pic2out, [Cmin Cmax]);
                  set(gca, 'YDir', 'normal');
                  axis image
              end
              
              figure(gcf);
          end % DoShowImage
          
          function obj = DoImageCorrelation(obj, varargin)
              MaxDistancePixels = ParseInputs('MaxDistancePix', inf, varargin);
              SubtractShotNoise = ParseInputs('Subtract shot noise', true, varargin);
              %check that the resolution is the same: assuming XY scan
              if ~(isempty(obj.PixSizeXYZ)),
                  if strcmp(obj.ScanType, 'XYscan'),
                      pixCols = obj.PixSizeXYZ(1);
                      pixRows = obj.PixSizeXYZ(2);
                  elseif strcmp(obj.ScanType, 'YZscan'),
                      pixCols = obj.PixSizeXYZ(2);
                      pixRows = obj.PixSizeXYZ(3);
                  elseif strcmp(obj.ScanType, 'ZXscan'),
                      pixCols = obj.PixSizeXYZ(3);
                      pixRows = obj.PixSizeXYZ(1);
                  else
                      warning('unrecognized scan type!');
                      pixCols = 1;
                      pixRows = 1;
                  end
                  
                  if (abs(pixCols - pixRows)/pixCols > 0.01),
                      warndlg('Pixels size in two directions  are not equal!')
                  end
              else
                  warning('No pixel size data: assuming unit size pixels');
                  pixCols = 1;
                  pixRows = 1;                  
              end
              
              isPic1 = ~isempty(obj.Pic1);
              isPic2 = ~isempty(obj.Pic2);
              
              if ~(isPic1 | isPic2), 
                  warning('No images in Pic1 and Pic2!');
                  return;
              end
               
              if isempty(obj.ROI.x),
                  if isPic1,
                    BW = ones(size(obj.Pic1(:, :, 1)));
                  else
                    BW = ones(size(obj.Pic2(:, :, 1)));
                  end
              else
                  if isPic1,
                    BW = roipoly(obj.Pic1(:, :, 1), obj.ROI.x, obj.ROI.y); % for normalization
                  else
                    BW = roipoly(obj.Pic2(:, :, 1), obj.ROI.x, obj.ROI.y); % for normalization                   
                  end
              end;
              
              if isinf(MaxDistancePixels), 
                  BWcut  = BW;
              else
                disk = strel('disk', MaxDistancePixels, 0);
                BWcut = imerode(BW, disk);
              end
              
              BWcorr = xcorr2(double(BWcut), double(BW));
              if isPic1, 
                BW3cut = repmat(BWcut, [1, 1, size(obj.Pic1, 3)]);
                Pic1cut = obj.Pic1.*BW3cut;
                if isPic2,
                    Pic2cut = obj.Pic2.*BW3cut;
                end
              else
                BW3cut = repmat(BWcut, [1, 1, size(obj.Pic2, 3)]);
                Pic2cut = obj.Pic2.*BW3cut;
              end
              
              
              %BWcorr = xcorr2(double(BW));
              for i = 1:size(BW3cut, 3),
                  if isPic1,
                     meanPic1 = xcorr2(double(BWcut), obj.Pic1(:, :, i))./BWcorr;
                     meanCut1 = xcorr2(Pic1cut(:, :, i), double(BW))./BWcorr;
                     Pic1corr = xcorr2(Pic1cut(:, :, i), obj.Pic1(:, :, i))./(BWcorr.*meanPic1.*meanCut1) - 1;
                     middle1 = round((size(Pic1corr)+1)/2);
                  end
                  
                  if isPic2,
                      meanPic2 = xcorr2(double(BWcut), obj.Pic2(:, :, i))./BWcorr;
                      meanCut2 = xcorr2(Pic2cut(:, :, i), double(BW))./BWcorr;
                      Pic2corr = xcorr2(Pic2cut(:, :, i), obj.Pic2(:, :, i))./(BWcorr.*meanPic2.*meanCut2) - 1;
                      middle2 = round((size(Pic2corr)+1)/2);
                  end
                  
                  if SubtractShotNoise,
                     %Pic1corr(middle1) = Pic1corr(middle1) - 
                  end
              end
              
              %get disk of different radii
              zeroCorr = zeros(size(Pic1corr));
              if isinf(MaxDistancePixels),
                  MD = min(size(BW))/2-1; %sqrt(sum(BW(:)));
              else 
                  MD = MaxDistancePixels;
              end
              r = 0:(2*MD);
              for j = r;
                  disk = strel('disk', j, 0).getnhood;
                  nrm(j+1) = sum(disk(:));
                  msk = zeroCorr;
                  
                  if isPic1,
                      msk(middle1(1)+((-j):j), middle1(2)+((-j):j)) = disk;
                      msk3 = repmat(msk, [1, 1, size(obj.Pic1, 3)]);
                      Pic1cr = Pic1corr;
                      Pic1cr(~msk3) = 0;
                      cumG1(j+1) = sum(Pic1cr(:));
                  end
                  
                  if isPic2,    
                      msk(middle2(1)+((-j):j), middle2(2)+((-j):j)) = disk;
                      msk3 = repmat(msk, [1, 1, size(obj.Pic2, 3)]);
                      Pic2cr = Pic2corr;
                      Pic2cr(~msk3) = 0;                  
                      cumG2(j+1) = sum(Pic2cr(:));
                  end
              end
              %organize structure
              ImageCorrelation.BW = BW;
              ImageCorrelation.BWcut = BWcut;
              ImageCorrelation.BWcorr = BWcorr;
              ImageCorrelation.Npix = diff([0 nrm]);
              if isPic1,
                ImageCorrelation.G1 = diff([0 cumG1])./ImageCorrelation.Npix;
              end
              
              if isPic2,
                ImageCorrelation.G2 = diff([0 cumG2])./ImageCorrelation.Npix;
              end
              
              ImageCorrelation.r = [0 (r(1:(end-1))+r(2:end))/2]*pixCols;
              if (~isprop(obj, 'ImageCorr')),
                obj.addprop('ImageCorr');
              end
              obj.ImageCorr = ImageCorrelation;
              obj.DoPlotImageCorrelation;
          end %DoImageCorrelation
          
          function obj = DoPlotImageCorrelation(obj, varargin)
              PlotType = ParseInputs('Plot Type', 'plot', varargin);
              subplot(1, 1, 1)
              switch lower(PlotType)
                  case {'plot','loglog', 'semilogx', 'semilogy'}
                      isG1G2 = isfield(obj.ImageCorr, {'G1', 'G2'});
                      
                      if isG1G2(1), % is G1
                        eval([lower(PlotType) '(obj.ImageCorr.r, obj.ImageCorr.G1, ''-o'')']);
                        hold all;
                      end
                      
                      if isG1G2(2)
                        eval([lower(PlotType) '(obj.ImageCorr.r, obj.ImageCorr.G2, ''-o'')']);
                      end
                      hold off;
                                               
                        %eval([lower(PlotType) '(obj.ImageCorr.r, [obj.ImageCorr.G1; obj.ImageCorr.G2], ''-o'')']);
                      
                      xlabel('r (\mu m)', 'FontSize', 18);
                      ylabel('G', 'FontSize', 18);
                      set(gca, 'FontSize', 16);
                      figure(gcf)                      
                  otherwise
              end
          end %DoPlotImageCorrelation
          
          function obj = DoLineCorrelation(obj, varargin)
              xcorrOption = ParseInputs('xcorrOption', 'biased', varargin);
               if ~(isempty(obj.PixSizeXYZ)),
                  if strcmp(obj.ScanType, 'XYscan'),
                      pixCols = obj.PixSizeXYZ(1);
                      pixRows = obj.PixSizeXYZ(2);
                  elseif strcmp(obj.ScanType, 'YZscan'),
                      pixCols = obj.PixSizeXYZ(2);
                      pixRows = obj.PixSizeXYZ(3);
                  elseif strcmp(obj.ScanType, 'ZXscan'),
                      pixCols = obj.PixSizeXYZ(3);
                      pixRows = obj.PixSizeXYZ(1);
                  else
                      warning('unrecognized scan type!');
                      pixCols = 1;
                      pixRows = 1;
                  end
               else
                  warning('unrecognized scan type!');
                  pixCols = 1;
                  pixRows = 1;
               end


              for i = 1:size(obj.Pic1, 3),
                  for j = 1:size(obj.Pic1, 1),
                      LineCorrelation.G1(j, :, i) = xcorr(obj.Pic1(j, :, i), xcorrOption);
                      LineCorrelation.G2(j, :, i) = xcorr(obj.Pic2(j, :, i), xcorrOption);
                      LineCorrelation.meanCount1(j, i) = mean(obj.Pic1(j, :, i));
                      LineCorrelation.meanCount2(j, i) = mean(obj.Pic2(j, :, i));
                      mid = (size(LineCorrelation.G1, 2) + 1)/2;
                      %subtract shot noise
                      LineCorrelation.G1(j, mid, i) = LineCorrelation.G1(j, mid, i) - LineCorrelation.meanCount1(j, i);
                      LineCorrelation.G2(j, mid, i) = LineCorrelation.G2(j, mid, i) - LineCorrelation.meanCount2(j, i);
                  end
              end
              %mid = (size(LineCorrelation.G1, 2) + 1)/2;
              LineCorrelation.G1 = LineCorrelation.G1(:, mid:end, :);
              LineCorrelation.G2 = LineCorrelation.G2(:, mid:end, :);
              LineCorrelation.aveG1 = mean(mean(LineCorrelation.G1),3);
              LineCorrelation.aveG2 = mean(mean(LineCorrelation.G2),3);
              LineCorrelation.normG1 = LineCorrelation.aveG1/LineCorrelation.aveG1(1);
              LineCorrelation.normG2 = LineCorrelation.aveG2/LineCorrelation.aveG2(1);
              LineCorrelation.r = pixRows*(0:(size(obj.Pic1, 2)-1));
%              LineCorrelation.x = pixRows*(-(size(obj.Pic1, 2)-1):(size(obj.Pic1, 2)-1));              
              if (~isprop(obj, 'LineCorrelation')),
                obj.addprop('LineCorrelation');
              end
              obj.LineCorrelation = LineCorrelation;
              subplot(1, 1, 1);
              plot(LineCorrelation.r, LineCorrelation.aveG1, ...
                  LineCorrelation.r, LineCorrelation.aveG2);
              xlabel('r (\mu m)', 'FontSize', 18);
              figure(gcf)
          end  
          
          function obj = DoImportDataFromPython(obj, folderPaths, varargin)
              OverwriteProps = ParseInputs('OverwriteProps', {}, varargin); % cell array of property-value pairs to assign to obj
              if ischar(folderPaths),
                  folderPaths = {folderPaths}; % convert to cell array for similar treatment
              end
              
              for i = 1:length(folderPaths),
                  fPi = folderPaths{i};
                  if exist(fPi, 'file') == 7, % first version all fields saved in separate files in the same directory
                      if folderPaths{i}(end) ~= filesep,
                          folderPaths{i} = [folderPaths{i} filesep];
                      end
                      fnames = dir([folderPaths{i} '*.mat']);
                      fnames = {fnames.name};
                      for j = 1:length(fnames),
                          fpathToFile = [folderPaths{i} fnames{j}];
                          tempS = load(fpathToFile);
                          if ~isprop(obj, tempS.StrucName),
                              obj.addprop(tempS.StrucName);
                              obj.(tempS.StrucName) = tempS
                          end
                          if i == 1,
                              if isfield(tempS, 'index'),
                                  obj.(tempS.StrucName)(tempS.index) = tempS;
                              else
                                  obj.(tempS.StrucName) = tempS;
                              end
                          else
                              obj.(tempS.StrucName)(end+1) = tempS;
                          end
                      end
                  elseif exist(fPi, 'file') == 2, %new versions: all data in the same file
                      tempS = load(fPi);
                      if i == 1,
                          if ~isprop(obj, tempS.StrucName),
                              obj.addprop(tempS.StrucName);
                          end
                          obj.(tempS.StrucName) = tempS;                            
                      else
                          obj.(tempS.StrucName)(i) = tempS;
                      end
                  end
                  
                  if i == 1, % second file might not have ROI
                      if isprop(obj, 'ROI'),
                          obj.X0 = obj.ROI(1).x0;
                          obj.Y0 = obj.ROI(1).y0;
                          obj.Z0 = obj.ROI(1).z0;
                      elseif isfield(tempS, 'ROI'),
                          if isfinite(tempS.ROI.x0),
                              obj.X0 = tempS.ROI(1).x0;
                              obj.Y0 = tempS.ROI(1).y0;
                              obj.Z0 = tempS.ROI(1).z0;
                          else
                              obj.X0 = 1;
                              obj.Y0 = 1;
                              obj.Z0 = 1;
                          end

                      else
                          obj.X0 = 1;
                          obj.Y0 = 1;
                          obj.Z0 = 1;
                      end
                  end

              end
              
              if ~isempty(OverwriteProps),
                  for i = 1:2:length(OverwriteProps),
                      newProp = OverwriteProps{i};
                      newVal = OverwriteProps{i+1};
                      if ~isprop(obj, newProp),
                          obj.addprop(newProp);
                      end
                      obj.(newProp) = newVal;                                                  
                  end
              end
              
              
              if isprop(obj, 'CellCoordImport'),
                  tempS = obj.('CellCoordImport');
                  for i = 1:length(tempS)
                      tempS(i).Xhsc = tempS(i).X_um/obj.PixSizeXYZ(1) + 1 - obj.X0;
                      tempS(i).Yhsc = tempS(i).Y_um/obj.PixSizeXYZ(2) + 1 - obj.Y0;
                      tempS(i).Zhsc = tempS(i).Z_um/obj.PixSizeXYZ(3) + 1 - obj.Z0;
                  end
                  obj.('CellCoordImport') = tempS; 
              elseif isprop(obj, 'PythonDataImport')
                  tempS = obj.('PythonDataImport');
                  for i = 1:length(tempS)
                      tempS(i).Xhsc = tempS(i).X_um/obj.PixSizeXYZ(1) + 1 - obj.X0;
                      tempS(i).Yhsc = tempS(i).Y_um/obj.PixSizeXYZ(2) + 1 - obj.Y0;
                      tempS(i).Zhsc = tempS(i).Z_um/obj.PixSizeXYZ(3) + 1 - obj.Z0;
                      for j = 1:length(tempS(i).DirectStretchList),
                          tempS(i).DirectStretchList{j}.Xhsc = ...
                              tempS(i).DirectStretchList{j}.X/obj.PixSizeXYZ(1) + 1 - obj.X0;
                          tempS(i).DirectStretchList{j}.Yhsc = ...
                              tempS(i).DirectStretchList{j}.Y/obj.PixSizeXYZ(2) + 1 - obj.Y0;
                      end
                      
                      for j = 1:length(tempS(i).ConfinementList),
                          tempS(i).ConfinementList{j}.Xhsc = ...
                              tempS(i).ConfinementList{j}.X/obj.PixSizeXYZ(1) + 1 - obj.X0;
                          tempS(i).ConfinementList{j}.Yhsc = ...
                              tempS(i).ConfinementList{j}.Y/obj.PixSizeXYZ(2) + 1 - obj.Y0;
                      end
                  end
                  obj.('PythonDataImport') = tempS; 
              end 
              
              if isprop(obj, 'DirectedStretch'),
                  tempS = obj.('DirectedStretch');
                  for i = 1:length(tempS)
                      tempS(i).Xhsc = tempS(i).X/obj.PixSizeXYZ(1) + 1 - obj.X0;
                      tempS(i).Yhsc = tempS(i).Y/obj.PixSizeXYZ(2) + 1 - obj.Y0;
                  end
                  obj.('DirectedStretch') = tempS;             
              end 
          end
          
          function obj = DoGeneralMethod(obj, varargin)
              [fpathToFile, isImportFromPython] = ParseInputs('Import From Python', '', varargin);
              if isImportFromPython,
                  tempS = load(fpathToFile);
                  if ~isprop(obj, tempS.StrucName),
                      obj.addprop(tempS.StrucName);
                  end
                  obj.(tempS.StrucName) = tempS;                  
              end
              
              [refPointInPix, isRefPointInPix] = ParseInputs('Add Reference Point in pxls', [], varargin);
              if isRefPointInPix,
                  obj.X0 = refPointInPix(1);
                  obj.Y0 = refPointInPix(2);
                  obj.Z0 = refPointInPix(3);
                  if isprop(obj, 'CellCoordImport'),
                      tempS = obj.('CellCoordImport');
                      tempS.Xhsc = tempS.X_um/obj.PixSizeXYZ(1) + 1 - obj.X0;
                      tempS.Yhsc = tempS.Y_um/obj.PixSizeXYZ(2) + 1 - obj.Y0;
                      tempS.Zhsc = tempS.Z_um/obj.PixSizeXYZ(3) + 1 - obj.Z0;
                      obj.('CellCoordImport') = tempS;
                  end
              end
                                
          end
          
         function [shiftI, shiftJ, shiftK, mm] = GetAlignShiftTwoStacks(obj, img1, img2, varargin)
            AlignmentType = ParseInputs('AlignmentType', 'ImageDifference', varargin);
             maxStep = ParseInputs('maxStep', [5 5 2], varargin); %ijk
            [minStep, isminStep] = ParseInputs('minStep', [-5 -5 -2], varargin); %ijk
            if ~isminStep,
                minStep = - maxStep;
            end
            
            if strcmp(AlignmentType, 'ImageDifference'),
                img1 = double(img1((-minStep(1)+1):(end - maxStep(1)), ...
                    (-minStep(2)+1):(end - maxStep(2)), (-minStep(3)+1):(end - maxStep(3))));
                img1 = img1/mean(img1(:));
                score = [];

                for i = minStep(1):maxStep(1),
                    for j = minStep(2):maxStep(2),
                        for k = minStep(3):maxStep(3),
                            img2shifted = double(img2((-minStep(1)+1 + i):(end - maxStep(1) + i), ...
                                (-minStep(2)+1 + j):(end - maxStep(2) + j), ...
                                (-minStep(3)+1 + k):(end - maxStep(3) + k)));
                            img2shifted = img2shifted/mean(img2shifted(:));
                            score(i - minStep(1) + 1, j - minStep(2) + 1, ...
                                k - minStep(3) + 1) = sum((img1(:) - img2shifted(:)).^2);
                        end
                    end
                end
                score =  score/prod(size(img1))
            elseif strcmp(AlignmentType, 'CorrCoeff'),
                score = [];
                for i = minStep(1):maxStep(1),
                    for j = minStep(2):maxStep(2),
                        for k = minStep(3):maxStep(3),
                            %[i j k]
                            img1shifted = double(img1((max((-i+1), 1):min(end, end-i)), ...
                                (max((-j+1), 1):min(end, end-j)), (max((-k+1), 1):min(end, end-k))));
                            img1shifted = img1shifted/mean(img1shifted(:));

                            img2shifted = double(img2(max(1, i+1):min((end + i), end), ...
                                max(1, j+1):min((end + j), end), max(1, k+1):min((end + k), end)));
                            img2shifted = img2shifted/mean(img2shifted(:));
                            score(i - minStep(1) + 1, j - minStep(2) + 1, ...
                                k - minStep(3) + 1) = -img1shifted(:)'*img2shifted(:)/prod(size(img1shifted));
                        end
                    end
                end
            end

%
            [mm, mi] = min(score(:)); 
            [imin, jmin, kmin] = ind2sub(size(score), mi);
            shiftI = imin + minStep(1) - 1;
            shiftJ = jmin + minStep(2) - 1;
            shiftK = kmin + minStep(3) - 1;
         end
         
         function shiftIJK = GetAlignShiftAllStacks(obj, varargin)
             Trange = ParseInputs('Trange', 1:obj.NoTimePoints, varargin);
             TimePointToAlignTo = ParseInputs('TimePointToAlignTo', Trange(1), varargin);
             ColorToAlignBy = ParseInputs('ColorToAlignBy', 1, varargin);
             Ch1 = obj.DoGetSubStack('Color', ColorToAlignBy, 'Time', TimePointToAlignTo);
             img1 = Ch1.Pic1;
             shiftIJK = [];
             for T = Trange,
                 disp(T)
                 Ch2 = obj.DoGetSubStack('Color', ColorToAlignBy, 'Time', T);
                 img2 = Ch2.Pic1;
                 [shiftI, shiftJ, shiftK, score] = obj.GetAlignShiftTwoStacks(img1, img2, varargin{:});
                 disp([shiftI, shiftJ, shiftK, score])
                 shiftIJK = [shiftIJK; [shiftI, shiftJ, shiftK, score]];
                 %img1 = img2;
             end
             
         end
         
         function objAligned = DoAlignByShiftMatrix(obj, shiftIJK, varargin)
            KeepPlanesAboveBelow = ParseInputs('KeepPlanesAboveBelow', 2, varargin);
            objAligned = ImageStack_Class;
            objAligned.Pic1 = cast(objAligned.Pic1, 'like', obj.Pic1);
            %copy properties from wholeStack
            props = properties(obj);
            for ii = 1:length(props),
                if ~strcmp(props{ii}, 'Pic1'),
                    if ~isprop(objAligned, props{ii}),
                        objAligned.addprop(props{ii});
                    end
                    objAligned.(props{ii}) = obj.(props{ii});
                end
            end
            Zhsc = [obj.PythonDataImport.Zhsc];
            Zhsc = Zhsc(:) - shiftIJK(1:length(Zhsc), 3);
            % plane range
            maxShift = max(shiftIJK(:, 1:3));
            minShift = min(shiftIJK(:, 1:3));
            minShift(3) = 1 - round((min(Zhsc)-KeepPlanesAboveBelow))
            maxShift(3) = obj.NoZplanes - round(max(Zhsc)+KeepPlanesAboveBelow)
            objAligned.NoZplanes = obj.NoZplanes + minShift(3) - maxShift(3);
            objAligned.NoPlanesPerTimePoint = objAligned.NoZplanes*objAligned.NoColors;
            objAligned.TotalPlanes = objAligned.NoPlanesPerTimePoint*objAligned.NoTimePoints;
            newSize = size(obj.Pic1)  + (minShift - maxShift);
            newSize(3) = objAligned.TotalPlanes;
            objAligned.Pic1 = zeros(newSize,  'like', obj.Pic1); %cast(NaN(size(obj.Pic1) - maxShift + minShift), 'like', obj.Pic1);
            planeStat = cell(objAligned.NoZplanes, 1);
            
            for planeNo = 1:obj.TotalPlanes,
                [C, Z, T] = obj.NoInStackToCZT(planeNo);
                shift = shiftIJK(T, :);
                Znew = Z + minShift(3)-shift(3);
                if (Znew > 0) & (Znew <= objAligned.NoZplanes),
                    planeNoAligned = objAligned.CZTtoNoInStack('Color', C, 'Z', Znew, 'Time', T);
                    objAligned.Pic1(:, :, planeNoAligned) = ...
                        obj.Pic1((1 - minShift(1) + shift(1)):(end - maxShift(1) + shift(1)),...
                            (1-minShift(2)+shift(2)):(end - maxShift(2) + shift(2)), planeNo);
                    planeStat{Znew} = [planeStat{Znew} T];
                end
            end
            
            commonTimePoints = 1:objAligned.NoTimePoints;
            for i  = 1:length(planeStat),
                planeStat{i} = unique(planeStat{i});
                commonTimePoints = intersect(commonTimePoints, planeStat{i});
            end
            objAligned.addprop('planeData');
            objAligned.planeData.timePoints = planeStat;
            objAligned.planeData.commonTimePoints = commonTimePoints;
                        
            if isprop(objAligned, 'PythonDataImport')
              tempS = objAligned.('PythonDataImport');
              Jhsc0 = 0;
              for i = 1:length(tempS)
                  Jhsc = (1:length(tempS(i).Xhsc)) + Jhsc0;
                  tempS(i).Xhsc = tempS(i).Xhsc(:) + minShift(2) - shiftIJK(Jhsc,2);
                  tempS(i).Yhsc = tempS(i).Yhsc(:) + minShift(1) - shiftIJK(Jhsc, 1);
                  tempS(i).Zhsc = tempS(i).Zhsc(:) + minShift(3) - shiftIJK(Jhsc, 3);
                  for j = 1:length(tempS(i).DirectStretchList),
                      jj = tempS(i).DirectStretchList{j}.('Stretch Indices')+1; % conversion of 0 based Python array to Matlab 1 based
                      jj = Jhsc0 + (jj(1):jj(2));
                      tempS(i).DirectStretchList{j}.Xhsc = ...
                          tempS(i).DirectStretchList{j}.Xhsc(:) + minShift(2) - shiftIJK(jj, 2);
                      tempS(i).DirectStretchList{j}.Yhsc = ...
                          tempS(i).DirectStretchList{j}.Yhsc(:) + minShift(1) - shiftIJK(jj, 1);
                  end

                  for j = 1:length(tempS(i).ConfinementList),
                      jj = tempS(i).ConfinementList{j}.('Stretch Indices');
                      jj = Jhsc0 + (jj(1):jj(2));

                      tempS(i).ConfinementList{j}.Xhsc = ...
                          tempS(i).ConfinementList{j}.Xhsc(:) + minShift(2) - shiftIJK(jj, 2);
                      tempS(i).ConfinementList{j}.Yhsc = ...
                          tempS(i).ConfinementList{j}.Yhsc(:) + minShift(1) - shiftIJK(jj, 1);
                  end
                  Jhsc0 = Jhsc(end);
              end
              objAligned.('PythonDataImport') = tempS; 
            end
            
            objAligned.addprop('shiftIJK');
            objAligned.shiftIJK = shiftIJK;

         end
              
    end
    
    methods (Static)
        function SeriesData = DoGetSeriesData(fpath)
            status = bfCheckJavaPath(1);
            assert(status, ['Missing Bio-Formats library. Either add bioformats_package.jar '...
                'to the static Java path or add it to the Matlab path.']);

            % Initialize logging
            bfInitLogging();

            % Get the channel filler
            r = bfGetReader(fpath, 1);
            numSeries = r.getSeriesCount();
            result = cell(numSeries, 2);

            globalMetadata = r.getGlobalMetadata();
            
            for s = 1:numSeries,
                r.setSeries(s - 1);
                SeriesData(s).numImages = r.getImageCount();
                SeriesData(s).seriesMetadata = r.getSeriesMetadata();
                SeriesData(s).sizeX = r.getSizeX();
                SeriesData(s).sizeY = r.getSizeY();
                disp(SeriesData(s).seriesMetadata)
            end
            r.close()
            
        end
        
   
    end
    methods (Access = protected)
        function obj = DoBinning1dim(obj, dim, binSize),
            BinSize = [1 1 1];
            BinSize(dim) = binSize;
            if (binSize > 1),
                pic1 = zeros(floor(size(obj.Pic1)./BinSize));
                JJ = (0:(size(pic1, dim)-1))*binSize;
                for i = 1:binSize,         
                    if (dim == 1), 
                        addPic = obj.Pic1(i + JJ, :, :);
                    elseif (dim == 2),
                        addPic = obj.Pic1(:, i + JJ, :);
                    else 
                        addPic = obj.Pic1(:, :, i + JJ);
                    end
                    pic1 = pic1 + addPic;
                end
            obj.Pic1 = pic1;
            obj.PixSizeXYZ(dim) = obj.PixSizeXYZ(dim)*binSize;
            obj.NoZplanes = size(obj.Pic1, 3);

            end
        end
        
        function [PixIdx, i, j, k, Vol, bound3Dprojection, A, iShell, jShell, kShell] = ...
                Apply3Dthreshold(obj, img3D, thr3D)
            cc = bwconncomp(imfill(img3D > thr3D, 'holes'));
            numPixels = cellfun(@numel,cc.PixelIdxList);
            [biggest,idx] = max(numPixels);
            [i, j, k] = ind2sub( size(img3D), cc.PixelIdxList{idx});
            %get the max contour of Pix3D
            img3D = zeros(size(img3D));
            img3D(cc.PixelIdxList{idx}) = 1;
            % find shell
            H = strel('disk', 1, 0);
            [iShell, jShell, kShell] = ind2sub(size(img3D), find((img3D - (imerode(img3D, H))>0.5)));
            
            img2D = (sum(img3D, 3) > 0.5);
            bound3Dprojection = bwboundaries(img2D);
            bound3Dprojection = bound3Dprojection{1};
            PixIdx = cc.PixelIdxList{idx};
            Vol = length(PixIdx);
            A = sum(img2D(:));
        end
        
%         function [CellParam, boundSoma, warningMes] = Apply2Dthreshold(obj, CellParamIn, img, thresh, se, rectSize, Xbw, Ybw)
        function [CellParam, boundSoma, warningMes] = Apply2Dthreshold(obj, CellParamIn, img, thresh, se, ROIx1bw, ROIy1bw)

            CellParam = CellParamIn;
            quantImg = imquantize(img, thresh);
            NoThresholds = length(thresh);
            [boundSoma, boundInd] = bwboundaries(imopen((quantImg == NoThresholds+1), se), 'noholes');
            if NoThresholds == 2, 
                boundLam = bwboundaries(imopen(((quantImg == NoThresholds)| (quantImg == NoThresholds+1)), se), 'noholes');
            end
            if isempty(boundSoma),
                CellParam.Zc = NaN;
                CellParam.widthZ_um = NaN;
                CellParam.brightnessZ = NaN;
                CellParam.backgroundZ = NaN;
                CellParam.contour.X = NaN;
                CellParam.contour.Y = NaN;
                CellParam.area = NaN;
                CellParam.Xc = NaN; %centroid
                CellParam.Yc = NaN; %centroid
                CellParam.perimeter = NaN;
                CellParam.D_um = NaN;
                CellParam.shapeFactor = NaN; 
                CellParam.Xcm = NaN;
                CellParam.Ycm = NaN;
                CellParam.polygeom.iner = [];
                CellParam.polygeom.cpmo = [];
                CellParam.Pix3D = [];
            
                return
            end
            
            [max_size, max_index] = max(cellfun('size', boundSoma, 1)); % find largest
            distFromCenter = cell2mat(cellfun(@mean, boundSoma, 'UniformOutput', false)) - repmat(size(img), length(boundSoma), 1)/2; % find closest to the center
            [~, closest_index] = min(sum(distFromCenter.^2, 2));
            brightest_ind = 0;
            brightest_mean =  mean(img(boundInd == 0));
            for bi = 1:length(boundSoma), % find brightest
                 newmean = mean(img(boundInd == bi));
                 if newmean > brightest_mean,
                     brightest_ind  = bi;
                     brightest_mean = newmean;
                 end
            end
            
            if ~((max_index == brightest_ind)&(max_index == closest_index)),
                if brightest_ind == closest_index,
                    max_index = closest_index;
                    warningMes = 'Warning: not the largest selection';
                elseif (max_index == brightest_ind),
                    warningMes = 'Warning: not the closest to center';
                elseif (max_index == closest_index),
                    warningMes = 'Warning: not the brightest selection';
                else 
                    max_index = closest_index;
                    warningMes = 'Warning: not the brightest or largest selection';
                end
            else
                warningMes = '';
            end
            
            boundSoma = boundSoma{max_index};
            if NoThresholds == 2,
                % find boundLam containing boundSoma
                for bL = 1:length(boundLam),
                    bLgood(bL) = any(inpolygon(boundSoma(:, 1), boundSoma(:, 2), boundLam{bL}(:, 1), boundLam{bL}(:, 2)));
                end
                bLgood = find(bLgood);
                if length(bLgood)~=1,
                    warning(['problem identifying soma at T = ' num2str(T)]);
                end
                boundSoma = boundLam{bLgood};
            end
  %             BW = poly2mask(boundSoma(:, 2), boundSoma(:, 1), size(img, 1), size(img, 2));
%             CellArea = sum(BW(:));
%             Xsoma = boundSoma(:,2)+round(Xbw - rectSize/2)-1;
%             Ysoma = boundSoma(:,1)+round(Ybw - rectSize/2)-1;
            Xsoma = boundSoma(:,2)+ROIx1bw-1;
            Ysoma = boundSoma(:,1)+ROIy1bw-1;

            BW = poly2mask(boundSoma(:, 2), boundSoma(:, 1), size(img, 1), size(img, 2));
            [geom, iner, cpmo] = polygeom(Xsoma, Ysoma);
            CellParam.contour.X = Xsoma;
            CellParam.contour.Y = Ysoma;
            CellParam.area = geom(1);
%             CellParam.Xc = sum(BW, 1)*(1:size(BW, 2))'/sum(BW(:))+round(Xbw - rectSize/2)-1;%geom(2); %centroid
%             CellParam.Yc = (1:size(BW, 1))*sum(BW, 2)/sum(BW(:))+round(Ybw - rectSize/2)-1; %geom(3); %centroid
            CellParam.Xc = sum(BW, 1)*(1:size(BW, 2))'/sum(BW(:))+ROIx1bw-1;%geom(2); %centroid
            CellParam.Yc = (1:size(BW, 1))*sum(BW, 2)/sum(BW(:))+ROIy1bw-1; %geom(3); %centroid

            CellParam.perimeter = geom(4);
            CellParam.D_um = geom(4)/pi*obj.PixSizeXYZ(1);
            CellParam.shapeFactor = geom(4)^2/(4*pi*geom(1)); 
            %imgMasked = img0.*BW;
            imgMasked = img.*BW;
%             CellParam.Xcm = sum(imgMasked, 1)*(1:size(imgMasked, 2))'/sum(imgMasked(:))+round(Xbw - rectSize/2)-1;
%             CellParam.Ycm = (1:size(imgMasked, 1))*sum(imgMasked, 2)/sum(imgMasked(:))+round(Ybw - rectSize/2)-1;
            CellParam.Xcm = sum(imgMasked, 1)*(1:size(imgMasked, 2))'/sum(imgMasked(:))+ROIx1bw-1;
            CellParam.Ycm = (1:size(imgMasked, 1))*sum(imgMasked, 2)/sum(imgMasked(:))+ROIy1bw-1;

            CellParam.polygeom.iner = iner;
            CellParam.polygeom.cpmo = cpmo;
        end

    end
    
end

