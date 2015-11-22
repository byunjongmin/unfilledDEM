function main
% @file main.m
% @brief Function to extract stream longitudinal profiles from unfilled
% DEMs
%
% @version 0.3.2 / 2015-11-22
% @author Jongmin Byun
%==========================================================================

%% Load DEM

clear all

% Constants
INPUT_DIR = '../data/input';
OUTPUT_DIR = '../data/output';
DEMFileName = 'n37_e128_1arc_v3.tif';
DEMFilePath = fullfile(INPUT_DIR,DEMFileName);
[DEM,R] = geotiffread(DEMFilePath);
% note that R means MapCellsReference or GeographicPostingsReference
% according to CoordinateSysteamType

% DEM basic properties
mRows = R.RasterSize(1,1);
nCols = R.RasterSize(1,2);
if strcmp(R.CoordinateSystemType,'planar')
    
    dX = R.CellExtentInWorldX;
    dY = R.CellExtentInWorldY;
    
else % strcmp(R,'geographic')
    
    CIR_EQUATOR = 40075; % circumference at equator [km]
    CIR_MERIDIAN = 40008; % meridional circumference [km]
    rLat = (R.LatitudeLimits(2) - R.LatitudeLimits(1)) * 0.5;
    
    dY = CIR_MERIDIAN * 1000 / 360 * R.SampleSpacingInLatitude;
    dX = CIR_EQUATOR * 1000 / 360 * cos(rLat * pi / 180) * R.SampleSpacingInLatitude;

end

DEM = double(DEM);

%% Define target drainage

% For debug
clear all
INPUT_DIR = '../data/input';
dataFileName = 'a_load_DEM_2015-11-22.mat';
dataFilePath = fullfile(INPUT_DIR,dataFileName);
load(dataFilePath);

% for the domain surrounded by null
nanMask = (DEM == 32767);
DEM(nanMask) = inf;
DEMArea = ~nanMask;

% note that you should check the location of the outlet of the test domain.
% It should be located on the boundary of a drainage. If it is within the
% drainage, you should modify the DEM

% extract the boundary of the target drainage
s = strel('square',3); % structural element when eroding image
eDEMArea = imerode(DEMArea,s); 
DEMAreaBnd = DEMArea & ~eDEMArea;
DEMAreaBndIdx = find(DEMAreaBnd);

% identify the coordinate of the outlet
DEMAreaBndElev = DEM(DEMAreaBndIdx);
[~,minElevIdx] = min(DEMAreaBndElev);
outletIdx = DEMAreaBndIdx(minElevIdx); % outlet on the boundary of drainage
[outletY,outletX] = ind2sub([mRows,nCols],outletIdx);

% Note that the elevation of the main outlet should be the lowest.
DEM(outletY,outletX) = min(DEM(:)) - 0.1;

% define target drainage
targetDrainage = (~nanMask);
targetDrainage(outletY,outletX) = false;

% Main body

% Smooth the imported DEM, until the flat is removed
orgDEM = DEM; % original DEM

WANT_SMOOTH = false;
fSize = 10;
dCon = 0.95; % decay constant
nSmooth = 150;
if WANT_SMOOTH == true
    
    % h = 1/9*ones(3); % mean filter
    % bell shaped weight
    h = ones(fSize*2+1,fSize*2+1);
    for i = 1:fSize
        h(i:(end+1)-i,i:(end+1)-i) = dCon^(fSize-i);
    end
        
    h = 0.25 * ones(7,7);
    h(2:end-1,2:end-1) = 0.5;
    h(3:end-2,3:end-2) = 1;
    h = 1/sum(sum(h)) * h;
    
    for j=1:nSmooth
        
        DEM = filter2(h,DEM);
        
    end
    
    figure(1); clf;
    imagesc(DEM);
    set(gca,'DataAspectRatio',[1 1 1]);
    colorbar;
    title('Smoothed Digital Elevation Model');

end
% Assign flow directions to the DEM using D8 algorithm
% Note that, to use CalcSDSFlow function, elevation of outer region of the
% target drainage should be inf.
[steepestDescentSlope,slopeAllNbr,SDSFlowDirection,SDSNbrY,SDSNbrX] ...
    = CalcSDSFlow(DEM,dX,dY);

% Process flat and sinks for flow to continue to move
flatRegMap = ProcessFlat(DEM,targetDrainage,slopeAllNbr);

rT = 0;
afterNFlat = inf;
while afterNFlat > 0

    fprintf('The remaining flat cell is %4.2f\n',afterNFlat);
    % smooth only flat. To prevent an infinite loop, change the size of
    % moving window when flat region does not reduce 
    if rT > 5
        nNbr = 5;
    else
        nNbr = 3;
    end    
    h = ones(nNbr,nNbr)/nNbr^2;
    smtDEM = filter2(h,DEM);
    DEM(flatRegMap == true) = smtDEM(flatRegMap == true);
    [steepestDescentSlope,slopeAllNbr,SDSFlowDirection,SDSNbrY,SDSNbrX] ...
        = CalcSDSFlow(DEM,dX,dY);
    flatRegMap = ProcessFlat(DEM,targetDrainage,slopeAllNbr);
    
    oldNFlat = afterNFlat;
    afterNFlat = numel(find(flatRegMap == true));
    if oldNFlat >= afterNFlat
        rT = rT + 1;
    else
        rT = 0;
    end
    
end

% for debug
figure(2); clf;

subplot(1,2,1);
imagesc(flatRegMap);
set(gca,'DataAspectRatio',[1 1 1]);
title('Distribution of Flat Region');

subplot(1,2,2);
diffDEM = orgDEM - DEM; % for debug
imagesc(diffDEM);
set(gca,'DataAspectRatio',[1 1 1]);
colorbar;
title('Difference in Elevation after Smoothing');

%% Remove isolated area

DEM_BW = false(mRows,nCols); % binary of DEM
DEM_BW(~isinf(DEM)) = true;
CC = bwconncomp(DEM_BW); % identify connected components of valid cells
if CC.NumObjects > 1
    DEM(CC.PixelIdxList{2:5}) = inf; % remove the connected components except for the target area
    targetDrainage(CC.PixelIdxList{2:end}) = false;
end

% for debug

% % For debug
% clear all
% INPUT_DIR = '../data/input';
% dataFileName = 'b_IS_IT_PART_2015-11-18.mat';
% dataFilePath = fullfile(INPUT_DIR,dataFileName);
% load(dataFilePath);

IS_IT_PART = false;
IS_BND_INF = true;

if IS_IT_PART == true
    
    % For debug for a part of DEM
    tYMin = 558; tYMax = 658;
    tXMin = 117; tXMax = 168;

    DEM = DEM(tYMin:tYMax,tXMin:tXMax);
    [mRows,nCols] = size(DEM);
    
    if IS_BND_INF == true
        
        % for the domain filled with only elevations
        nanMask = true(mRows,nCols);
        nanMask(2:mRows-1,2:nCols-1) = false;
        DEMArea = ~nanMask;
        DEM(nanMask) = inf;

        % extract the boundary of the target drainage
        s = strel('square',3); % structural element when eroding image
        eDEMArea = imerode(DEMArea,s); 
        DEMAreaBnd = DEMArea & ~eDEMArea;
        DEMAreaBndIdx = find(DEMAreaBnd);

        % identify the coordinate of the outlet
        DEMAreaBndElev = DEM(DEMAreaBndIdx);
        [~,minElevIdx] = min(DEMAreaBndElev);
        outletIdx = DEMAreaBndIdx(minElevIdx); % outlet on the boundary of drainage
        [outletY,outletX] = ind2sub([mRows,nCols],outletIdx);

        % Note that the elevation of the main outlet should be the lowest.
        DEM(outletY,outletX) = min(DEM(:)) - 0.1;
        
        % define target drainage
        targetDrainage = (~nanMask);
        targetDrainage(outletY,outletX) = false;
        
    else
        
        % for the domain filled with only elevations
        nanMask = true(mRows,nCols);
        nanMask(2:mRows-1,2:nCols-1) = false;
        DEMArea = ~nanMask;
        DEM(nanMask) = 0;
        
        % define target drainage
        targetDrainage = (~nanMask);
        
    end
    
    [steepestDescentSlope,slopeAllNbr,SDSFlowDirection,SDSNbrY,SDSNbrX] ...
        = CalcSDSFlow(DEM,dX,dY);
    % flatRegMap = ProcessFlat(DEM,targetDrainage,slopeAllNbr);
    
    % for debug
    figure(3); clf;
    set(gcf, 'Color',[1,1,1]);

    imagesc(DEM);
    set(gca,'DataAspectRatio',[1 1 1]);
    colorbar;
    title('Digital Elevation Model');
    
end

% Identify depressions and their outlets, then modify each depression
% outlet's flow direction to go downstream when flows are overspilled
% over the depression

% Note : To use ProcessSink function, make a mask excepting the main outlet
% and outer region of the target drainage
[flood ... % flooded region map
,m1SDSNbrY,m1SDSNbrX ... % modified for flooed region and its outlet
,m2SDSNbrY,m2SDSNbrX ... % modified only for flooded region's outlets 
,fldRegID ... % flooded region ID
,nFldRegCells ... % each flooded region cells number
,subFldRegOutlet ... % sub-flooded region outlet location
,subFldRegOutInfo ... % information on sub-flooded region outlets
,subFldRegID ... % sub-flooded region ID
,regionalMin ... % regional minimum map
,sharedOutlet] ... % shared outlet map
    = ProcessSink(DEM,targetDrainage,slopeAllNbr ...
        ,SDSNbrY,SDSNbrX,SDSFlowDirection);

    
% frequency distribution of the types of sub-flooded region's outlet
figure(4); clf;
set(gcf, 'Color',[1,1,1]);
h = histogram(subFldRegOutlet(subFldRegOutlet > 0));
xlabel('Type of Sub-flooded Region Outlet');
ylabel('Frequency');
grid on

% set boundary
fXMin = 1; fXMax = nCols;
fYMin = 1; fYMax = mRows;

figure(4); clf;
set(gcf, 'Color',[1,1,1]);

imagesc(fldRegID(fYMin:fYMax,fXMin:fXMax));
colorbar;
set(gca,'DataAspectRatio',[1 1 1]);
title('Flooded Region ID');

figure(5); clf;
set(gcf, 'Color',[1,1,1]);

imagesc(subFldRegID(fYMin:fYMax,fXMin:fXMax));
colorbar;
set(gca,'DataAspectRatio',[1 1 1]);
title('Sub-flooded Region ID');

figure(6); clf;
set(gcf, 'Color',[1,1,1]);

imagesc(subFldRegOutlet(fYMin:fYMax,fXMin:fXMax));
colorbar;
set(gca,'DataAspectRatio',[1 1 1]);
title('Sub-flooded Region Outlet');

% Modify flow direction of the cells within each depression

% Assign flow direction to the cells in each depression region
[m2SDSNbrY,m2SDSNbrX ... % flow direction modified cells along the path to each regional minima
,mFlowDir_SubFldReg ... % flow direction modified cell within a sub-flooded region 
,mFlowDir_Saddle ... % flow direction modified cell on a saddle
,subFldRegTree ...
,fldRegInfo] ... % informaiton on flooded regions: ID, outlet index, number of sub-flooded region, outlet's elevation
    = AssignFlowDirInFldReg(m2SDSNbrY,m2SDSNbrX,subFldRegOutInfo,DEM ...
    ,slopeAllNbr,regionalMin,fldRegID,subFldRegID,sharedOutlet);

% display the tree DB of flooded region
disp(subFldRegTree.tostring);

%% Calculate upstream cells number of all cells within every depression

% for debug
clear all
INPUT_DIR = '../data/input';
dataFileName = 'a_CalcUpstreamCells_2015-11-20_2.mat';
dataFilePath = fullfile(INPUT_DIR,dataFileName);
load(dataFilePath);

% A. Calculate the number of upstream cells out of flooded region
nUpstreamCellsWithFldReg = CalcUpstreamCellsWithFldReg(DEM,targetDrainage ...
    ,flood,m1SDSNbrY,m1SDSNbrX,fldRegID,nFldRegCells);

% B. Calculate flow accumulation within each flooded region
nUpstreamCells ...
    = CalcUpstreamCellsInSubFldReg(subFldRegTree,mFlowDir_SubFldReg ...
    ,fldRegID,nUpstreamCellsWithFldReg,fldRegInfo ...
    ,m2SDSNbrY,m2SDSNbrX,subFldRegID,DEM,regionalMin);   

% Visualization for the number of upstream cells
figure(7);
set(gcf,'Color',[1 1 1])

subplot(1,3,1)
imagesc(nUpstreamCellsWithFldReg(2:end-1,2:end-1));
title('Upstream Cells No. with Depressions');
axis image
set(gca,'YTick',[],'XTick' ,[])
colormap(flipud(colormap(gray)))
colorbar

subplot(1,3,2)
imagesc(nUpstreamCells(2:end-1,2:end-1));
title('Upstream Cells No.');
axis image
set(gca,'YTick',[],'XTick' ,[])
colormap(flipud(colormap(gray)))
colorbar

subplot(1,3,3)
imagesc(log(nUpstreamCells(2:end-1,2:end-1)));
title('Upstream Cells No. [Log]');
axis image
set(gca,'YTick',[],'XTick' ,[])
colormap(flipud(colormap(gray)))
colorbar

figure(8)
set(gcf,'Color',[1 1 1])
imagesc(nUpstreamCells(2:end-1,2:end-1));
title('Upstream Cells No.');
axis image
set(gca,'YTick',[],'XTick' ,[])
colormap(flipud(colormap(gray)))
colorbar

%% extract and analyze stream longitudinal profiles

%% draw a stream longitudinal profile on the interesting stream path

% % for debug
% clear all
% INPUT_DIR = '../data/input';
% dataFileName = 'a_CalcUpstreamCells_2015-11-18.mat';
% dataFilePath = fullfile(INPUT_DIR,dataFileName);
% load(dataFilePath);

figure(2); clf;
set(gcf, 'Color',[1,1,1]);

imagesc(DEM);
set(gca,'DataAspectRatio',[1 1 1]);
colorbar;
title('Digital Elevation Model');

% identify stream path using given initial and end point coordinates
initY = 617; initX = 719;
endY = 852; endX = 299;
[streamPath,distFromInit] = RecordStrPath(initY,initX,endY,endX ...
    ,mRows,nCols,m2SDSNbrY,m2SDSNbrX,dY,dX);

% stream longitudinal profile using raw DEM
elev = DEM(streamPath(:));
nCells = numel(streamPath);

% draw figure
figure(9); clf;
set(gcf, 'Color',[1,1,1]);

subplot(3,1,1)
stairs(distFromInit,elev);

% grid on
title('Longitudinal Stream Profile')
xlabel('Distance From Divide [m]')
ylabel('Elevation [m]')
xlim([0 distFromInit(end)])
ylim([min(elev),max(elev)])

% smoothing stream profiles

considerNbrForElev = [50,150]; % number of considering neighbor cells
% note that you can compare the differences due to the number of
% considering neighbor cells
nCNbr = numel(considerNbrForElev); % number of cases

smoothedElev = SmoothingElev(considerNbrForElev,distFromInit,DEM,streamPath);

subplot(3,1,2)

cc = jet(nCNbr);
for ithLine = 1:nCNbr
    plot(distFromInit,smoothedElev(:,ithLine),'color',cc(ithLine,:));
    hold on
end

grid on
title('Longitudinal Stream Profile')
xlabel('Distance From Divide [m]')
ylabel('Elevation [m]')
xlim([0 distFromInit(end)])
ylim([min(elev),max(elev)])

% draw stream gradient

considerNbrForSlope = 45;
chosenSizeForSlope = 1;
nCNbrForSlope = numel(considerNbrForSlope);

slopePerDist = CalcSlope(considerNbrForSlope,distFromInit,smoothedElev ...
    ,chosenSizeForSlope);

subplot(3,1,3)

cc = jet(nCNbrForSlope);
for ithLine = 1:nCNbrForSlope
    plot(distFromInit,slopePerDist(:,ithLine),'color',cc(ithLine,:));
    hold on
end

grid on
title('Stream Gradient Profile')
xlabel('Distance From Initiaion [m]')
ylabel('Slope')
xlim([0 distFromInit(end)])
ylim([min(slopePerDist(:,1)),max(slopePerDist(:,1))])

% draw corrected upstream area profiles on the interesting stream paths

% choose a stream gradient profile
chosenSlopeForUpArea_Slope = 1;

% calculate upstream area
upstreamAreaProf = nUpstreamCells(streamPath(:)) .* dX .* dY;

% correction of upstream area prof
prevUpArea = upstreamAreaProf(1);
initUpArea = inf;
isInitFirst = true;

for ithCell = 2:nCells
    
    if upstreamAreaProf(ithCell) < prevUpArea
        
        if isInitFirst == true
        
            initUpArea = prevUpArea;
            isInitFirst = false;
            
        end
        
        upstreamAreaProf(ithCell) = initUpArea;        
        
    end
    
    if upstreamAreaProf(ithCell) > initUpArea
        
        isInitFirst = true;
        
    end       
    
    prevUpArea = upstreamAreaProf(ithCell);

end

figure(12); clf;
set(gcf, 'Color',[1,1,1]);

subplot(2,1,1)
plot(distFromInit,upstreamAreaProf);

set(gca,'YScale','log');
grid on
title('Upstream Area Profile')
xlabel('Distance From Divide [m]')
ylabel('Upstream Area [m^2]')
xlim([0 distFromInit(end)])
ylim([min(upstreamAreaProf),max(upstreamAreaProf)])

% draw upslope area - slope relationship

subplot(2,1,2)
scatter(upstreamAreaProf,slopePerDist(:,chosenSlopeForUpArea_Slope) ...
    ,'Marker','*');

set(gca,'XScale','log','YScale','log');
grid off
title('Upstream Area - Slope')
xlabel('Upstream Area [m^2]')
ylabel('Slope')
xlim([upstreamAreaProf(1),upstreamAreaProf(end)])

% export stream gradient profile to ArcGIS

fileName = strcat(num2str(initY),'_',num2str(initX),'_',num2str(endY),'_',num2str(endX),'_profile_gradient');

ithProfileDistMap = nan(mRows,nCols);
for ithCell = 1:nCells
    ithProfileDistMap(streamPath(ithCell)) = slopePerDist(ithCell);    
end
% IPDFDFilePath = fullfile(OUTPUT_DIR,fileName);
IPDFDFilePath = fullfile('/Users/cyberzen/Dropbox/temp',fileName);
ExportRasterAsArcGrid(IPDFDFilePath,ithProfileDistMap,R.DeltaX ...
    ,R.XLimWorld(1,1),R.YLimWorld(1,1));

