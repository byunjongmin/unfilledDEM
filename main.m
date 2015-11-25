function main
% @file main.m
% @brief Function to extract stream longitudinal profiles from unfilled
% DEMs
%
% @version 0.4.0 / 2015-11-25
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

figure(1); clf;
imagesc(DEM);
set(gca,'DataAspectRatio',[1 1 1]);
colorbar;
title('Smoothed Digital Elevation Model');

%% Check and remove flat cells in DEM

% Assign flow directions to the DEM using D8 algorithm
% Note that, to use CalcSDSFlow function, elevation of outer region of the
% target drainage should be inf.
nanMask = (DEM == 32767 | DEM == -9999);
DEM(nanMask) = inf;
[steepestDescentSlope,slopeAllNbr,SDSFlowDirection,SDSNbrY,SDSNbrX] ...
    = CalcSDSFlow(DEM,dX,dY);

% identify flat cells in DEM
flatRegMap = ProcessFlat(DEM,~nanMask,slopeAllNbr);

orgDEM = DEM; % original DEM
afterNFlat = inf;
while afterNFlat > 0

    fprintf('The number of remaining flat cells is %4.0f\n',afterNFlat);
        
    % add small random elevation values to DEM
    randE = rand(mRows,nCols);
    DEM = DEM + randE;

    % update elevation values only for the flat cells
    DEM(flatRegMap == true) = DEM(flatRegMap == true);

    [steepestDescentSlope,slopeAllNbr,SDSFlowDirection,SDSNbrY,SDSNbrX] ...
        = CalcSDSFlow(DEM,dX,dY);

    flatRegMap = ProcessFlat(DEM,~nanMask,slopeAllNbr);

    % update the number of flat cells in DEM
    afterNFlat = numel(find(flatRegMap == true));
    
end

% for debug
diffDEM = orgDEM - DEM;

figure(2); clf;
imagesc(diffDEM);
set(gca,'DataAspectRatio',[1 1 1]);
colorbar;
title('Flat cells in DEM');

%% Define target drainage

% for a test domain
IS_IT_PART = false;
if IS_IT_PART == true

    % coordinates of the test domain
    tYMin = 558; tYMax = 658;
    tXMin = 117; tXMax = 168;

    % for debug
    figure(3); clf;
    set(gcf, 'Color',[1,1,1]);

    imagesc(DEM);
    set(gca,'DataAspectRatio',[1 1 1]);
    colorbar;
    title('DEM of the test domain');

end

% define a target drainage according to the type of DEM
IS_ISLAND = false;
IS_BND_INF = true;
LOCATE_OUTLET = true;

bndMask = true(mRows,nCols);
bndMask(2:mRows-1,2:nCols-1) = false;

if IS_ISLAND == true

    % for the domain surrouned by null values
    targetDrainage = (~nanMask); % target drainage
    DEM(~targetDrainage) = inf;

else % IS_ISLAND == false

    % for the domain filled with only elevations
    targetDrainage = (~nanMask); % target drainage
    DEM(bndMask) = inf;

end

if IS_BND_INF == false

    % for the lower boundary
    DEM(bndMask) = min(min(DEM(targetDrainage))) - 1;

end

if LOCATE_OUTLET == true

    % locate the outlet of the target drainage

    % extract the boundary of the target drainage
    s = strel('square',3); % structural element when eroding image
    dilatedTarget = imerode(targetDrainage,s); 
    targetBnd = targetDrainage & ~dilatedTarget;
    targetBndIdx = find(targetBnd);

    % identify the coordinate of the outlet
    targetBndElev = DEM(targetBndIdx);
    [~,minElevIdx] = min(targetBndElev);
    outletIdx = targetBndIdx(minElevIdx); % outlet on the boundary of drainage
    [outletY,outletX] = ind2sub([mRows,nCols],outletIdx);

    % Note that the elevation of the main outlet should be the lowest.
    DEM(outletY,outletX) = min(DEM(:)) - 0.1;

    % locate outlet of target drainage
    targetDrainage(outletY,outletX) = false;
    
end

[steepestDescentSlope,slopeAllNbr,SDSFlowDirection,SDSNbrY,SDSNbrX] ...
    = CalcSDSFlow(DEM,dX,dY);

%% Process sinks
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

