function main
% @file main.m
% @brief Function to extract stream longitudinal profiles from unfilled
% DEMs
%
% @version 0.5.1 / 2015-11-26
% @author Jongmin Byun
%==========================================================================

%% Load DEM

clear all

% Constants
INPUT_DIR = '../data/input';
OUTPUT_DIR = '../data/output';
DEMFileName = 'tareaR50m.tif';
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

% mark a mask of nan
nanMask = (DEM == 32767 | DEM == -9999); % avoid null values of DEM
nanMask(DEM == 0) = true; % include 0 as null value
orgDEM = DEM; % original DEM

% draw DEM
figure(1); clf;

DEM(nanMask) = nan;
imagesc(DEM);
DEM(nanMask) = orgDEM(nanMask);

set(gca,'DataAspectRatio',[1 1 1]);
colorbar;
title('Digital Elevation Model');

%% Remove flat cells in DEM by adding random elevation values and smoothing
% using a specific filter

% % for the MATLAB without MappingToolbox
% dataFileName = 'a_load_DEM_2015-11-26.mat';
% dataFilePath = fullfile(INPUT_DIR,dataFileName);
% load(dataFilePath);

% add small random elevation values to DEM

randE = rand(mRows,nCols);
DEM = DEM + randE * 0.00001;

% smooth DEM
% bell shaped weight window
fSize = 2; % radius of window
transRatio = 1; % ratio of transfering values to next neighbors
h = ones(fSize*2-1,fSize*2-1);
for i = 1:fSize
    h(i:end-(i-1),i:end-(i-1)) = transRatio^(fSize-i);
end
h = 1/sum(h(:)) * h;
% smooth using the filter
nSmooth = 1;
for i=1:nSmooth
    DEM = filter2(h,DEM);
end

% check whether flat cells are removed

% assign flow directions to the DEM using D8 algorithm
% note that, for the calculation of flow direction for the target domain,
% fill nan mask area with lower elevation values
minElev = min(DEM(~nanMask));
DEM(nanMask) = minElev - 1;
[steepestDescentSlope,slopeAllNbr,SDSFlowDirection,SDSNbrY,SDSNbrX] ...
    = CalcSDSFlow(DEM,dX,dY);
DEM(nanMask) = minElev + 1;

% make a map of flat cells
flatRegMap = ProcessFlat(DEM,~nanMask,slopeAllNbr);

nFlatCells = sum(flatRegMap(:));
if nFlatCells > 0
    error('Error. There remain flat cells in DEM\n');
end

% display difference between the original and smoothed dEM
diffDEM = orgDEM - DEM;

figure(2); clf;

subplot(1,2,1);
imagesc(flatRegMap(fSize*2:end-fSize*2,fSize*2:end-fSize*2));
set(gca,'DataAspectRatio',[1 1 1]);
colorbar;
title('Flat cells in DEM');

subplot(1,2,2);
diffDEM(nanMask) = nan;
imagesc(diffDEM(fSize*2:end-fSize*2,fSize*2:end-fSize*2));
set(gca,'DataAspectRatio',[1 1 1]);
colorbar;
title('Difference between the original and smoothed DEM');

%% for a test domain
IS_IT_PART = false;
if IS_IT_PART == true

    % coordinates of the test domain
    tYMin = 558; tYMax = 658;
    tXMin = 117; tXMax = 168;
    
    % cut DEM
    orgDEM = DEM(tYMin:tYMax,tXMin:tXMax); % original DEM
    DEM = DEM(tYMin:tYMax,tXMin:tXMax); % smoothed DEM
    
    [mRows,nCols] = size(DEM);
    
    % mark a mask of nan
    nanMask = (DEM == 32767 | DEM == -9999); % avoid null values of DEM
    nanMask(DEM == 0) = true; % include 0 as null value
    
    % for debug
    figure(3); clf;
    set(gcf, 'Color',[1,1,1]);

    DEM(nanMask) = nan;
    imagesc(DEM);
    DEM(nanMask) = orgDEM(nanMask);
    
    set(gca,'DataAspectRatio',[1 1 1]);
    colorbar;
    title('DEM of the test domain');
    
end

%% Define target drainage and treat boundaries of DEM
% note that you can choose whether the boundaries are higher or lower

% first, include the area out of nanMask
targetDrainage = (~nanMask); % target drainage

IS_BND_LOWER = false;

if IS_BND_LOWER == true
    
    DEM(~targetDrainage) = min(DEM(targetDrainage)) - 0.1;
    
else % IS_BND_LOWER == false

    % pick the outlet of the target drainage

    DEM(~targetDrainage) = inf;
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

%% Assign flow direction to the target area in DEM using D8 algorithm

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

%% Modify flow direction of the cells within each depression

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

% % for debug
% clear all
% INPUT_DIR = '../data/input';
% dataFileName = 'a_CalcUpstreamCells_2015-11-20_2.mat';
% dataFilePath = fullfile(INPUT_DIR,dataFileName);
% load(dataFilePath);

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
% set(gca,'YTick',[],'XTick' ,[])
colormap(flipud(colormap(gray)))
colorbar

subplot(1,3,2)
imagesc(nUpstreamCells(2:end-1,2:end-1));
title('Upstream Cells No.');
axis image
% set(gca,'YTick',[],'XTick' ,[])
colormap(flipud(colormap(gray)))
colorbar

subplot(1,3,3)
imagesc(log(nUpstreamCells(2:end-1,2:end-1)));
title('Upstream Cells No. [Log]');
axis image
% set(gca,'YTick',[],'XTick' ,[])
colormap(flipud(colormap(gray)))
colorbar

figure(8)
set(gcf,'Color',[1 1 1])
imagesc(nUpstreamCells(2:end-1,2:end-1));
title('Upstream Cells No.');
axis image
% set(gca,'YTick',[],'XTick' ,[])
colormap(flipud(colormap(gray)))
colorbar

%% Extract and analyze stream longitudinal profiles

%% Draw a stream longitudinal profile on the interesting stream path

% for debug
clear all
INPUT_DIR = '../data/input';
dataFileName = 'SRTM_1_dCon_0.5_20151126.mat';
dataFilePath = fullfile(INPUT_DIR,dataFileName);
load(dataFilePath);

% slopePerDistTotal = []; % for total sloperPerDist statistics

% draw DEM
figure(9); clf;
set(gcf, 'Color',[1,1,1]);

imagesc(DEM);
set(gca,'DataAspectRatio',[1 1 1]);
colorbar;
title('Digital Elevation Model');

% input the initiation and end point for a profile
initY = 1753; initX = 3220;
endY = 3181; endX = 1951;

% set the number of considered neighbor
considerNbrForSlope = 3:10:500;
chosenProfile = 1; % for smoothing

% chose range for identify knickpoint
initX_Knick = 16;
endX_Knick = 31;

% choose a stream gradient profile for corrected upstream area profiles
chosenSlopeForUpArea_Slope = 25;

% binned area-slope relationship
% make elevation values with fixed vertical interval.
% note: if you want to see subtle change in relatively flat region, lower
% contour interval
contInterval = 1; % contour interval
% determine logarithmic binned average slopes
logBinSize = 0.01;

[upstreamAreaProf,slopePerDist] ...
    = AnalyzeStreamProfile(initY,initX,endY,endX ...
    ,chosenProfile,considerNbrForSlope,initX_Knick,endX_Knick ...
    ,chosenSlopeForUpArea_Slope,contInterval,logBinSize...
    ,fSize,orgDEM,nUpstreamCells,mRows,nCols,m2SDSNbrY,m2SDSNbrX,dY,dX);

ithStreamGradientMap = nan(mRows,nCols);

for i = 1:nCells
    
    ithStreamGradientMap(streamPath(i)) = slopePerDist(i,chosenSlopeForUpArea_Slope);    

end

%% regression analysis

% within a specific reach
minArea = 1*10^6;
maxArea = 2*10^9;
rangeIdx = find(upstreamAreaProf > minArea & upstreamAreaProf < maxArea);

polyfit(log(upstreamAreaProf(rangeIdx)),log(slopePerDist(rangeIdx,chosenSlopeForUpArea_Slope)),1)

% for chosen locations
cX = [3.483*10^7,8.733*10^7,1.764*10^8,5.717*10^8,1.130*10^9];
cY = [0.04147,0.01987,0.01445,0.005069,0.003747];

polyfit(log(cX),log(cY),1)

%% export stream gradient profile to ArcGIS

fileName = strcat(num2str(initY),'_',num2str(initX),'_',num2str(endY),'_',num2str(endX),'_profile_gradient');
% IPDFDFilePath = fullfile(OUTPUT_DIR,fileName);
IPDFDFilePath = fullfile('/Volumes/DATA_BAK/WORKSPACE/Project/NewKnickPoints/Raster',fileName);

% ithStreamGradientMap(isnan(ithStreamGradientMap)) = 32767;
[~,R] = geotiffread(DEMFilePath);
geotiffwrite(IPDFDFilePath,ithStreamGradientMap,R);

