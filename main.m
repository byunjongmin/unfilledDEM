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

% draw DEM
figure(9); clf;
set(gcf, 'Color',[1,1,1]);

imagesc(DEM);
set(gca,'DataAspectRatio',[1 1 1]);
colorbar;
title('Digital Elevation Model');

% input the initiation and end point for a profile

initY = 717; initX = 2121;
endY = 2995; endX = 1691;

% identify stream path using given initial and end point coordinates
[streamPath,distFromInit] = RecordStrPath(initY,initX,endY,endX ...
    ,mRows,nCols,m2SDSNbrY,m2SDSNbrX,dY,dX);

% stream longitudinal profile using raw DEM
% note that we use the original DEM
profElev = orgDEM(streamPath(:));
nCells = numel(streamPath);

% draw figure
figure(10); clf;
set(gcf, 'Color',[1,1,1]);

subplot(3,1,1)
stairs(distFromInit,profElev);

% grid on
title('Longitudinal Stream Profile')
xlabel('Distance From Divide [m]')
ylabel('Elevation [m]')
xlim([0 distFromInit(end)])
ylim([min(profElev),max(profElev)])

%% smoothing stream profiles
% note that, if you input multiple values, you can compare the differences
% due to the number of considering neighbor cells

considerNbrForProf = fSize*2-1; % number of considering neighbor cells
nCNbr = numel(considerNbrForProf); % number of cases

smoothedProfElev = SmoothingElev(considerNbrForProf,distFromInit,DEM,streamPath);

subplot(3,1,2)

cc = jet(nCNbr);
for ithLine = 1:nCNbr
    plot(distFromInit,smoothedProfElev(:,ithLine),'color',cc(ithLine,:));
    hold on
end

grid on
title('Smoothed Longitudinal Stream Profile')
xlabel('Distance From Divide [m]')
ylabel('Elevation [m]')
xlim([0 distFromInit(end)])
ylim([min(profElev),max(profElev)])

%% Draw stream gradient

% chose base elevation of profile
inputProfElev = smoothedProfElev;
chosenProfile = 1;

% set the number of considered neighbor
considerNbrForSlope = 3:10:1000;
nCNbrForSlope = numel(considerNbrForSlope);

slopePerDist = CalcSlope(considerNbrForSlope,distFromInit,inputProfElev ...
    ,chosenProfile);

subplot(3,1,3)

cc = jet(nCNbrForSlope);
for ithLine = 1:nCNbrForSlope
    plot(distFromInit,slopePerDist(:,ithLine),'color',cc(ithLine,:));
    hold on
end

grid on
title('Stream Gradient')
xlabel('Distance From Initiaion [m]')
ylabel('Slope')
xlim([0 distFromInit(end)])
ylim([min(slopePerDist(:,1)),max(slopePerDist(:,1))])

% simple statistics of slopePerDist
minSlopePerDist = zeros(1,nCNbrForSlope);
maxSlopePerDist = zeros(1,nCNbrForSlope);
meanSlopePerDist = zeros(1,nCNbrForSlope);
stdSlopePerDist = zeros(1,nCNbrForSlope);

for i=1:nCNbrForSlope
    
    minSlopePerDist(1,i) = min(slopePerDist(:,i));
    maxSlopePerDist(1,i) = max(slopePerDist(:,i));
    meanSlopePerDist(1,i) = mean(slopePerDist(:,i));
    stdSlopePerDist(1,i) = std(slopePerDist(:,i));
    
end
    
figure(11); clf;

subplot(2,2,1)
scatter(considerNbrForSlope .* dX,minSlopePerDist);

grid on
title('Change in Minimum Stream Gradient')
xlabel('Distance (Dd)')
ylabel('Stream Gradient')

subplot(2,2,2)
scatter(considerNbrForSlope .* dX,maxSlopePerDist);

grid on
title('Change in Maximum Stream Gradient')
xlabel('Distance (Dd)')
ylabel('Stream Gradient')

subplot(2,2,3)
scatter(considerNbrForSlope .* dX,meanSlopePerDist);

grid on
title('Change in Mean Stream Gradient')
xlabel('Distance (Dd)')
ylabel('Stream Gradient')

subplot(2,2,4)
scatter(considerNbrForSlope .* dX,stdSlopePerDist);

grid on
title('Change in Standard Deviation Stream Gradient')
xlabel('Distance (Dd)')
ylabel('Standard Deviation')


%% draw corrected upstream area profiles on the interesting stream paths

% choose a stream gradient profile
chosenSlopeForUpArea_Slope = 25;

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

figure(11); clf;
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

%% export stream gradient profile to ArcGIS

ithStreamGradientMap = nan(mRows,nCols);

for i = 1:nCells
    
    ithStreamGradientMap(streamPath(i)) = slopePerDist(i,chosenSlopeForUpArea_Slope);    

end

fileName = strcat(num2str(initY),'_',num2str(initX),'_',num2str(endY),'_',num2str(endX),'_profile_gradient');
% IPDFDFilePath = fullfile(OUTPUT_DIR,fileName);
IPDFDFilePath = fullfile('/Volumes/DATA_BAK/WORKSPACE/Project/NewKnickPoints/Raster',fileName);

% ithStreamGradientMap(isnan(ithStreamGradientMap)) = 32767;
[~,R] = geotiffread(DEMFilePath);
geotiffwrite(IPDFDFilePath,ithStreamGradientMap,R);

%% Analyze slope-area relationship according to fixed vertical interval
% note: errors in DEM produce many depressions and negative slope values
% along a stream profile. Therefore, this algorithm concentrates on
% making a stream profile without depression. Using this stream profile,
% it calculates new slope values and log binnined average slope values.
% However using a smoothed stream proile as an input data makes analysis
% better. Because the smoothed stream profile has less depressions than
% raw stream profile data.

% make elevation values with fixed vertical interval.
% note: if you want to see subtle change in relatively flat region, lower
% contour interval
contInterval = 1; % contour interval
% fixed vertical interval
ipElevMin = 10 * floor(0.1 * profElev(end));
ipElevMax = 10 * ceil(0.1 * profElev(1));
ipElev = ipElevMax:-contInterval:ipElevMin;

% interpolated distance values according to fixed vertical interval
tStrProfElev = profElev;
% remove depressions in an input stream profile.
nCells = numel(distFromInit);
for a = nCells:-1:2
    if tStrProfElev(a-1) <= tStrProfElev(a)
        tStrProfElev(a-1) = tStrProfElev(a) + 0.000001;
    end
end
ipDistFromInit = interp1(tStrProfElev,distFromInit,ipElev);

% draw new stream profile with both elevation with fixed vertical interval
% and interpolated distance.

figure(12); clf;

% subplot(5,1,1)
% % for debug
% % plot(distFromInit,tStrProfElev,'o',ipDistFromInit,ipElev);
% plot(ipDistFromInit,ipElev);
% 
% % grid on
% title('Longitudinal Stream Profile')
% xlabel('Distance From Divide [m]')
% ylabel('Elevation [m]')
% xlim([0 max(ipDistFromInit)])
% ylim([min(ipElev),max(ipElev)])

% calculate new slope values using interpolated elevation and distance data
% note: basically new slope means local slope rather than 'trend slope'.
nCells2 = numel(ipElev);
newSlope = zeros(nCells2,1);
for a = 1
    newSlope(a) = (ipElev(a) - ipElev(a+1)) ...
        / (ipDistFromInit(a+1) - ipDistFromInit(a));
end
for a = 2:nCells2-1
    newSlope(a) = (ipElev(a-1) - ipElev(a+1)) ...
        / (ipDistFromInit(a+1) - ipDistFromInit(a-1));
end
for a = nCells2
    newSlope(a) = (ipElev(a-1) - ipElev(a)) ...
        / (ipDistFromInit(a) - ipDistFromInit(a-1));
end

% % draw newly calculated slope vs distance relationship
% subplot(5,1,2)
% plot(ipDistFromInit,newSlope,'*')
% 
% % grid on
% xlabel('Distance From Initiaion [m]')
% ylabel('Slope')
% xlim([0 max(ipDistFromInit)])
% ylim([min(newSlope),max(newSlope)])

% draw newly calculated upstream area vs elevation relationship
ipUpstreamArea = interp1(tStrProfElev,upstreamAreaProf,ipElev);

subplot(3,1,1)
plot(ipDistFromInit,ipUpstreamArea);

set(gca,'YScale','log');

% grid on
title('Upstream Area Profile')
xlabel('Distance From Divide [m]')
ylabel('Upstream Area [m^2]')
xlim([0 max(ipDistFromInit)])
ylim([min(ipUpstreamArea),max(ipUpstreamArea)])

% draw newly calculated upstream area vs slope relationship
subplot(3,1,2)
scatter(ipUpstreamArea,newSlope,'m+')

set(gca,'XScale','log','YScale','log');
% grid on
title('Upstream Area - Slope')
xlabel('Upstream Area [m^2]')
ylabel('Slope')
xlim([min(ipUpstreamArea),max(ipUpstreamArea)])

% determine logarithmic binned average slopes

logBinSize = 0.01;

minPower = floor(log10(min(ipUpstreamArea))/logBinSize) * logBinSize;
maxPower = ceil(log10(max(ipUpstreamArea))/logBinSize) * logBinSize;
diffPower = maxPower - minPower;

cellsNo3 = floor(diffPower/logBinSize);

tempArea = zeros(1,cellsNo3);
tempSlope = zeros(1,cellsNo3);

for a = 1:cellsNo3
    
    logAreaMin = minPower + (a-1) * logBinSize;
    logAreaMax = logAreaMin + logBinSize;
    
    tempArea(a) = 10^(logAreaMin + logBinSize * 0.5);
    
    j = find(ipUpstreamArea >= 10^logAreaMin ...
        & ipUpstreamArea < 10^logAreaMax);
    k = size(nonzeros(j));
    if k(1) > 0
        tempSlope(a) = 10^(mean(log10(newSlope(j))));
    end
end

% extract non-zero elements
lbarea = tempArea(find(tempSlope));
lbslope = tempSlope(find(tempSlope));

subplot(3,1,3)
loglog(lbarea,lbslope,'rs','MarkerFaceColor','r','MarkerSize',3)
% grid on
title('Upstream Area - Slope')
xlabel('Upstream Area [m^2]')
ylabel('Slope')
xlim([min(ipUpstreamArea),max(ipUpstreamArea)])

