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

% subplot(1,3,1)
% imagesc(nUpstreamCellsWithFldReg);
% title('Upstream Cells No. with Depressions');
% axis image
% % set(gca,'YTick',[],'XTick' ,[])
% colormap(flipud(colormap(gray)))
% colorbar
% 
% subplot(1,3,2)
% imagesc(nUpstreamCells);
% title('Upstream Cells No.');
% axis image
% % set(gca,'YTick',[],'XTick' ,[])
% colormap(flipud(colormap(gray)))
% colorbar

% subplot(1,3,3)
imagesc(log(nUpstreamCells));
title('Upstream Cells No. [Log]');
axis image
% set(gca,'YTick',[],'XTick' ,[])
colormap(flipud(colormap(gray)))
colorbar

figure(8)
set(gcf,'Color',[1 1 1])
imagesc(nUpstreamCells);
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

%% input the initiation and end point for a profile
% Note : X, Y coordinate can be obtained from the Figure 7 of upstream cells
% number.

% initYList = [1468, 3437, 3285, 3385, 3065, 2821, 2834, 1472, 839, 786, 654, 980, 1112, 1070, 1469, 2283];
% initXList = [3186, 3228, 3041, 3027, 3238, 3271, 3287, 2999, 2419, 2294, 2099, 1669, 1429, 1127, 1007, 183];
% endY = 3181; endX = 1951;

% % Imgyecheon, NHi
% initYList = 1468; initXList = 3186;
% endY = 3181; endX = 1951;

% Goljicheon, GJ
initYList = 2834; initXList = 3287;
endY = 1886; endX = 2605;

% % Dongdaecheon, DD
% initYList = 2821; initXList = 3271;
% endY = 2251; endX = 2408;

% % Dongnamcheon, DN
% initYList = 3065; initXList = 3238;
% endY = 2471; endX = 2256;

% % Ocdongcheon, OCD
% initYList = 3437; initXList = 3228;
% endY = 3181; endX = 1953;

% % Deokgucheon, OCDd
% initYList = 3285; initXList = 3041;
% endY = 3181; endX = 1953;

% % Deokgeocheon ??
% initYList = 1112; initXList = 1429;
% endY = 2796; endX = 1219;

% % Naericheon, OCDn
% initYList = 3385; initXList = 3027;
% endY = 3181; endX = 1953;

% % Daegicheon, DG
% initYList = 1472; initXList = 2999;
% % endY = 1528; endX = 2680;
% endY = 1546; endX = 2667; % for remove boundary condition

% % Songcheon, SC
% initYList = 839; initXList = 2419; % reference
% endY = 1885; endX = 2604;

% Songcheon Addition
% initYList = 886; initXList = 2390; % 1
% initYList = 972; initXList = 2516; % 2
% initYList = 1030; initXList = 2376; % 3
% initYList = 1007; initXList = 2676; % 4
% initYList = 917; initXList = 2338; % 5
% initYList = 1416; initXList = 2288; % 6
% endY = 1387; endX = 2556;

% % Songcheon Set
% initYList = [839,886,972,1030,1007,917,1416];
% initXList = [2419,2390,2516,2376,2676,2338,2288];
% endY = 1387; endX = 2556;

% % Odaecheon from Durobong, ODd
% initYList = 654; initXList = 2099;
% endY = 2028; endX = 2340;

% % Odaecheon, Noinbong, ODn
% initYList = 786; initXList = 2294;
% endY = 2028; endX = 2340;

% % Soksacheon, PCs
% initYList = 980; initXList = 1669;
% endY = 2796; endX = 1219;

% % Heungjungcheon, PCh
% initYList = 1070; initXList = 1127;
% endY = 2796; endX = 1219;

% % Jucheongang, JC
% initYList = 1469; initXList = 1007;
% endY = 2989; endX = 1694;

% % Ganglimcheon, JCg
% initYList = 2283; initXList = 183;
% endY = 2989; endX = 1694;

%% analyze stream profile and make map for slope, distance downstream and
% relative slope.

% set the number of considered neighbor
considerNbrForSlope = 3:5:100;
% for the selection among the many smoothed profiles
% note: if you want one smoothed profile in the figure, just 1.
chosenProfile = 1;

% choose range for the calculation of the realtive slope (Rd)
initX_Knick = 4; % local gradient
endX_Knick = 10; % trend gradient

% for the selection of a stream gradient profile for the corrected
% upstream area profiles
% note that large value would be better for area slope curve in log-log
% scale.
chosenSlopeForUpArea_Slope = 10;

% binned area-slope relationship
% make elevation values with fixed vertical interval.
% note: if you want to see subtle change in relatively flat region, lower
% contour interval
contInterval = 1; % contour interval
% determine logarithmic binned average slopes
logBinSize = 0.01;

% define variables
streamGradientMap = nan(mRows,nCols);
distFromInitMap = nan(mRows,nCols);
RdMap = nan(mRows,nCols);
nInit = numel(initYList);

for i=1:nInit
    
    initY = initYList(i);
    initX = initXList(i);

    [upstreamAreaProf,slopePerDist,Rd,streamPath,distFromInit,area_bin,slope_bin] ...
        = AnalyzeStreamProfile(initY,initX,endY,endX ...
        ,chosenProfile,considerNbrForSlope,initX_Knick,endX_Knick ...
        ,chosenSlopeForUpArea_Slope,contInterval,logBinSize...
        ,fSize,orgDEM,nUpstreamCells,mRows,nCols,m2SDSNbrY,m2SDSNbrX,dY,dX);

    nLength = numel(streamPath);
    
    % stream gradient map
    for j=1:nLength
        if isnan(streamGradientMap(streamPath(j))) ...
            || slopePerDist(j,chosenSlopeForUpArea_Slope) < streamGradientMap(streamPath(j))
            streamGradientMap(streamPath(j)) = slopePerDist(j,chosenSlopeForUpArea_Slope);
        end
    end
    
    % distance from init map
    for j=1:nLength
        if isnan(distFromInitMap(streamPath(j))) ...
            || distFromInitMap(streamPath(j)) < distFromInit(j)
            distFromInitMap(streamPath(j)) = distFromInit(j);
        end
    end
    
    
    % relative slope
    for j=1:nLength
        if isnan(RdMap(streamPath(j))) || RdMap(streamPath(j)) < Rd(j)
            RdMap(streamPath(j)) = Rd(j);
        end
    end
    
end

%% Determine the time period during KZ due to orogenic events had passed
% from the outlet of study area to watershed boundary

% Variables
% upstreamAreaProf: upstream area [m^2]
% slopePerDist: stream gradient []
mC = 2; % determine m, exponent of area
nC = 3; % determine n, exponent of slope
uC = 2; % choose uplift rate
mKm = 2; % choose scenario for representative K
sC = 4; % choose slope values

if nC == 1
    n = 1; % plucking
elseif nC == 2
    n = 5/3; % abrasion
elseif nC ==3
    n = 7/3; % cavitation
end

if mC == 1
    m = 0.4*n; % exponent of area
else
    m = 0.6*n;
end

if uC == 1
    upliftRate = 0.00013; % uplift rate [m/yr]
elseif uC == 2 % averaged uplift rate
    upliftRate = 0.0001757; % uplift rate [m/yr]
else
    upliftRate = 0.000269; % uplift rate [m/yr]
end
    
slp = slopePerDist(:,sC); % slope using 7th variable in considerNbrForSlope 

% Determine K

K = upliftRate ./ (upstreamAreaProf.^m .* slp.^n); % [m^(1-2m)/yr]

% distribution of K along the profile
figure(20);

subplot(2,1,1)
plot(distFromInit,K)
xlim([min(distFromInit),max(distFromInit)]);
% ylim([0 1*10^-5])
xlabel('Distance from channel head');
ylabel('K');
title('Distribution of K along the profile')
grid on

subplot(2,1,2)
plot(distFromInit,slp)
xlim([min(distFromInit),max(distFromInit)]);
xlabel('Distance from channel head');
ylabel('Stream gradient');
title('Distribution of stream gradient along the profile')

% statistics of K
negativeMask = slp < 0;
infMask = isinf(K);
artificialArea = false(numel(K),1);
% artificialArea(690:799) = true;

% condition output
fprintf('=============================================\n');
fprintf('Calculation of K conditions: \n');
fprintf('Uplift rate is                     %12.10f.\n',upliftRate);
fprintf('m is                               %12.10f.\n',m);
fprintf('n is                               %12.10f.\n',n);
fprintf('---------------------------------------------\n');

Kmean = mean(K(~negativeMask & ~infMask & ~artificialArea));
fprintf('Mean K is                          %12.10f.\n',Kmean);
Kstd = std(K(~negativeMask & ~infMask & ~artificialArea));
fprintf('Standard deviation of K is         %12.10f.\n',Kstd);

Kmedian = median(K(~negativeMask & ~infMask & ~artificialArea));
fprintf('Median K is                        %12.10f.\n',Kmedian);
Kmad = 1.4826 * mad(K(~negativeMask & ~infMask & ~artificialArea),1); % median absolute deviation
fprintf('Median absolute deviation of K is  %12.10f.\n',Kmad);

% remove outlier
if mKm == 1 %
    outlierMask = K > (Kmean + 2.5*Kstd) | K < -(Kmean + 2.5*Kstd);
    repK_no = mean(K(~(outlierMask | negativeMask | infMask)));
else
    outlierMask = K > (Kmedian + 2.5*Kmad) | K < -(Kmedian + 2.5*Kmad);
    repK_no = median(K(~(outlierMask | negativeMask | infMask)));
end
fprintf('Mean K without outlier is          %12.10f.\n',repK_no);

% Calculate the period for the knickpoint passing through the study area

% variables
nodeLength = zeros(numel(K),1); % length of nodes along a profile [m]
nodeLength(2:end) = distFromInit(2:end) - distFromInit(1:end-1);

% Calculate the time for knickzone passing across a node [yr]
mK = 2; % choose scenario
if mK == 1
    % 1st scenario : using representative K value
    % dt = nodeLength ./ (Kmean .* upstreamAreaProf.^m);
    % dt = nodeLength ./ (Kmean_no .* upstreamAreaProf.^m);
    dt = nodeLength ./ (Kmedian .* upstreamAreaProf.^m);
else
    % 2nd scenario : based on each K value
    K(outlierMask | negativeMask | infMask) = repK_no;
    dt = nodeLength ./ (K .* upstreamAreaProf.^m);
end
dtFromOutlet = flipud(dt);
accDtFromOutlet = cumsum(dtFromOutlet);
accDt = max(accDtFromOutlet);

fprintf('---------------------------------------------\n');
fprintf('Calculation of passing time condition : \n');
if mKm == 1
    fprintf('Input mean K   in NaNs\n');
else
    fprintf('Input median K in NaNs\n');
end
fprintf('---------------------------------------------\n');

fprintf('Estimated time for knickpoint passing is %12.1f.\n',accDt);

figure(21)
distFromOutlet = flipud(max(distFromInit) - distFromInit);
plot(distFromOutlet,accDtFromOutlet);
xlim([min(distFromOutlet),max(distFromOutlet)]);
xlabel('Distance from outlet');
ylabel('Passing time (Year)');
title('Distribution of passing time along the profile')

%% export Rd Map

figure(122); clf;
set(gcf, 'Color',[1,1,1]);

imagesc(RdMap);
set(gca,'DataAspectRatio',[1 1 1]);
colorbar;
title('Relative Stream Gradient Map');

[~,R] = geotiffread('D:\WorkSpace\RawData\SRTM\ver_2\SRTM1\n37_e128_1arc_v3.tif');
% Note that Rd should be multiplied with 10^11 because of the limitation of
% ArcGIS Double Floating data type with the maximum precision 12.
% So after exporting, you should convert the data type of RdMap from doulbe
% to integer, and then convert RdMap into polyline vector using Raster to
% Polyline tool, and add a field for Rd in the attribute table of the
% polyline vector 
geotiffwrite('D:\WorkSpace\Project\NewKnickPoints\Raster\SC_all.tif',RdMap.*10^11,R);


%% statistics of Rd for all profiles

nanSlpMap = isnan(streamGradientMap);
meanSlp = mean(streamGradientMap(~nanSlpMap));
stdSlp = std(streamGradientMap(~nanSlpMap));

oRdMap = RdMap;
nanMRdMap = isnan(oRdMap);
meanRd = mean(oRdMap(~nanMRdMap));
stdRd = std(oRdMap(~nanMRdMap));

mRdMap = RdMap;
mRdMap(RdMap > 10^-5 | RdMap < -10^-5) = nan;
nanMRdMap = isnan(mRdMap);
meanRd = mean(mRdMap(~nanMRdMap));
stdRd = std(mRdMap(~nanMRdMap));

meanRd + 0.5*stdRd
meanRd - 0.5*stdRd
meanRd + 1.5*stdRd
meanRd - 1.5*stdRd
meanRd + 2.5*stdRd
meanRd - 2.5*stdRd

%% regression analysis

% within a specific reach
minArea = 1*10^6;
maxArea = 100*10^6;

% % for not binned
% rangeIdx = find(upstreamAreaProf > minArea & upstreamAreaProf < maxArea);
% polyfit(log(upstreamAreaProf(rangeIdx)),log(slopePerDist(rangeIdx,chosenSlopeForUpArea_Slope)),1)

% for binned
rangeIdx = find(area_bin > minArea & area_bin < maxArea);
% polyfit(log(area_bin(rangeIdx)),log(slope_bin(rangeIdx)),1)
ft = fittype('poly1');
[result,gof] = fit(log(area_bin(rangeIdx))',log(slope_bin(rangeIdx))',ft)

% % for chosen locations
% cX = [3.483*10^7,8.733*10^7,1.764*10^8,5.717*10^8,1.130*10^9];
% cY = [0.04147,0.01987,0.01445,0.005069,0.003747];
% 
% polyfit(log(cX),log(cY),1)

%% make distibution maps

figure(21); clf;
set(gcf, 'Color',[1,1,1]);

imagesc(streamGradientMap);
set(gca,'DataAspectRatio',[1 1 1]);
colorbar;
title('Stream Gradient Map');

figure(22); clf;
set(gcf, 'Color',[1,1,1]);

imagesc(RdMap);
set(gca,'DataAspectRatio',[1 1 1]);
colorbar;
title('Relative Stream Gradient Map');

figure(23); clf;
set(gcf, 'Color',[1,1,1]);

imagesc(distFromInitMap);
set(gca,'DataAspectRatio',[1 1 1]);
colorbar;
title('Distance Downstream');


%% export stream gradient profile to ArcGIS

[~,R] = geotiffread('/Volumes/DATA_BAK/WORKSPACE/RawData/SRTM/ver_2/SRTM1/n37_e128_1arc_v3.tif');
mNUpstreamCells = nUpstreamCells;
mNUpstreamCells(isnan(distFromInitMap)) = nan;
geotiffwrite('/Volumes/DATA_BAK/WORKSPACE/Project/NewKnickPoints/Raster/n37_e128_Area',mNUpstreamCells.*dX.*dY./10^6,R);
geotiffwrite('/Volumes/DATA_BAK/WORKSPACE/Project/NewKnickPoints/Raster/n37_e128_ProfSlp_8',streamGradientMap.*10^8,R);
geotiffwrite('/Volumes/DATA_BAK/WORKSPACE/Project/NewKnickPoints/Raster/n37_e128_Rd_11',RdMap.*10^11,R);
geotiffwrite('/Volumes/DATA_BAK/WORKSPACE/Project/NewKnickPoints/Raster/n37_e128_DownDist',distFromInitMap,R);

