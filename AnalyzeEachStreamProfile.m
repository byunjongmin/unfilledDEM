function AnalyzeEachStreamProfile(outputDir,initY,initX,endY,endX ...
    ,considerNbrForElev,considerNbrForSlope,nameStrPath)
% @file AnalyzeEachStreamProfile.m
% @brief AnalyzeEachStreamProfile
% @detail Analyze a part of a specific stream profile which begins and ends
%   with given initial and end point. The exact row and column number of
%   initial and end point can be identified in ArcGIS using arrays made by
%   "FindRowColNum.py" script. This function makes stream longitudinal
%   profiles and draws upslope area - slope relationship as well.
%   Additionally it produces the path map of a chosen stream. So you can
%   export the path map into the ArcGIS and identify its location.
%
% @param[in] outputDir name of the directory which stores output data
% @param[in] initY row number of initial point. (e.g. channel head)
% @param[in] initX column number of initial point
% @param[in] endY row number of end point
% @param[in] endX column number of end point. (e.g. channel outlet)
% @param[in] nameStrPath name of exported file
% @param[in] considerNbrForElev number of associated neighbours calculating
%   smoothed elevation
% @param[in] considerNbrForSlope number of associated neighbours
%   calculating slope
%
% @note Version History
%   2013-09-15: Ver. 0.11
%
% @author Jong-Min Byun
%==========================================================================

%% 0. prepare function

% constants
OUTPUT_DIR = '../output';

% prepare output sub-directory
OUTPUT_SUBDIR_PATH = fullfile(OUTPUT_DIR,outputDir);

% load the analyzed data
matFileName = strcat('PostMakeStreamNetwork_',outputDir,'.mat');
matFilePath = fullfile(OUTPUT_SUBDIR_PATH,matFileName);
load(matFilePath);

% load original DEM
[originalDEM,R] = geotiffread(DEMFileName);
baseDEM = originalDEM;

%% 1. identify and export stream path

% identify stream path using given initial and end points coordinates
[streamPath,distFromInit] = RecordStrPath(initY,initX,endY,endX ...
    ,mRows,nCols,SDSNbrY,SDSNbrX,dY,dX);

% export distance from divide map for ArcGIS
intPathDistFromDivideMap = nan(mRows,nCols);

nCells = numel(streamPath);
for ithCell = 1:nCells
   
    intPathDistFromDivideMap(streamPath(ithCell)) = distFromInit(ithCell);    
    
end

IPDFDFileName = strcat(nameStrPath,'_DistFromDivide_',outputDir);
IPDFDFilePath = fullfile(OUTPUT_SUBDIR_PATH,IPDFDFileName);
ExportRasterAsArcGrid(IPDFDFilePath,intPathDistFromDivideMap,R.DeltaX ...
    ,R.XLimWorld(1,1),R.YLimWorld(1,1));

%% 2. longitudinal elevation profile
% draw a stream longitudinal profile on the interesting stream path
% Note: In this procedure, you should find the range affected by artifical
% blocks along chosen stream paths

% stream longitudinal profile using raw DEM
elev = baseDEM(streamPath(:));

% smoothed stream longitudinal profile
smoothedElev = SmoothingElev(considerNbrForElev,distFromInit,baseDEM ...
    ,streamPath);

% draw figure
hf1 = figure(1);
subplot(5,1,1);
plot(distFromInit,elev); hold on
plot(distFromInit,smoothedElev,'m'); hold on

% grid on
title('Longitudinal Stream Profile')
xlabel('Distance From Divide [m]')
ylabel('Elevation [m]')
xlim([0 distFromInit(end)])
ylim([min(elev),max(elev)])


%% 3. longitudinal slope profile

% calculate slope values
slopePerDist = CalcSlope(considerNbrForSlope,distFromInit,smoothedElev);

% prepare figure
subplot(5,1,2)
plot(distFromInit,slopePerDist,'*');

% grid on
xlabel('Distance From Initiaion [m]')
ylabel('Slope')
xlim([0 distFromInit(end)])
ylim([min(slopePerDist),max(slopePerDist)])

%% 4. Draw corrected upstream area profiles on the interesting stream paths

% calculate upstream area
upstreamAreaProf = upstreamCellsNo(streamPath(:)) .* dX .* dY;

% make correction within each depressions
% note: upstream area 가 작은 곳을 만나면 이전 값을 대입할 것

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

% draw figure
subplot(5,1,3);
plot(distFromInit,upstreamAreaProf);

set(gca,'YScale','log');

% grid on
title('Upstream Area Profile')
xlabel('Distance From Divide [m]')
ylabel('Upstream Area [m^2]')
xlim([0 distFromInit(end)])
ylim([min(upstreamAreaProf),max(upstreamAreaProf)])

%% 5. draw upslope area - slope relationship

% prepare figure
subplot(5,1,4)
scatter(upstreamAreaProf,slopePerDist,'Marker','*');

% draw profile
set(gca,'XScale','log','YScale','log');
% grid off
title('Upstream Area - Slope')
xlabel('Upstream Area [m^2]')
ylabel('Slope')
xlim([upstreamAreaProf(1),upstreamAreaProf(end)])

% determine logarithmic binned average slopes

logBinSize = 0.1;

minPower = floor(log10(min(upstreamAreaProf))/logBinSize) * logBinSize;
maxPower = ceil(log10(max(upstreamAreaProf))/logBinSize) * logBinSize;
diffPower = maxPower - minPower;

cellsNo3 = floor(diffPower/logBinSize);

tempArea = zeros(1,cellsNo3);
tempSlope = zeros(1,cellsNo3);

for a = 1:cellsNo3
    
    logAreaMin = minPower + (a-1) * logBinSize;
    logAreaMax = logAreaMin + logBinSize;
    
    tempArea(a) = 10^(logAreaMin + logBinSize * 0.5);
    
    j = find(upstreamAreaProf >= 10^logAreaMin ...
        & upstreamAreaProf < 10^logAreaMax);
    k = size(nonzeros(j));
    if k(1) > 0
        tempSlope(a) = 10^(nanmean(log10(slopePerDist(j))));
    end
end

% extract non-zero elements
lbarea = tempArea(find(tempSlope));
lbslope = tempSlope(find(tempSlope));

subplot(5,1,5)
loglog(lbarea,lbslope,'rs','MarkerFaceColor','r','MarkerSize',3)
% grid on
title('Upstream Area - Slope')
xlabel('Upstream Area [m^2]')
ylabel('Slope')
xlim([min(upstreamAreaProf),max(upstreamAreaProf)])

%% 6. analyze upslope area - slope relationship

AnalyzeSlopeArea(smoothedElev,distFromInit,upstreamAreaProf)


function AnalyzeSlopeArea(streamProfElev,distFromInit,upstreamAreaProf)
%
% function AnalyzeSlopeArea
% 
% 	analyze slope-area relationship using newly calculated slope vaules.
%
% input variables
%
%	streamProfElev: stream profile elevation
%	distFromInit: distance from initial point
%   upstreamAreaProf: upstream area
%	
% Output Variable
%

% calculate new slope values according to fixed vertical interval.
% note: Errors in DEM produce many depressions and negative slope values
%   along a stream profile. Therefore, this algorithm concentrates on
%   making a stream profile without depression. Using this stream profile,
%   it calculates new slope values and log binnined average slope values.
%   However using a smoothed stream proile as an input data makes analysis
%   better. Because the smoothed stream profile has less depressions than
%   raw stream profile data.

% make elevation variable with fixed vertical interval.
% note: if you want to see subtle change in relatively flat region, lower
%   contour interval

contInterval = 1; % contour interval

ipElevMin = 10 * floor(0.1 * streamProfElev(end));
ipElevMax = 10 * ceil(0.1 * streamProfElev(1));
ipElev = ipElevMax:-contInterval:ipElevMin;

% calculate interpolated 'distance from init' value.

tStrProfElev = streamProfElev;
% remove depressions in an input stream profile.
nCells = numel(distFromInit);
for a = nCells:-1:2
    if tStrProfElev(a-1) <= tStrProfElev(a)
        tStrProfElev(a-1) = tStrProfElev(a) + 0.000001;
    end
end

nCells2 = numel(ipElev);
ipDistFromInit = interp1(tStrProfElev,distFromInit,ipElev);

% draw new stream profile with both elevation with fixed vertical interval
% and interpolated distance.

hf10 = figure(10);
subplot(5,1,1)
% for debug
% plot(distFromInit,tStrProfElev,'o',ipDistFromInit,ipElev);
plot(ipDistFromInit,ipElev);

% grid on
title('Longitudinal Stream Profile')
xlabel('Distance From Divide [m]')
ylabel('Elevation [m]')
xlim([0 max(ipDistFromInit)])
ylim([min(ipElev),max(ipElev)])

% calculate new slope values using interpolated elevation and distance data
% note: Basically new slope means local slope rather than 'trend slope'.

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

% draw newly calculated slope vs distance relationship

subplot(5,1,2)
plot(ipDistFromInit,newSlope,'*')

% grid on
xlabel('Distance From Initiaion [m]')
ylabel('Slope')
xlim([0 max(ipDistFromInit)])
ylim([min(newSlope),max(newSlope)])

% draw newly calculated upstream area vs elevation relationship

ipUpstreamArea = interp1(tStrProfElev,upstreamAreaProf,ipElev);

subplot(5,1,3)
plot(ipDistFromInit,ipUpstreamArea);

set(gca,'YScale','log');

% grid on
title('Upstream Area Profile')
xlabel('Distance From Divide [m]')
ylabel('Upstream Area [m^2]')
xlim([0 max(ipDistFromInit)])
ylim([min(ipUpstreamArea),max(ipUpstreamArea)])

% draw newly calculated upstream area vs slope relationship

subplot(5,1,4)
scatter(ipUpstreamArea,newSlope,'m+')

set(gca,'XScale','log','YScale','log');
% grid on
title('Upstream Area - Slope')
xlabel('Upstream Area [m^2]')
ylabel('Slope')
xlim([min(ipUpstreamArea),max(ipUpstreamArea)])

% determine logarithmic binned average slopes

logBinSize = 0.1;

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

subplot(5,1,5)
loglog(lbarea,lbslope,'rs','MarkerFaceColor','r','MarkerSize',3)
% grid on
title('Upstream Area - Slope')
xlabel('Upstream Area [m^2]')
ylabel('Slope')
xlim([min(ipUpstreamArea),max(ipUpstreamArea)])



function slopePerDist = CalcSlope(considerNbr,distFromInit,streamProfElev)
%
% function CalcSlope
% 
%   Calculate slope according to the input window size and base DEM over
%   the defined stream path
%
% Input Variables
%
%   considerNbr: array consisting of number of neighbours
%   distFromInit: distance from initial point
%   streamProfElev: stream profile elevation smoothed according to window sizes
%   
% Output Variable
%
%   slopePerDist
%
%
 
% define variables
nCells = numel(distFromInit);
slopePerDist = zeros(nCells,1);
 
% calculate slope values
 
% remove depressions in an input stream profile.
for a = nCells:-1:2
    if streamProfElev(a-1) <= streamProfElev(a)
        streamProfElev(a-1) = streamProfElev(a) + 0.000001;
    end
end

% first node: forward difference
firstElev = streamProfElev(1);
secondElev = streamProfElev(2);
distBtwCells = distFromInit(2) - distFromInit(1);
 
slopePerDist(1) = (firstElev - secondElev) / distBtwCells;
if slopePerDist(1) <= 0.000001 * 2 / distBtwCells
    slopePerDist(1) = nan;
end
 
% end node: backward difference
firstElev = streamProfElev(nCells-1);
secondElev = streamProfElev(nCells);
distBtwCells = distFromInit(nCells) - distFromInit(nCells-1);
 
slopePerDist(nCells) = (firstElev - secondElev) / distBtwCells;
if slopePerDist(nCells) <= 0.000001 * 2 / distBtwCells
    slopePerDist(nCells) = nan;
end
 
% boundary nodes: central difference with variable ranges
if considerNbr > 1
 
    % upper boundary
    for ithCell = 2:considerNbr
 
        firstEdge = ithCell - (ithCell - 1);
        secondEdge = ithCell + (ithCell - 1);
 
        firstElev = streamProfElev(firstEdge);
        secondElev = streamProfElev(secondEdge); 
 
        distBtwCells = distFromInit(secondEdge) - distFromInit(firstEdge);
        slopePerDist(ithCell) = (firstElev - secondElev) / distBtwCells;
        
        if slopePerDist(ithCell) <= 0.000001 *  (secondEdge - firstEdge) / distBtwCells
            slopePerDist(ithCell) = nan;
        end
 
    end
 
    % lower boundary    
    for ithCell = nCells-considerNbr+1:nCells-1
 
        firstEdge = ithCell - (nCells - ithCell);
        secondEdge = ithCell + (nCells - ithCell);
 
        firstElev = streamProfElev(firstEdge);
        secondElev = streamProfElev(secondEdge); 
 
        distBtwCells = distFromInit(secondEdge) - distFromInit(firstEdge);
        slopePerDist(ithCell) = (firstElev - secondElev) / distBtwCells;
        
        if slopePerDist(ithCell) <= 0.000001 * (secondEdge - firstEdge) / distBtwCells
            slopePerDist(ithCell) = nan;
        end
 
    end
 
end
 
% mid nodes: central difference 
for ithCell = considerNbr+1:nCells-considerNbr
 
    firstEdge = ithCell - considerNbr;
    secondEdge = ithCell + considerNbr;
 
    firstElev = streamProfElev(firstEdge,end);
    secondElev = streamProfElev(secondEdge,end); 
 
    distBtwCells = distFromInit(secondEdge) - distFromInit(firstEdge);
    slopePerDist(ithCell) = (firstElev - secondElev) / distBtwCells;
    
    if slopePerDist(ithCell) <= 0.000001 * 10 * (secondEdge - firstEdge) / distBtwCells
        slopePerDist(ithCell) = nan;
    end
 
end


function ExportRasterAsArcGrid(outFileName,inputMap,cellSize ...
            ,xllCorner,yllCorner)
% 
% function
%   Export array as an ArcGrid ascii type
%
%

% Write (+ overwrite)
fid = fopen([outFileName, '.asc'], 'w');

% Size of input data
[mapRows, mapCols] = size(inputMap);

fprintf(fid, 'ncols         %d\n', mapCols);
fprintf(fid, 'nrows         %d\n', mapRows);
fprintf(fid, 'xllcorner     %14.8f\n', xllCorner);
fprintf(fid, 'yllcorner     %14.8f\n', yllCorner);
fprintf(fid, 'cellsize      %d\n', cellSize);
fprintf(fid, 'NODATA_value  %d\n', -9999);

inputMap(isnan(inputMap)) = -9999;

for iCols = 1 : mapRows
    fprintf(fid, '%8.4f ', inputMap(iCols, :));
    fprintf(fid, '\n');
end

fclose(fid);

function [streamPath,distFromInit] = RecordStrPath(initY,initX ...
    ,endY,endX,mRows,nCols,SDSNbrY,SDSNbrX,dY,dX)
%
% function RecordStrPath
%
% Outline
% 	Identify and record stream path using the input initial and end point 
%   coordinates. Meanwhile it calculates the distance from the inititial
%   point.
%

initIdx = sub2ind([mRows,nCols],initY,initX);
endIdx = sub2ind([mRows,nCols],endY,endX);

% Identify the path determined by the new algorithm
streamPath = initIdx;
distFromInit = 0;
accDist = 0;
ithCellIdx = initIdx;
while ithCellIdx ~= endIdx
    
    % for debugging
    % [tmpY,tmpX] = ind2sub([mRows,nCols],ithCellIdx)
    
    nextCellIdx = sub2ind([mRows,nCols] ...
        ,SDSNbrY(ithCellIdx),SDSNbrX(ithCellIdx));
    
    streamPath = [streamPath; nextCellIdx];
    
    % measure the distance between nodes  
    if abs(ithCellIdx - nextCellIdx) == mRows ...
            || abs(ithCellIdx - nextCellIdx) == 1
        
        accDist = accDist + dY;
        
    else
        
        accDist = accDist + sqrt(dY^2 + dX^2);
    
    end
    
    ithCellIdx = nextCellIdx;

    distFromInit = [distFromInit; accDist];
    
end


function smoothedElev = SmoothingElev(considerNbr,distFromInit,baseDEM ...
    ,streamPath)
%
% function SmoothingElev
% 
% 	Calculate smoothed elevation according to the input window size and
%   base DEM over the defined stream path
%
% Input Variables
%
%	considerNbr: array consisting of number of neighbours
%	distFromInit: distance from initial point
%	baseDEM: base DEM for calcuation
%	streamPath: extracted stream path
%	
% Output Variable
%
%	smoothedElev
%
%

% define variables
nCells = numel(distFromInit);
smoothedElev = zeros(nCells,1);

% calculate smoothed elevation value

% first node
smoothedElev(1) = double(baseDEM(streamPath(1)));

% last node
smoothedElev(nCells) = double(baseDEM(streamPath(nCells)));

% boundary nodes
if considerNbr > 1

    for ithCell = 2:considerNbr

        firstEdge = ithCell - (ithCell - 1);
        secondEdge = ithCell + (ithCell - 1);

        smoothedElev(ithCell) ...
            = median(baseDEM(streamPath(firstEdge:secondEdge)));

    end

    for ithCell = nCells-considerNbr+1:nCells-1

        firstEdge = ithCell - (nCells - ithCell);
        secondEdge = ithCell + (nCells - ithCell);

        smoothedElev(ithCell) ...
            = median(baseDEM(streamPath(firstEdge:secondEdge)));

    end

end

% mid nodes
for ithCell = considerNbr+1:nCells-considerNbr

    firstEdge = ithCell - considerNbr;
    secondEdge = ithCell + considerNbr;

    smoothedElev(ithCell) ...
        = median(baseDEM(streamPath(firstEdge:secondEdge)));

end