function [upstreamAreaProf,slopePerDist,Rd,streamPath,distFromInit] ...
    = AnalyzeStreamProfile(initY,initX,endY,endX ...
    ,chosenProfile,considerNbrForSlope,initX_Knick,endX_Knick ...
    ,chosenSlopeForUpArea_Slope,contInterval,logBinSize...
    ,fSize,orgDEM,nUpstreamCells,mRows,nCols,m2SDSNbrY,m2SDSNbrX,dY,dX)
%
% function AnalyzeStreamProfile
% 
% 	Draw longitudinal profile, calculate relative slope according to given
%   window sizes, draw area-slope relationship and log-binned relationship
%
%	
% Output Variable
%
%	upstreamAreaProf
%   slopePerDist
%
% Version: 0.1.0 / 2016-01-05
%

%% identify stream path using given initial and end point coordinates
[streamPath,distFromInit] = RecordStrPath(initY,initX,endY,endX ...
    ,mRows,nCols,m2SDSNbrY,m2SDSNbrX,dY,dX);

% stream longitudinal profile using raw DEM
% note that we use the original DEM
profElev = orgDEM(streamPath(:));
nCells = numel(streamPath);

% draw figure
figure(10); clf;
set(gcf, 'Color',[1,1,1]);

subplot(4,1,1)
stairs(distFromInit,profElev);

% grid on
set(gca,'FontSize',13,'fontWeight','bold')
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

smoothedProfElev = SmoothingElev(considerNbrForProf,distFromInit,orgDEM,streamPath);

subplot(4,1,2)

cc = jet(nCNbr);
for ithLine = 1:nCNbr
    plot(distFromInit,smoothedProfElev(:,ithLine),'color',cc(ithLine,:));
    hold on
end

grid on
set(gca,'FontSize',13,'fontWeight','bold')
title('Smoothed Longitudinal Stream Profile')
xlabel('Distance From Divide [m]')
ylabel('Elevation [m]')
xlim([0 distFromInit(end)])
ylim([min(profElev),max(profElev)])

%% Draw stream gradient

% chose base elevation of profile
inputProfElev = smoothedProfElev;
nCNbrForSlope = numel(considerNbrForSlope);
% 너무 크게 구간을 잡으면 연산이 안됨.
[slopePerDist,slopePerDist_TF] ...
    = CalcSlope(considerNbrForSlope,distFromInit,inputProfElev,chosenProfile);

subplot(4,1,3)

cc = jet(nCNbrForSlope);
for ithLine = 1:nCNbrForSlope
    RdIdx = slopePerDist_TF(:,ithLine) == true;
    plot(distFromInit(RdIdx),slopePerDist(RdIdx,ithLine),'color',cc(ithLine,:));
    hold on
end

grid on
set(gca,'FontSize',13,'fontWeight','bold')
title('Stream Gradient')
xlabel('Distance From Initiaion [m]')
ylabel('Slope')
xlim([0 distFromInit(end)])
ylim([min(slopePerDist(:,1)),max(slopePerDist(:,1))])

%% Identify Knickpoints: Relative Slopes (Rd)

Rd = nan(nCells,1); % relative slope
rangeX = initX_Knick:endX_Knick;
for ithCell = 1:nCells
    
    rXIdx = find(slopePerDist(ithCell,slopePerDist_TF(ithCell,rangeX) == true));
    if numel(rXIdx) > 1
        
        y = slopePerDist(ithCell,rangeX(rXIdx));
        x = considerNbrForSlope(rangeX(rXIdx)) .* 2 .* double(dX);
        c = polyfit(x,y,1);
    
        Rd(ithCell,1) = -c(1);
    end

end

subplot(4,1,4)

plot(distFromInit,Rd);

grid on
set(gca,'FontSize',13,'fontWeight','bold')
title('Relative Slope')
xlabel('Distance From Divide [m]')
ylabel('Relative Slope [m^-1]')
xlim([0 distFromInit(end)])
ylim([-1E-5,1E-5])


% simple statistics of slopePerDist
slopePerDistMin = zeros(1,nCNbrForSlope);
slopePerDistMax = zeros(1,nCNbrForSlope);
slopePerDistMean = zeros(1,nCNbrForSlope);
slopePerDistStd = zeros(1,nCNbrForSlope);

for i=1:nCNbrForSlope
    
    RdIdx = slopePerDist_TF(:,i) == true;
    slopePerDistMin(1,i) = min(slopePerDist(RdIdx,i));
    slopePerDistMax(1,i) = max(slopePerDist(RdIdx,i));
    slopePerDistMean(1,i) = mean(slopePerDist(RdIdx,i));
    slopePerDistStd(1,i) = std(slopePerDist(RdIdx,i));
    
end

% slopePerDistTotal = [slopePerDistTotal...
%     ;slopePerDistMin; slopePerDistMax; slopePerDistMean; slopePerDistStd];

figure(11); clf;

subplot(2,2,1)
scatter(considerNbrForSlope .* dX,slopePerDistMin);

grid on
set(gca,'FontSize',13,'fontWeight','bold')
title('Change in Minimum Stream Gradient')
xlabel('Distance (Dd)')
ylabel('Stream Gradient')

subplot(2,2,2)
scatter(considerNbrForSlope .* dX,slopePerDistMax);

grid on
set(gca,'FontSize',13,'fontWeight','bold')
title('Change in Maximum Stream Gradient')
xlabel('Distance (Dd)')
ylabel('Stream Gradient')

subplot(2,2,3)
scatter(considerNbrForSlope .* dX,slopePerDistMean);

grid on
set(gca,'FontSize',13,'fontWeight','bold')
title('Change in Mean Stream Gradient')
xlabel('Distance (Dd)')
ylabel('Stream Gradient')

subplot(2,2,4)
scatter(considerNbrForSlope .* dX,slopePerDistStd);

grid on
set(gca,'FontSize',13,'fontWeight','bold')
title('Change in Standard Deviation Stream Gradient')
xlabel('Distance (Dd)')
ylabel('Standard Deviation')


%% draw corrected upstream area profiles on the interesting stream paths

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
title('Upstream Contributing Area - Distance')
xlabel('Distance From Divide [m]')
ylabel('Area [m^2]')
xlim([0 distFromInit(end)])
ylim([min(upstreamAreaProf),max(upstreamAreaProf)])
set(gca,'FontSize',18);

% draw upslope area - slope relationship

subplot(2,1,2)
scatter(upstreamAreaProf,slopePerDist(:,chosenSlopeForUpArea_Slope) ...
    ,'Marker','*');

set(gca,'XScale','log','YScale','log');
grid off
title('Upstream Contributing Area - Slope')
xlabel('Area [m^2]')
ylabel('Slope')
xlim([upstreamAreaProf(1),upstreamAreaProf(end)])
set(gca,'FontSize',18);

%% binned area-slope relationship: analyze slope-area relationship
% according to fixed vertical interval
% note: errors in DEM produce many depressions and negative slope values
% along a stream profile. Therefore, this algorithm concentrates on
% making a stream profile without depression. Using this stream profile,
% it calculates new slope values and log binnined average slope values.
% However using a smoothed stream proile as an input data makes analysis
% better. Because the smoothed stream profile has less depressions than
% raw stream profile data.

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

figure(13); clf;

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
area_bin = tempArea(find(tempSlope));
slope_bin = tempSlope(find(tempSlope));

subplot(3,1,3)
loglog(area_bin,slope_bin,'rs','MarkerFaceColor','r','MarkerSize',3)
% grid on
title('Upstream Area - Slope')
xlabel('Upstream Area [m^2]')
ylabel('Slope')
xlim([min(ipUpstreamArea),max(ipUpstreamArea)])