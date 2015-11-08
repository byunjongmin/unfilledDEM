function upstreamCellsNo = CalcUpstreamCellsWithFldReg(DEM,chosenWatershed ...
    ,flood,SDSNbrY,SDSNbrX,floodedRegionIndex,floodedRegionCellsNo)
%
% AccumulateUpstreamFlow
%

% DEM 기본 속성
[mRows,nCols] = size(DEM);

Y = mRows - 2;          % DEM 외곽 경계를 제외한 Y축 크기
X = nCols - 2;          % DEM 외곽 경계를 제외한 X축 크기

% Y_TOP_BND = 1;          % DEM 외곽 위 경계 Y 좌표값
% Y_BOTTOM_BND = mRows;   % DEM 외곽 아래 경계 Y 좌표값
Y_INI = 2;              % DEM 영역 Y 시작 좌표값
Y_MAX = Y+1;            % DEM 영역 Y 마지막 좌표값
% 
% X_LEFT_BND = 1;         % DEM 외곽 좌 경계 X 좌표값
% X_RIGHT_BND = nCols;    % DEM 외곽 우 경계 X 좌표값
X_INI = 2;              % DEM 영역 X 시작 좌표값
X_MAX = X+1;            % DEM 영역 X 마지막 좌표값

[sArrayX,sArrayY] = meshgrid(X_INI:X_MAX,Y_INI:Y_MAX);
vectorY = reshape(sArrayY,[],1);
vectorX = reshape(sArrayX,[],1);

% 변수 정의
FLOODED = 2;
upstreamCellsNo = zeros(mRows,nCols);

% 셀 정렬 1차: 고도
% 주의: flood 셀과 특정 유역 외곽은 제외함. flooded region은 한꺼번에 처리함
orgDEM = DEM;
DEM(flood == FLOODED | ~chosenWatershed) = -9999;

vectorDEM = reshape(DEM(Y_INI:Y_MAX,X_INI:X_MAX),[],1);

% 셀 정렬 2차: 같은 고도일 경우 outlet 가장 늦게 처리되어야 함
fldRegOutlet = floodedRegionCellsNo > 0;
vectorFldRegOutlet = reshape(fldRegOutlet(Y_INI:Y_MAX,X_INI:X_MAX),[],1);

% 셀 정렬 3차: 같은 고도, 같은 outlet일 경우 flooded region 평균 고도가 높은
% 순서로 정렬
outletIdx = find(floodedRegionIndex < 0);
meanElevFldReg = zeros(mRows,nCols);
for ithOutlet = 1:numel(outletIdx)
    
    ithFldRegID = -floodedRegionIndex(outletIdx(ithOutlet));
    ithFldRegIdx = floodedRegionIndex == ithFldRegID;
    meanElevFldReg(outletIdx(ithOutlet)) = mean(orgDEM(ithFldRegIdx));
    
end
vectorMeanElevFldReg = reshape(meanElevFldReg(Y_INI:Y_MAX,X_INI:X_MAX),[],1);

sortedDEMYX = [vectorY,vectorX,vectorDEM,vectorFldRegOutlet,vectorMeanElevFldReg];

sortedDEMYX = sortrows(sortedDEMYX,[-3,4,-5]);

consideringCellsNo = find(vectorDEM > -9999);

consideringCellsNo = size(consideringCellsNo,1);

% (높은 고도 순으로) 개별 셀의 상부 유역 누적 셀 개수를 구함
for ithCell = 1:consideringCellsNo

    % i번째 셀의 좌표
    ithCellY = sortedDEMYX(ithCell,1);
    ithCellX = sortedDEMYX(ithCell,2);
    
    % 상부 유역으로부터의 누적 셀 개수를 유향을 따라 다음 셀에 분배함
    % * 개요: 상부 유역으로부터의 누적 셀 개수에 1을 더하고, 이를 유향을 따라
    %   다음 셀에 분배함
    
    if floodedRegionIndex(ithCellY,ithCellX) >= 0
    
        % 상부 유역으로부터의 누적 셀 개수에 1을 더함
        upstreamCellsNo(ithCellY,ithCellX) ...
            = upstreamCellsNo(ithCellY,ithCellX) + 1;
        
    else
        
        % 상부 유역으로부터의 누적 셀 개수에 1과 flooded region 셀 수를 더함
        upstreamCellsNo(ithCellY,ithCellX) ...
            = upstreamCellsNo(ithCellY,ithCellX) + 1 ...
            + floodedRegionCellsNo(ithCellY,ithCellX);
        
    end
    
    % 누적 셀 개수를 유향을 따라 다음 셀에 분배함
    
    % 하류 셀의 좌표와 선형 색인을 구함
    downStreamNbrY = SDSNbrY(ithCellY,ithCellX);
    downStreamNbrX = SDSNbrX(ithCellY,ithCellX);
    
    % * 주의 : 하류 셀이 flooded region에 해당한다면 상부 유역으로부터의 누적
    %   셀 개수를 유출구의 누적 셀 개수에 반영함
    if flood(downStreamNbrY,downStreamNbrX) == FLOODED
        
       outletY = SDSNbrY(downStreamNbrY,downStreamNbrX);
       outletX = SDSNbrX(downStreamNbrY,downStreamNbrX);
       
       upstreamCellsNo(outletY,outletX) ...
        = upstreamCellsNo(outletY,outletX) ...
        + upstreamCellsNo(ithCellY,ithCellX);
        
    else
    
        upstreamCellsNo(downStreamNbrY,downStreamNbrX) ...
            = upstreamCellsNo(downStreamNbrY,downStreamNbrX) ...
            + upstreamCellsNo(ithCellY,ithCellX);
    
    end
    
end

% flooded region에 해당하는 셀들의 누적 셀 개수를 유출구 보다 1개 작게 설정함
% * 주의: 이는 그래프에서 큰 차이가 나지 않도록 하는 요령에 불과함

% flooded region 유출구 좌표와 개수를 구함
[tmpOutletY,tmpOutletX] = find(floodedRegionIndex < 0);
floodedRegionsNo = size(tmpOutletY,1);

% 각 flooded region의 누적 셀 개수를 구한다.
for ithFloodedRegion = 1:floodedRegionsNo

    %  현재 처리할 flooded region의 유출구 좌표
    outletY = tmpOutletY(ithFloodedRegion,1);
    outletX = tmpOutletX(ithFloodedRegion,1);

    upstreamCellsNo(floodedRegionIndex ...
        == - floodedRegionIndex(outletY,outletX)) ...
        = upstreamCellsNo(outletY,outletX) - 1;

end % ithFloodedRegion = 1:
