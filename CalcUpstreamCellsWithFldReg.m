function upstreamCellsNo = CalcUpstreamCellsWithFldReg(DEM,chosenWatershed ...
    ,flood,SDSNbrY,SDSNbrX,floodedRegionIndex,floodedRegionCellsNo)
%
% AccumulateUpstreamFlow
%

% DEM �⺻ �Ӽ�
[mRows,nCols] = size(DEM);

Y = mRows - 2;          % DEM �ܰ� ��踦 ������ Y�� ũ��
X = nCols - 2;          % DEM �ܰ� ��踦 ������ X�� ũ��

% Y_TOP_BND = 1;          % DEM �ܰ� �� ��� Y ��ǥ��
% Y_BOTTOM_BND = mRows;   % DEM �ܰ� �Ʒ� ��� Y ��ǥ��
Y_INI = 2;              % DEM ���� Y ���� ��ǥ��
Y_MAX = Y+1;            % DEM ���� Y ������ ��ǥ��
% 
% X_LEFT_BND = 1;         % DEM �ܰ� �� ��� X ��ǥ��
% X_RIGHT_BND = nCols;    % DEM �ܰ� �� ��� X ��ǥ��
X_INI = 2;              % DEM ���� X ���� ��ǥ��
X_MAX = X+1;            % DEM ���� X ������ ��ǥ��

[sArrayX,sArrayY] = meshgrid(X_INI:X_MAX,Y_INI:Y_MAX);
vectorY = reshape(sArrayY,[],1);
vectorX = reshape(sArrayX,[],1);

% ���� ����
FLOODED = 2;
upstreamCellsNo = zeros(mRows,nCols);

% �� ���� 1��: ��
% ����: flood ���� Ư�� ���� �ܰ��� ������. flooded region�� �Ѳ����� ó����
orgDEM = DEM;
DEM(flood == FLOODED | ~chosenWatershed) = -9999;

vectorDEM = reshape(DEM(Y_INI:Y_MAX,X_INI:X_MAX),[],1);

% �� ���� 2��: ���� ���� ��� outlet ���� �ʰ� ó���Ǿ�� ��
fldRegOutlet = floodedRegionCellsNo > 0;
vectorFldRegOutlet = reshape(fldRegOutlet(Y_INI:Y_MAX,X_INI:X_MAX),[],1);

% �� ���� 3��: ���� ��, ���� outlet�� ��� flooded region ��� ���� ����
% ������ ����
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

% (���� �� ������) ���� ���� ��� ���� ���� �� ������ ����
for ithCell = 1:consideringCellsNo

    % i��° ���� ��ǥ
    ithCellY = sortedDEMYX(ithCell,1);
    ithCellX = sortedDEMYX(ithCell,2);
    
    % ��� �������κ����� ���� �� ������ ������ ���� ���� ���� �й���
    % * ����: ��� �������κ����� ���� �� ������ 1�� ���ϰ�, �̸� ������ ����
    %   ���� ���� �й���
    
    if floodedRegionIndex(ithCellY,ithCellX) >= 0
    
        % ��� �������κ����� ���� �� ������ 1�� ����
        upstreamCellsNo(ithCellY,ithCellX) ...
            = upstreamCellsNo(ithCellY,ithCellX) + 1;
        
    else
        
        % ��� �������κ����� ���� �� ������ 1�� flooded region �� ���� ����
        upstreamCellsNo(ithCellY,ithCellX) ...
            = upstreamCellsNo(ithCellY,ithCellX) + 1 ...
            + floodedRegionCellsNo(ithCellY,ithCellX);
        
    end
    
    % ���� �� ������ ������ ���� ���� ���� �й���
    
    % �Ϸ� ���� ��ǥ�� ���� ������ ����
    downStreamNbrY = SDSNbrY(ithCellY,ithCellX);
    downStreamNbrX = SDSNbrX(ithCellY,ithCellX);
    
    % * ���� : �Ϸ� ���� flooded region�� �ش��Ѵٸ� ��� �������κ����� ����
    %   �� ������ ���ⱸ�� ���� �� ������ �ݿ���
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

% flooded region�� �ش��ϴ� ������ ���� �� ������ ���ⱸ ���� 1�� �۰� ������
% * ����: �̴� �׷������� ū ���̰� ���� �ʵ��� �ϴ� ��ɿ� �Ұ���

% flooded region ���ⱸ ��ǥ�� ������ ����
[tmpOutletY,tmpOutletX] = find(floodedRegionIndex < 0);
floodedRegionsNo = size(tmpOutletY,1);

% �� flooded region�� ���� �� ������ ���Ѵ�.
for ithFloodedRegion = 1:floodedRegionsNo

    %  ���� ó���� flooded region�� ���ⱸ ��ǥ
    outletY = tmpOutletY(ithFloodedRegion,1);
    outletX = tmpOutletX(ithFloodedRegion,1);

    upstreamCellsNo(floodedRegionIndex ...
        == - floodedRegionIndex(outletY,outletX)) ...
        = upstreamCellsNo(outletY,outletX) - 1;

end % ithFloodedRegion = 1:
