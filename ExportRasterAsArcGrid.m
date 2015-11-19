function ExportRasterAsArcGrid(outFileName,inputMap,cellSize,xllCorner,yllCorner)
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