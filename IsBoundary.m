function isTrue = IsBoundary(y,x,Y_TOP_BND,Y_BOTTOM_BND,X_LEFT_BND,X_RIGHT_BND)
%
% function IsBoundary
%
% =========================================================================
%> @section INTRO IsBoundary
%>
%> - 입력되는 셀의 좌표(y,x)가 모형 외곽 경계에 위치하는지를 확인하는 함수
%>
%> @version 0.8
%> @callgraph
%> @callergraph
%>
%> @retval isTrue       : 입력되는 셀의 좌표가 외곽 경계에 위치하는가를 표시하는 변수
%>
%> @param y             : 입력되는 셀의 Y 좌표값
%> @param x             : 입력되는 셀의 X 좌표값
%> @param Y_TOP_BND     : 모형 외곽 위 경계 Y 좌표값
%> @param Y_BOTTOM_BND  : 모형 외곽 아래 경계 Y 좌표값
%> @param X_LEFT_BND    : 모형 외곽 좌 경계 X 좌표값
%> @param X_RIGHT_BND   : 모형 외곽 우 경계 X 좌표값
% =========================================================================

if (x == X_RIGHT_BND) || (x==X_LEFT_BND)
    
    isTrue = true;

elseif (y == Y_BOTTOM_BND) || (y==Y_TOP_BND)
    
    isTrue = true;

else
    isTrue = false;

end