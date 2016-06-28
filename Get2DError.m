function [error, xx, yy] = Get2DError(FncType, ErrType, Reconstruction, x, y)
%%
%
%
%
%
%
%
%


% Initialization
FncType = lower(FncType);
ErrType = lower(ErrType);

%need function
[~, fxy] = Get2DFourierCoefficients(FncType, 1, 1);

%use meshgrid for plotting
[xx, yy] = meshgrid(x,y);

abs_error = abs(Reconstruction - fxy(xx,yy));

switch ErrType
    case('absolute')
        error = abs_error;
    
    case('2norm')
        error = abs_error;
                

        
end

return;
        
        

