function [error, xx, yy] = Get2DError(fxy, ErrType, Reconstruction, x, y)
%%
%
%
%
%
%
%
%


% Initialization
ErrType = lower(ErrType);

%use meshgrid for plotting
[xx, yy] = meshgrid(x,y);

abs_error = abs(Reconstruction - fxy(xx,yy));

switch ErrType
    case('2norm')
        h = x(2) - x(1);
        error = sqrt(h) * 
    
    case('infinity')
        error = abs_error;
                

        
end

return;
        
        

