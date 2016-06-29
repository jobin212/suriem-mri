function [error] = Get2DError(ErrType, abs_error, x, y)
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

error = 0;

switch ErrType
    case('2norm')
        h = x(2) - x(1);
        error = sqrt(h) * norm(abs_error(:));
    
    case('infinity')
 
        error = max(max(abs_error));
                
end

return;
        
        

