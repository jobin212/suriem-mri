function [error] = GetError(ErrType, abs_error, x)
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
        error = h * norm(abs_error);
    
    case('infinity')
 
        error = max(abs_error);
                
end

return;
        
        

