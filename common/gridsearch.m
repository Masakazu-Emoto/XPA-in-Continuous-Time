function ix = gridsearch(x0, xgrid)
%% gridserach.m : This code is for finding the index for a value with the given grid
%
%% INPUTS
%    x0         : a value
%    xgrid      : the given grid
%
%% OUTPUTS
%    ix         : the index for x0 with the given grid

    intx = max(size(xgrid));
    ix = 0;

    for jx = 1:intx

        if (x0 < xgrid(jx))
            break;
        end
        ix = ix + 1;
        
    end

    ix = min(max(1,ix),intx - 1);
    
end