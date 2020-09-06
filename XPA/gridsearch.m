%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid Search
% Masakazu Emoto @ Kobe-University
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ix = gridsearch(x0, xgrid)

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
