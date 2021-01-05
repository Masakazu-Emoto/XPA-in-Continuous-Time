%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "Applying the Explicit Aggregation Algorithm to Heterogeneous Agent Models in Continuous Time."
% By Masakazu Emoto and Takeki Sunakawa
% This code is for finding the grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masakazu EMOTO @ Kobe univerisity 2020/10/22 
% Address : masakazu.emoto@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
