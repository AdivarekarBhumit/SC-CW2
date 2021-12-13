
function indx = index3DFD( m, i,j,k )
% function indx = index3DFD( m, i,j,k )
%    Compute vector index on an mxmxm grid
%    with row-based numbering and known
%    boundary data
%    (system is (m-2)^3 interior points)
% Input: FD grid dimension m in each direction
%        (i,j,k) coordinate of the point
% Output: indx - location in the nonlinear system

indx = (k-2)*(m-2)^2 + (j-2)*(m-2) + (i-1);
    
end
