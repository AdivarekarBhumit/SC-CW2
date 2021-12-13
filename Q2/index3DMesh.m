function indx = index3DMesh( m, i,j,k )
% function indx = index3DMesh( m, i,j,k )
%   Global index for a m^3 3d mesh
% Input: size m in each direction
%        coordinate indices (i,j,k)
% Ouput: vector indx

indx = (k-1)*m^2 + (j-1)*m + (i-1);

end