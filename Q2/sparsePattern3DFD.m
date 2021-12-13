
function [ A, nz ] = sparsePattern3DFD( n )
% function [ irow,jcol,nz ] = sparsePattern3DFD( n )
%   Set 3D FD sparsity pattern
% Input: mesh dimension n (same in x,y,z dir)
% Output: coordinate format sparse matrix A, number of non-zeros nz

nz = 7*(n-2)^3; % slight over-estimate
irow=zeros(nz,1);
jcol=zeros(nz,1);

nz = 0;

% interior mesh loop to build sparsity
for k=2:n-1
    for j=2:n-1
        for i=2:n-1
            indx = index3DFD( n, i,j,k );
            if( k>2 )
                nz=nz+1; irow(nz) = indx; jcol(nz) = index3DFD( n, i,j,k-1 );
            end
            if( j>2 )
                nz=nz+1; irow(nz) = indx; jcol(nz) = index3DFD( n, i,j-1,k );
            end
            if( i>2 )
                nz=nz+1; irow(nz) = indx; jcol(nz) = index3DFD( n, i-1,j,k );
            end
            nz = nz + 1; irow(nz) = indx; jcol(nz) = indx;
            if( i<n-1 )
                nz=nz+1; irow(nz) = indx; jcol(nz) = index3DFD( n, i+1,j,k );
            end
            if( j<n-1 )
                nz=nz+1; irow(nz) = indx; jcol(nz) = index3DFD( n, i,j+1,k );
            end
            if( k<n-1 )
                nz=nz+1; irow(nz) = indx; jcol(nz) = index3DFD( n, i,j,k+1 );
            end
        end
    end
end

irow = irow(1:nz); % reallocate to correct size
jcol = jcol(1:nz);
mat = ones(nz,1);
A=sparse(irow,jcol,mat);

end

