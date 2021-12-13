function v = vectorTOmesh( u, m )
% function v = vectorTOmesh( u, m )
%   Convert a vector on the interior of the m^3 mesh
%   to a vector on the whole mesh
% Input: vector u(1:(m-2)^3)
%        mesh size m in each direction
% Ouput: vector v(1:m^3)

v = zeros(m^3,1);

for k=2:m-1
    for j=2:m-1
        for i=2:m-1
            v(index3DMesh( m, i,j,k )) = u(index3DFD( m, i,j,k ));
        end
    end
end

end
