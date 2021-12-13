function f = fdm3D( u )
% Nonlinear FD model of the PDE
%   function f = fdm3d( u )
% Input: u - current interior solution vector
% Output: f - computed interior residual vector

global m h x y z

v = vectorTOmesh( u, m );
f = zeros((m-2)^3,1);

% fdm approx for interior
for k = 2:m-1
    for j = 2:m-1
        for i = 2:m-1
            indx = index3DMesh( m, i,j,k );
            
            % (u[i+1]-u[i])        
            dx1 = v(index3DMesh(m,i+1,j,k)) - v(indx);
            dx2 = v(index3DMesh(m,i-1,j,k)) - v(indx);
            dy1 = v(index3DMesh(m,i,j+1,k)) - v(indx);
            dy2 = v(index3DMesh(m,i,j-1,k)) - v(indx);
            dz1 = v(index3DMesh(m,i,j,k+1)) - v(indx);
            dz2 = v(index3DMesh(m,i,j,k-1)) - v(indx);

            %  1 + ((u[i+1]+u[i])/2)^2 
            qx1 = 1 + ( 0.5*( v(index3DMesh(m,i+1,j,k)) + v(indx) ) )^2;
            qx2 = 1 + ( 0.5*( v(indx) + v(index3DMesh(m,i-1,j,k)) ) )^2;
            qy1 = 1 + ( 0.5*( v(index3DMesh(m,i,j+1,k)) + v(indx) ) )^2;
            qy2 = 1 + ( 0.5*( v(indx) + v(index3DMesh(m,i,j-1,k)) ) )^2;
            qz1 = 1 + ( 0.5*( v(index3DMesh(m,i,j,k+1)) + v(indx) ) )^2;
            qz2 = 1 + ( 0.5*( v(indx) + v(index3DMesh(m,i,j,k-1)) ) )^2;   
            
            findx = index3DFD( m, i,j,k );
            f(findx) = (1/h^2)*( qx1*dx1 + qx2*dx2 + qy1*dy1 + qy2*dy2 + qz1*dz1 + qz2*dz2 ) ...
                        - rhs(x(i,j,k),y(i,j,k),z(i,j,k));          
        end
    end
end

end
