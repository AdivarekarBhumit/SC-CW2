function [ loc, nf ] = sparseJac( P,n )
% function [ loc, nf ] = sparseJac( P,n )
%   Compute the nf function calls required
%   to build a Jacobian matrix with sparsity 
%   pattern given by P and size nxn
%   Simple greedy algorithm - not optimal
% Input: sparse matrix structure P - same locations as J but matrix of ones
%        full matrix dimension n
% Output: list of required function calls
%         loc(*,1) function call number
%         loc(*,2) equation number
%         nf - number of function calls required

% function index for each dof
loc = zeros(n,2);
loc(1:n,2)=1:n;

nf = 1;
F(:,1)=zeros(n,1);
F(:,1) = F(:,1) + full(P(:,1));
loc(1,1) = 1;
%fprintf(' Column %d in new F() %d\n',1,nf);

for j = 2:n
    found = 0;
    for k = 1:nf
        G = F(:,k)+full(P(:,j));
        if( max(G) == 1 )
            F(:,k) = G;
            loc(j,1) = k;
            %fprintf(' Column %d in F() %d\n',j,k);
            found = k;
            break;
        end
    end
    if( found == 0 )
        nf = nf + 1;
        F(:,nf) = zeros(n,1);
        F(:,nf) = F(:,nf) + full(P(:,j));
        loc(j,1) = nf;
        %fprintf(' Column %d in new F() %d\n',j,nf);
    end
end

% sort function index to dof list for each function call
loc = sortrows(loc,1);

end

