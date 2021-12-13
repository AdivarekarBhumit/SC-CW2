function x = directLU(A,b)

[L,U,P] = lu(A);
z = forwardSubstitution(L, P*b, size(b,1));
x = backSubstitution(U,z,size(z,1));
