function F = EstimateFundamentalMatrix(x1, x2)
%% EstimateFundamentalMatrix
% Estimate the fundamental matrix from two image point correspondences 
% Inputs:
%     x1 - size (N x 2) matrix of points in image 1
%     x2 - size (N x 2) matrix of points in image 2, each row corresponding
%       to x1
% Output:
%    F - size (3 x 3) fundamental matrix with rank 2

N= size(x1, 1);
A=zeros(N,9);
for n = 1:N
    A(n,:)= [x1(n,1)*x2(n,1) x1(n,1)*x2(n,2) x1(n,1)...
             x1(n,2)*x2(n,1) x1(n,2)*x2(n,2) x1(n,2)...
             x2(n,1) x2(n,2) 1];
end

[u,s,v]=svd(A);
x = v(:,9);
F= reshape(x,3,3);
%ensure F is rank 2.
[u,s,v]=svd(F);
s(3,3)=0;
F = u * s * v';
F=F/norm(F,'fro');
