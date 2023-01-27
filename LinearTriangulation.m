function X = LinearTriangulation(K, C1, R1, C2, R2, x1, x2)
%% LinearTriangulation
% Find 3D positions of the point correspondences using the relative
% position of one camera from another
% Inputs:
%     C1 - size (3 x 1) translation of the first camera pose
%     R1 - size (3 x 3) rotation of the first camera pose
%     C2 - size (3 x 1) translation of the second camera
%     R2 - size (3 x 3) rotation of the second camera pose
%     x1 - size (N x 2) matrix of points in image 1
%     x2 - size (N x 2) matrix of points in image 2, each row corresponding
%       to x1
% Outputs: 
%     X - size (N x 3) matrix whos rows represent the 3D triangulated
%       points


N = size(x1,1);
X = zeros(N,3);
for n=1:N
    %for one point:
    P1 = K * [R1 -1*R1*C1];
    P2 = K * [R2 -1*R2*C2];
    x1x = [0 -1 x1(n,2) ; 1 0 -x1(n,1) ; -x1(n,2) x1(n,1) 0 ];
    x2x = [0 -1 x2(n,2) ; 1 0 -x2(n,1) ; -x2(n,2) x2(n,1) 0 ];
    A = [x1x*P1; x2x*P2];
    [u, s, v] = svd(A);
    X_ = v(:,4);
    X_ = X_ / X_(4);
    X(n,:) = X_(1:3)';
end
