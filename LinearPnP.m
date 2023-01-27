function [C, R] = LinearPnP(X, x, K)
%% LinearPnP
% Getting pose from 2D-3D correspondences
% Inputs:
%     X - size (N x 3) matrix of 3D points
%     x - size (N x 2) matrix of 2D points whose rows correspond with X
%     K - size (3 x 3) camera calibration (intrinsics) matrix
% Outputs:
%     C - size (3 x 1) pose transation
%     R - size (3 x 1) pose rotation
%
% IMPORTANT NOTE: While theoretically you can use the x directly when solving
% for the P = [R t] matrix then use the K matrix to correct the error, this is
% more numeically unstable, and thus it is better to calibrate the x values
% before the computation of P then extract R and t directly

N = size(x,1);
K_1 = inv(K);
xK = [x ones(N,1)] * K_1';
A_ = zeros(N*3,12);
%Build the A matrix in blocks
for i=1:N
    Xtd = [X(i,:) 1];
    A= zeros(3,12);
                               A(1,5:8) = -1*Xtd*xK(i,3); A(1,9:12) =    Xtd*xK(i,2);
    A(2,1:4) =    Xtd*xK(i,3);                            A(2,9:12) = -1*Xtd*xK(i,1);
    A(3,1:4) = -1*Xtd*xK(i,2); A(3,5:8) =    Xtd*xK(i,1);
    A_(i*3-2:i*3,:) = A;
end
%find null space of A_
[u,d,v]=svd(A_);
P=v(:,12);
P=reshape(P,[4,3])';

R = P(:,1:3);
t = P(:,4);
[u,d,v]=svd(R);
R = u*v';
t = t/d(1,1);
if (det(u*v') <0)
    R = R*-1;
    t = t*-1;
end
C = -1* R' * t;
