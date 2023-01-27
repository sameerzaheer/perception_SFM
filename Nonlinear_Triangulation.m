function X = Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)
%% Nonlinear_Triangulation
% Refining the poses of the cameras to get a better estimate of the points
% 3D position
% Inputs: 
%     K - size (3 x 3) camera calibration (intrinsics) matrix
%     x
% Outputs: 
%     X - size (N x 3) matrix of refined point 3D locations 
N = size(X0,1);
X = zeros(N,3);
for n=1:N
    X(n,:) = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, ...
        x1(n,:)', x2(n,:)', x3(n,:)', X0(n,:)');
end

end

function X = Single_Point_Nonlinear_Triangulation(K, C1, R1, C2, R2, C3, R3, x1, x2, x3, X0)
    Xprev=X0;
    maxIters = 100;
    for i=1:maxIters
        
        [J1, uvw1] = Jacobian_Triangulation(C1, R1, K, Xprev);
        [J2, uvw2] = Jacobian_Triangulation(C2, R2, K, Xprev);
        [J3, uvw3] = Jacobian_Triangulation(C3, R3, K, Xprev);
        J = [J1; J2; J3];
        b = [x1; x2; x3];
        fX = [uvw1(1)/uvw1(3); uvw1(2)/uvw1(3); ...
              uvw2(1)/uvw2(3); uvw2(2)/uvw2(3); ...
              uvw3(1)/uvw3(3); uvw3(2)/uvw3(3)];
          
        dX = (inv(J' * J)) *J' * (b - fX);
        X = Xprev + dX;
        %error sufficiently low, so break
        if all((X - Xprev) < 0.1*ones(3,1))
            break;
        end
        Xprev = X;
    end    

end


function [J, uvw] = Jacobian_Triangulation(C, R, K, X)
%returns 2x3 matrix.
    uvw = K * R * (X - C);
    
    f = K(1,1); px = K(1,3); py = K(2,3);
    u_ =  [R(1,1)*f+R(3,1)*px R(1,2)*f+R(3,2)*px R(1,3)*f+R(3,3)*px];
    v_ =  [R(2,1)*f+R(3,1)*py R(2,2)*f+R(3,2)*px R(2,3)*f+R(3,3)*px];
    w_ =  [R(3,1) R(3,2) R(3,3)];
    
    J =     (u_ * uvw(3) - w_ * uvw(1))/(uvw(3)^2);
    J = [J; (v_ * uvw(3) - w_ * uvw(2))/(uvw(3)^2)];
  
end
