function phi = EVOLUTION_CV(phi0, mu, nu, lambda, delta_t, epsilon, numIter, g);
%   evolution_withoutedge(I, phi0, mu, nu, lambda_1, lambda_2, delta_t, delta_h, epsilon, numIter);
%   input: 
%       I: input image
%       phi0: level set function to be updated
%       mu: weight for length term
%       nu: weight for area term, default value 0
%       lambda:  weight for 
%       delta_t: time step
%       epsilon: parameter for computing smooth Heaviside and dirac function
%       numIter: number of iterations
%       g: 
%   output: 
%       phi: updated level set function
%  
%   created on 04/26/2004
%   author: Chunming Li
%   email: li_chunming@hotmail.com
%   Copyright (c) 2004-2006 by Chunming Li


 % 镜像边缘延拓
phi = BoundMirrorExpand(phi0)
g = BoundMirrorExpand(g);
for k = 1 : numIter
    phi = BoundMirrorEnsure(phi);
    g = BoundMirrorEnsure(g);
    delta_h = Delta(phi,epsilon);
    % 求那个散度项
    Curv = curvature(phi);
    
    % 拉普拉斯算子
    laplace_phi = del2(phi);
    
    [gx, gy] = gradient(g);
    [grad_phix, grad_phiy] = gradient(phi);

    grad_norm = sqrt(grad_phix.^2 + grad_phiy.^2 + 1e-10);
    grad_phixn = grad_phix./grad_norm;
    grad_phiyn = grad_phiy./grad_norm;

    % updating the phi function
    % div(uF) = udiv(F) + F . grad(u)
    phi=phi+delta_t*(mu*(4*laplace_phi - Curv) + lambda*delta_h.*(g.*Curv+gx.*grad_phixn+gy.*grad_phiyn) + nu*g.*delta_h);    
end
phi = BoundMirrorShrink(phi); % 去掉延拓的边缘

