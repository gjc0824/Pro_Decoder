function Delta_h = Delta(phi, epsilon)
%   Delta(phi, epsilon) compute the smooth Dirac function
%  
%   created on 04/26/2004
%   author: Chunming Li
%   email: li_chunming@hotmail.com
%   Copyright (c) 2004-2006 by Chunming Li

Delta_h_1 = (0.5 / epsilon) * (1 + cos(pi*phi/epsilon));

Delta_h = Delta_h_1 .* ((phi <= epsilon) & (phi >= -epsilon));
