function hessian=computeHessian(z,sigma,lambda,data)

%COMPUTEHESSIAN - Evaluate the Hessian of the Lagrangian at z
%
% Syntax: hessian=computeHessian(z,sigma,lambda,data)
%
% Inputs:
%    z - NLP variable 
%    sigma - Cost function factor
%    lambda - Vector of Lagrange multipliers
%    data - Structure with data needed to evaluate CostFn
%
% Outputs:
%    hessian - Sparse lower triangular segment of the hessian at z 
%
% Other m-files required: multipleShooting.m or directCollocation.m
% Subfunctions: none
% MAT-files required: none
%
% Copyright (C) 2010 Paola Falugi, Eric Kerrigan and Eugene van Wyk. All Rights Reserved.
% This code is published under the BSD License.
% Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) 5 May 2010
% iclocs@imperial.ac.uk


%------------- BEGIN CODE --------------

data.sigma=sigma;data.lambda=lambda;

switch data.options.transcription
    
    case {'multiple_shooting'}
        
        hessian=multipleShooting('hessian',z,data);
   
    otherwise % Direct Transcription Method
        
        hessian=directCollocation('hessian',z,data);

         
end

 %------------- END OF CODE --------------                     
