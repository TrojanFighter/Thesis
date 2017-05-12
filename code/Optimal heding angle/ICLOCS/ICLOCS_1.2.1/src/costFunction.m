function J=costFunction(z,data)
%COSTFUNCTION - Evaluate the cost function
%
% Syntax: J=costFunction(z,data)
%
% Inputs:
%    z    - NLP variable 
%    data - Structure with data needed to evaluate CostFn
%
% Outputs:
%    J - Scalar cost evaluated at z
%
% Other m-files required: multipleShooting.m or directTranscription.m
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

switch data.options.transcription
    
    case {'multiple_shooting'}
        
        J=multipleShooting('cost',z,data);
      
    % Direct Collocation
    case{'discrete','euler','trapezoidal','hermite'}; 
        
        J=directCollocation('cost',z,data);
    
    % Pseudospectral Collocation

end

%------------- END OF CODE --------------