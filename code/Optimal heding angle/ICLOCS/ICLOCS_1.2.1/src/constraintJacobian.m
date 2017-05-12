function jacConst=constraintJacobian(z,data)

%CONSTRAINTJACOBIAN - Evaluate the Jacobian of the constraints
%
% Syntax: jacConst=constraintJacobian(z,data)
%
% Inputs:
%    z - NLP variable 
%    data - Structure with data needed to evaluate CostFn
%
% Outputs:
%    jacConst - Sparse matrix with jacobian of the constraints at z
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
        
         jacConst=multipleShooting('jacConst',z,data);
   
        
  
    otherwise % Direct Transcription Method
        
      
     jacConst=directCollocation('jacConst',z,data);
            
         
         
end

 %------------- END OF CODE --------------                     






