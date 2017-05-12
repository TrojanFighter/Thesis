function  structure=getStucture(functions,sizes,method,data)

%GETSTRUCTURE - Generate sparsity templates. 
%
% Syntax:  structure=getStucture(functions,sizes,method,data)
%
% Inputs: Defined in transcribeOCP.m
%
% Outputs:
%    structure - Structure containing the sparsity templates
%
% Copyright (C) 2010 Paola Falugi, Eric Kerrigan and Eugene van Wyk. All Rights Reserved.
% This code is published under the BSD License.
% Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) 5 May 2010
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------


disp('Determine sparsity structure');

% Get dimensions
[nt,np,n,m,ng,nb]=deal(sizes{1:6});


% Preallocate memory 
dLdx=ones(1,n);dEdx0=ones(1,n);dEdxf=ones(1,n);dfdx=ones(n);


dLdu=ones(1,m);dfdu=ones(n,m);dEdu0=ones(1,m);dEduf=ones(1,m);


dLdp=zeros(1,np);
dfdp=zeros(n,np);
dEdp=zeros(1,np);
dbdp=zeros(nb,np);


if nt                                      
dLdtf=1;
else
 dLdtf=[];
end;





if ng
 dgdx=ones(ng,n);
 dgdu=ones(ng,m);
 if np
  dgdp=ones(ng,np);
 else
  dgdp=zeros(ng,np);  
 end
else
 dgdx=zeros(ng,n);
 dgdu=zeros(ng,m);
 dgdp=zeros(ng,np);
end
 dgdt=zeros(ng,nt);
 
 if nt&ng                                        
  dgdt(:,1)=ones(ng,nt);
end
 

if nb
  dbdx0=ones(nb,n);
  dbdxf=ones(nb,n);  
  dbdu0=ones(nb,m);
  dbduf=ones(nb,m);
  if np
    dbdp=ones(nb,np);  
  else
    dbdp=zeros(nb,np);  
  end
  if nt
    dbdtf=ones(nb,nt);
  else
    dbdtf=zeros(nb,nt);    
  end   
else
  dbdx0=zeros(nb,n);
  dbdxf=zeros(nb,n);
  dbdu0=zeros(nb,m);
  dbduf=zeros(nb,m);
  dbdp=zeros(nb,np);
  dbdtf=zeros(nb,nt);  
 end
    


if np
  dLdp=ones(1,np);
  dfdp=ones(n,np);
  dEdp=ones(1,np);
else  
  dLdp=zeros(1,np);
  dfdp=zeros(n,np);
  dEdp=zeros(1,np);
end   
    


% 
%  Sturctures for the jacobians 
%---------------------------------------


structure.dfdx=spones(dfdx);
structure.dLdx=spones(dLdx);

structure.dEdx0=spones(dEdx0);
structure.dEdxf=spones(dEdxf);
structure.dfdu=spones(dfdu);


structure.dLdu=spones(dLdu);

structure.dEdu0=spones(dEdu0);
structure.dEduf=spones(dEduf);

  
structure.dfdp=spones(dfdp);structure.dLdp=spones(dLdp);
structure.dEdp=spones(dEdp);




structure.dLdtf=spones(dLdtf);
structure.dgdt=spones(dgdt);
if nt                                      % Get structure of E wrt tf
dEdtf = 1;
else
dEdtf=[];
end;
structure.dEdtf=spones(dEdtf);



structure.dgdx=spones(dgdx);structure.dgdu=spones(dgdu);
structure.dgdp=spones(dgdp);structure.dgdt=spones(dgdt);
structure.dbdx0=spones(dbdx0);structure.dbdxf=spones(dbdxf);
structure.dbdu0=spones(dbdu0);structure.dbduf=spones(dbduf);
structure.dbdp=spones(dbdp);structure.dbdtf=spones(dbdtf);



%------------- END OF CODE --------------

