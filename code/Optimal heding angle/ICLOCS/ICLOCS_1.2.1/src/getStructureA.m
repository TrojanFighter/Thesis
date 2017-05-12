function  structure=getStuctureA(functions,sizes,data)

%GETSTRUCTUREA - Generate sparsity templates when the analytic option has been selected
%
% Syntax:  structure=getStructureA(functions,sizes,data)
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


% Check the structure of the derivatives for the stage cost L and the terminal cost E

[L,E,f,g,b]=deal(functions{:});

[dL,dE]=gradCost(L,rand(1,n),rand(1,n),rand(1,m),rand(1,m),rand(1,np),0,E,rand(n,1),rand(n,1),rand(m,1),rand(m,1),rand(np,1),rand(nt),data);
[df,dg,db]=jacConst(f,g,rand(1,n),rand(1,m),rand(1,np),0,b,rand(n,1),rand(n,1),rand(m,1),rand(m,1),rand(np,1),rand(nt),data.t0,data);
structure.dE.flag=dE.flag;
structure.dL.flag=dL.flag;
structure.df.flag=df.flag;
structure.dg.flag=dg.flag;
structure.db.flag=db.flag;


if dE.flag==1
if nt&&(~isempty(dE.dtf)); 
  structure.dEdtf = 1;
else
  structure.dEdtf=zeros(1,nt);  
end
if np&&(~isempty(dE.dp)); 
  structure.dEdp=ones(1,np);
else   
  structure.dEdp=zeros(1,np);  
end
if ~isempty(dE.dx0);
  structure.dEdx0=ones(1,n);
else
  structure.dEdx0=zeros(1,n);  
end
if ~isempty(dE.du0)
  structure.dEdu0=ones(1,m);
else
  structure.dEdu0=zeros(1,m);
end
if ~isempty(dE.duf); 
 structure.dEduf=ones(1,m);
else 
  structure.dEduf=zeros(1,m);
end    
if ~isempty(dE.dxf)  
 structure.dEdxf=ones(1,n);
else
 structure.dEdxf=zeros(1,n);
end
else
structure.dEdtf=spones(nt*ones(1,nt));
structure.dEdp=spones(np*ones(1,np));
structure.dEdu0=ones(1,m); 
structure.dEdx0=ones(1,n);
structure.dEduf=ones(1,m);
structure.dEdxf=ones(1,n);
end    


% Preallocate memory

dLdx=ones(1,n);dfdx=ones(n);
dLdu=ones(1,m);dfdu=ones(n,m);
dLdp=zeros(1,np);
dfdp=zeros(n,np);
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
 
 if nt&&ng                                        
  dgdt(:,1)=ones(ng,nt);
 end

if nb
   %%%%%%%%%%%%%%%%%%%%%% 
    
 if db.flag==1
  if nt&&(~isempty(db.dtf)); 
    structure.dbdtf =ones(nb,1);
  else
    structure.dbdtf=zeros(nb,nt);  
  end
  if np&&(~isempty(db.dp)); 
   structure.dbdp=ones(nb,np);
  else   
   structure.dbdp=zeros(nb,np);  
  end
if ~isempty(db.dx0);
  structure.dbdx0=ones(nb,n);
else
  structure.dbdx0=zeros(nb,n);  
end
if ~isempty(db.du0)
  structure.dbdu0=ones(nb,m);
else
  structure.dbdu0=zeros(nb,m);
end
if ~isempty(db.duf); 
 structure.dbduf=ones(nb,m);
else 
  structure.dbduf=zeros(nb,m);
end    
if ~isempty(db.dxf)  
 structure.dbdxf=ones(nb,n);
else
 structure.dbdxf=zeros(nb,n);
end
else
structure.dbdtf=spones(nt*ones(nb,nt));
structure.dbdp=spones(np*ones(nb,np));
structure.dbdu0=ones(nb,m); 
structure.dbdx0=ones(nb,n);
structure.dbduf=ones(nb,m);
structure.dbdxf=ones(nb,n);
end    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  structure.dbdx0=zeros(nb,n);
  structure.dbdxf=zeros(nb,n);
  structure.dbdu0=zeros(nb,m);
  structure.dbduf=zeros(nb,m);
  structure.dbdp=zeros(nb,np);
  structure.dbdtf=zeros(nb,nt);  
end
    


if np
  dLdp=ones(1,np);
  dfdp=ones(n,np);
else  
  dLdp=zeros(1,np);
  dfdp=zeros(n,np);
end   
    


% 
%  Sturctures for the jacobians 
%---------------------------------------


structure.dfdx=spones(dfdx);
structure.dLdx=spones(dLdx);


structure.dfdu=spones(dfdu);
structure.dLdu=spones(dLdu);



  
structure.dfdp=spones(dfdp);structure.dLdp=spones(dLdp);



structure.dLdtf=spones(dLdtf);
structure.dgdt=spones(dgdt);




structure.dgdx=spones(dgdx);structure.dgdu=spones(dgdu);
structure.dgdp=spones(dgdp);structure.dgdt=spones(dgdt);



%------------- END OF CODE --------------

