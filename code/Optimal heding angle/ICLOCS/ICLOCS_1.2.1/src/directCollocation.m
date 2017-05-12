function solution=directCollocation(required,z,data)
%DIRECTCOLLOCATION - Generate the cost, constraint and gradient
%information for the direct transcription formulation
%
% Syntax:  solution=directCollocation(required,z,data)
%
% Inputs:
%    required - Flag that determines what to compute for current z
%    z - Current NLP variable
%    data - Structure of data required to compute
%
% Outputs:
%    solution - Data structure containing the solution
%
% Other m-files required: gradientFD, jacobianFD, hessian FD (optional)
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


global sol;


% Define some useful variables
[nt,np,n,m,ng,nb,M,N,ns]=deal(data.sizes{1:9});

% Get function definitions
[L,E,f,g,b]=deal(data.functions{:});
vdat=data.data;

% Get matrices
mp=data.map;

%--------------------------------------------------------------------------
% Extract and format vectors from NLP variable
%--------------------------------------------------------------------------

% Extract states and inputs from z and reshape for function evaluations
X=reshape(mp.Vx*z,n,M)';
usp=reshape(mp.Vu*z,m,N)';
U=kron(usp,ones((M-1)/N,1));
U=[U;U(end,:)];


% Extract design parameters if specified and convert to cells
if np; P=reshape(repmat(z(nt+1:nt+np),M,1),np,M)';else P=spalloc(M,0,1);end
% check this line

% Construct time vector
if nt; tf=z(1); else tf=data.tf; end
t0=data.t0;
k0=data.k0;


% if strcmp(data.options.transcription,'discrete'); tf=1; t0=0; end
T=[0;cumsum(data.tau)]*data.Nm/ns;
t=(tf-t0)*T+k0;

% Extract x0,u0,xf,uf,p
p=z(nt+1:nt+np);
x0=z(nt+np+1:nt+np+n);
u0=z(nt+np+(M-1)/N*n+1:nt+np+(M-1)/N*n+m);
xf=z(end-n+1:end);
uf=z(end-m-n+1:end-n);

% Format reference inputs and states if applicable
if ~isempty(data.references.xr);
    Xr=data.references.xr;
else Xr=[];
end

if ~isempty(data.references.ur);
    Ur=data.references.ur;
else Ur=[];
end



%--------------------------------------------------------------------------
% Return the relevant data
%--------------------------------------------------------------------------



switch required
    
    case{'cost'}
    
 %   
   snm=ones(M,1); 
   sol.cost= mp.w'*((tf-t0)*L(X,Xr,U,Ur,P,t,vdat).*snm)+E(x0,xf,u0,uf,p,tf,vdat);
   solution=sol.cost;
   
    
    case{'gradCost'}
    
    [sol.gradCost,sol.JL]=gradientCost(L,X,Xr,U,Ur,P,T,E,x0,xf,u0,uf,p,tf,data);
    solution=sol.gradCost;
      
    case{'const'}

        
    sol.const=[(x0-data.x0t)*data.cx0;mp.A*mp.Vx*z+mp.B*reshape((tf-t0)*f(X,U,P,t,vdat)',M*n,1);
               reshape(g(X,U,P,t,vdat)',M*ng,1);
               b(x0,xf,u0,uf,p,tf,vdat)];
    
      
    solution=sol.const;
        
    case{'jacConst'}
        
        
        if strcmp(data.options.derivatives,'numeric')
            sol.jacConst=jacobianFD(f,g,X,U,P,T,b,x0,xf,u0,uf,p,tf,data);
        end
        
        if strcmp(data.options.derivatives,'analytic')
            [df,dg,db]=jacConst(f,g,X,U,P,t,b,x0,xf,u0,uf,p,tf,t0,data);
            [sol.jacConst,sol.Jf]=jacConstz(df,dg,g,f,X,U,P,T,db,b,x0,xf,u0,uf,p,tf,data);
        end
          
    solution=sol.jacConst;

   
    case{'hessian'}
        
      if strcmp(data.options.derivatives,'analytic')
            [Lzz,Ezz,fzz,gzz,bzz]=hessianAN(L,f,g,sol.Jf,sol.JL,X,U,P,T,E,b,x0,...
                                                    xf,u0,uf,p,tf,data);
               hessc=data.sigma*(Lzz+Ezz)+fzz+gzz+bzz;
               hessc(end-n+1:end,end-n-m+1:end-n)=hessc(end-n+1:end,end-n-m+1:end-n)+hessc(end-n-m+1:end-n,end-n+1:end)';
               sol.hessian=tril(hessc);                         
        else
          if strcmp(data.options.hessianFD,'central')
            sol.hessian=hessianCD(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,tf,data); 
          else
            sol.hessian=hessianFD(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,tf,data);  
          end
        end

    solution=sol.hessian;
    


end

