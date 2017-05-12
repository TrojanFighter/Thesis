function [infoNLP,data]=transcribeOCP(problem,guess,options)
%TRANSCRIBEOCP - Process information from 'problem', 'guess' and 'options' for NLP solver
%Specifically:
%Error checking of function definitions and bounds
%Define bounds for NLP variable + continuity, path and boundary constraints
%Format matrices for direct transcription method(if required)
%Generate initial guess for optimization
%Generate structure of the jacobian of the constraints(if required)
%Construct optimal finite-difference pertubation sets(if required)
%
% Syntax:  [infoNLP,data]=transcribeOCP(problem,guess,options)
%
% Inputs:
%    problem - Optimal control problem definition
%    guess   - Guess used to generate starting point for optimization
%    options - Settings input in file settings.m
%
% Outputs:
%    infoNLP - Information required by the NLP solver
%    data - Data passed to the functions evaluated during optimization
%
% Other m-files required: getStructure, getStructureA, getPertubations.
% Subfunctions: checkErrors, transcriptionMatrix
% MAT-files required: none
%
% Copyright (C) 2010 Paola Falugi, Eric Kerrigan and Eugene van Wyk. All Rights Reserved.
% This code is published under the BSD License.
% Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) 5 May 2010
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------
 
% Define (and assign) some parameters
%---------------------------------------


n=length(problem.states.x0l);              % Number of states
m=length(problem.inputs.ul);               % Number of inputs
np=length(problem.parameters.pl);          % Number of free parameters
ng=length(problem.constraints.gl);         % Number of path constraints
nb=length(problem.constraints.bl);         % Number of boundary constraints
N=problem.inputs.N;                        % Number of control actions
M=options.nodes;                           % Number of mesh nodes


if N==0;N=M-1;end % Default N=M-1.

% Multiple-shooting: Number of integration nodes in the interval t=[t0,tf] is N+1
if strcmp(options.transcription,'multiple_shooting')
   M=N+1;
 end




% The final time for discrete time systems is imposed equal to 1 in order 
% to use a unform formulation of the optimization problem
if (strcmp(options.transcription,'discrete'))
   problem.time.tf_min=1;
   problem.time.tf_max=1;
end 
 

% Get bounds for final time and check if time is free or fixed
% nt=0 when the final time is not a variable otherwise nt=1.
tfl=problem.time.tf_min; tfu=problem.time.tf_max;
if tfl==tfu;data.tf=tfl; nt=0;tfl=[];tfu=[];
else nt=1; end


% Set the initial time and other parameters to adjust the temporal
% scale

data.t0=problem.time.t0;
if (strcmp(options.transcription,'discrete'))
     data.k0=problem.time.t0;
     problem.time.t0=0;
     data.Nm=N;
   else
    data.k0=problem.time.t0;
    data.Nm=1;
end

% Set flag if there are more mesh points than integration steps
% The hermite method doubles  the number of integration nodes 
if (strcmp(options.transcription,'hermite'));
   M=2*M-1;ns=2;
else ns=1;
end

data.sizes={nt,np,n,m,ng,nb,M,N,ns};

% Get state and input bounds
xl=problem.states.xl(:); xu=problem.states.xu(:);
ul=problem.inputs.ul(:); uu=problem.inputs.uu(:);

% Choose the smallest constraint set for the terminal bounds
xfl=problem.states.xfl(:);xfu=problem.states.xfu(:);
xf_l=xl;xf_l(xl<xfl)=xfl(xl<xfl);
xf_u=xu;xf_u(xu>xfu)=xfu(xu>xfu);

% Choose the smallest constraint set for the initial state bounds
x0l=problem.states.x0l(:);x0u=problem.states.x0u(:);
x0_l=xl;
x0_l(xl<x0l)=x0l(xl<x0l);
x0_u=xu;
x0_u(xu>x0u )=x0u(xu>x0u);

if isempty(problem.states.x0)
  cx0=0;
  data.x0=x0l;
  data.x0t=(x0_l+(x0_u-x0_l).*rand(n,1));
else  
  cx0=1;
  data.x0t=problem.states.x0.';
  data.x0=problem.states.x0;
end 
data.cx0=cx0;




% Get parameter bounds
if np; pl=problem.parameters.pl(:); pu=problem.parameters.pu(:);
else pl=[];pu=[]; end

% Get bounds for the constraint functions
if ng; gl=problem.constraints.gl(:); gu=problem.constraints.gu(:);
else gl=[];gu=[]; end

if nb; bl=problem.constraints.bl(:);bu=problem.constraints.bu(:);
else bl=[];bu=[]; end


% Store some matrices in data structure

data.options=options;
data.functions=problem.functions;
data.data=problem.data;


if isempty(data.options.perturbation.J)
  data.options.perturbation.J=(eps/2)^(1/3);
end
if isempty(data.options.perturbation.H)    
  data.options.perturbation.H=(8*eps)^(1/3);
end    


% Reference trajectory for the cost function; it is assigned as the set-point
% but the user can be redefine it if necessary

data.references.xr=repmat(problem.setpoints.states,M,1);   
data.references.ur=repmat(problem.setpoints.inputs,M,1);   


% Error checking function definitions
%---------------------------------------
checkErrors(data.sizes,problem.functions,data.data);


% Define bounds for the NLP variable 
%---------------------------------------
% A) direct multiple shooting              : z = [tf p x(0) u(0) x(1) u(1) ... x(M)]'
%
% B) direct transcription (N = M-1) :
%    z = [tf p x(0) u(0) x(1) u(1) ... x(M-1) u(M-1) x(M) ]' 
%    In the trapezoidal method the required u(M) is imposed equal to u(M-1)
%
% C) direct transcription (N < M-1) : z = [tf p x(0) x(1) ...x((M-1)/N-1) u(0) x((M-1)/N)...x(M-1) u(N-1) x(M)]'

switch options.transcription
    case{'multiple_shooting'}
        if N<1; error('Number of control actions incorrect'); end

        nx=n*M;                         % Number of unknown states
        nu=N*m;                         % Number of unknown controls
        nz=nt+np+nx+nu;                 % Length of the optimization variable
        infoNLP.zl=[tfl; pl;  kron(ones(N,1),[xl(:);ul(:)]);xf_l(:)];
        infoNLP.zu=[tfu; pu;  kron(ones(N,1),[xu(:);uu(:)]);xf_u(:)];
        data.Nm=1;
        data.options.ipopt.hessian_approximation='limited-memory';
        

    case{'discrete','euler','trapezoidal','hermite'}
        if mod(M-1,N)||(N>M)              % Check for errors
            error('# integration steps + 1 not divisible or less than N');
        end

        
        nk=M-1;                      % Number of integration steps
        nx=M*n;                      % Number of unknown states  
        nu=N*m;                      % Number of unknown controls
	nz=nt+np+nx+nu;              % Length of the optimization variable

        xpl=repmat(xl(:),(M-1)/N,1);xpu=repmat(xu(:),(M-1)/N,1);
       
        infoNLP.zl=[tfl; pl;repmat([xpl(:);ul(:)],N,1);xf_l(:)];
        infoNLP.zu=[tfu; pu;repmat([xpu(:);uu(:)],N,1);xf_u(:)];
       
  

    otherwise;disp('Unknown method. Check spelling');
end


infoNLP.zl(nt+np+1:nt+np+n)=x0_l;    % Prune bounds for initial conditions
infoNLP.zu(nt+np+1:nt+np+n)=x0_u;

% Define bounds for constraint functions
%---------------------------------------
% Constraints: c = [c0... cF g(0)....g(ng-1) bo, ,bo(nb-1)]' 
%              ck -> defects  for k=0, ..., F=M-1 
%              gk -> general path constraints g(x(k),u(k),p,k) for
%              k=0,...,ng-1
%              bo(k) -> nb boundary conditions b(x(0),x(f),u(0),u(f),p,t)
% 

if ~strcmp(options.transcription,'multiple_shooting')
infoNLP.cl=[kron(ones(M,1),zeros(n,1));kron(ones(M,1),gl(:));bl(:)];
infoNLP.cu=[kron(ones(M,1),zeros(n,1));kron(ones(M,1),gu(:));bu(:)];
else
infoNLP.cl=[kron(ones(M,1),zeros(n,1));kron(ones(M-1,1),gl(:));bl(:)];
infoNLP.cu=[kron(ones(M,1),zeros(n,1));kron(ones(M-1,1),gu(:));bu(:)];   
end
infoNLP.cl(1:n)=-eps;
infoNLP.cu(1:n)=eps;




% Extract sparsity structures
%---------------------------------------




if ~strcmp(options.derivatives,'analytic')
  sparsity=getStructure(problem.functions,data.sizes,options.transcription,data.data);
else
    % The structure of the derivatives is determined only considering a
    % a fixed structure
    sparsity=getStructureA(problem.functions,data.sizes,data);
end
data.sparsity=sparsity; 



if options.tau==0; 
    tau=ns*ones(M-1,1)/(M-1); 
else
   texst=ones(ns,1);
   tau=kron(options.tau,texst);
end
if abs(sum(tau)-ns)>sqrt(eps);error('Time vector (tau) should sum to 1');end
data.tau=tau;






% Format direct transcription matrices
%---------------------------------------
% Format matrices for continuity constraints: c(z)=[A.Vx.z+B.F(z)]
% and generate the quadrature vector w(tau)

if ~strcmp(options.transcription,'multiple_shooting')

    % Generate matrices for transcription method
    data=transcriptionMatrix(options.transcription,nx,nu,tau,data);
    [data.FD.vector,data.FD.index]=getPertubations(sparsity,data.sizes,data);

    % Exploit structure to generate FD vectors
%---------------------------------------


        %
        % Generate structure of Jacobian for the constraints for direct
        % collocations methods
        %---------------------------------------
        %

        Fxu=[kron(speye((M-1)/N),sparsity.dfdx),repmat(sparsity.dfdu,(M-1)/N,1)];
        dfz=[ones((M-1)*n,nt) repmat(sparsity.dfdp,(M-1),1) kron(speye(N),Fxu)];
        dfz=[dfz zeros((M-1)*n,n); ones(n,nt), sparsity.dfdp zeros(n,size(dfz,2)-m-nt-np) sparsity.dfdu sparsity.dfdx];

        Gxu=[kron(speye((M-1)/N),sparse(sparsity.dgdx)),repmat(sparse(sparsity.dgdu),(M-1)/N,1)];
        dgz=[repmat(sparsity.dgdt,M-1,1) repmat(sparsity.dgdp,M-1,1) kron(speye(N),Gxu)];
        [rdgz,cdgz]=size(dgz);
        dgz=[dgz zeros(rdgz,n); sparsity.dgdt sparsity.dgdp zeros(size(sparsity.dgdx,1),cdgz-m-nt-np) sparsity.dgdu sparsity.dgdx];
        A=data.map.A;B=data.map.B;Vx=data.map.Vx;
        
               
        if N==1
        jS= [[zeros(n,nt), zeros(n,np), eye(n), zeros(n,nx+nu-n)]*cx0;...
            A*Vx+B*dfz;...
            dgz;...
            [sparsity.dbdtf sparsity.dbdp sparsity.dbdx0 zeros(nb,(((M-1)/N-1))*n),...
            sparsity.dbdu0 zeros(nb,nx+nu-2*(n+m)-(((M-1)/N-1))*n),...
            sparsity.dbdxf]];
        else
           jS=[[zeros(n,nt), zeros(n,np), eye(n) zeros(n,nx+nu-n)]*cx0;...
             A*Vx+B*dfz;
            dgz;...
            [sparsity.dbdtf sparsity.dbdp sparsity.dbdx0 zeros(nb,(((M-1)/N-1))*n),...
            sparsity.dbdu0 zeros(nb,nx+nu-2*(n+m)-(((M-1)/N-1))*n) sparsity.dbduf,...
            sparsity.dbdxf]]; 
        end
        
        
       data.jacStruct=spones(jS);

       
       %
        % Generate structure of Hessian for the cost function
        %---------------------------------------
        
        Lxu=kron(speye(N),[kron(speye((M-1)/N),sparsity.dLdx),repmat(sparsity.dLdu,(M-1)/N,1)]);
        Lz=[ones(M-1,nt) repmat(sparsity.dLdp,M-1,1) Lxu];
        [rLz,cLz]=size(Lz);
        Lz=[Lz zeros(M-1,n); ones(1,nt) ones(1,np) zeros(1,cLz-m-nt-np) sparsity.dLdu sparsity.dLdx];
        
        Ez=sparse(1,data.FD.index.Ey,1,1,nt+np+n*M+m*N);
        [ib,jb,sb]=find(data.FD.index.b);
        data.costStruct.B=sparse(ib,sb,1,nb,nt+np+n*M+m*N);
        data.costStruct.E=Ez;
        data.costStruct.L=Lz;

        data.hessianStruct=spalloc(M*(n+m),M*(n+m),M*(n+m)*(n+m));
        data.hessianStruct=tril(Lz'*Lz+Ez'*Ez+(data.jacStruct'*data.jacStruct));
        data.funcs.hessianstructure  = @hessianstructure;
        data.funcs.hessian           = @computeHessian;

  else   % if options.transcription,'multiple_shooting'
 
    % Create matrix to extract x vector from z
    Vx=[zeros(nx,nt+np) [kron(speye(N),[speye(n) zeros(n,m)]); zeros(n,(n+m)*N)] [zeros((M-1)*n,n);speye(n)]];
    xV=Vx\speye(n*M);
     

   % Create matrix to extract u vector from z
   Vu=[zeros(nu,nt+np) kron(speye(N),[zeros(m,n) speye(m)])  zeros(nu,n)]; 
   uV=Vu\speye(m*N);
   
   
   data.map.Vu=Vu; 
   data.map.Vx=Vx;     
   data.map.xV=xV;
   data.map.uV=uV;

  
% Exploit structure to generate FD vectors 
%---------------------------------------

    [data.FD.vector,data.FD.index]=getPertubations(sparsity,data.sizes,data);  
    Ez=sparse(1,data.FD.index.Ey,1,1,nt+np+n*M+m*N);
    data.costStruct.E=Ez;
    [ib,jb,sb]=find(data.FD.index.b);
    data.costStruct.B=sparse(ib,sb,1,nb,nt+np+n*M+m*N);
    
   
        %
        % Generate structure of Jacobian for the constraints for the
        % multiple-shooting
        %---------------------------------------
        %
         
       dcdz=[zeros(n,nt), zeros(n,np), eye(n), zeros(n,nx+nu-n)]*cx0;     % Dynamic constraints
       dgdz=[];   % Path Constraints
       dbdz=[]; 
       nm=n+m; 
    
       
       if ng
          dgdz=[dgdz; sparsity.dgdt,  sparsity.dgdp, sparsity.dgdx, sparsity.dgdu, zeros(ng,nm*N-m)];
       end
       
       
       for i=0:M-3 % Compute costs,constraints and sensitivities M-1 times
         k=i*nm+1;  
            dcdz=[dcdz;zeros(n,nt), sparsity.dfdp, zeros(n,k-1), ones(n,nm), -eye(n), zeros(n,(nm)*(M-i-2))];
          if ng
            dgdz=[dgdz; sparsity.dgdt,  sparsity.dgdp,  zeros(ng,k-1+nm), sparsity.dgdx, sparsity.dgdu, zeros(ng,nm*(N-i-1)-m)];
          end
       end
       dcdz=[dcdz;zeros(n,nt), sparsity.dfdp, zeros(n,(M-2)*nm), ones(n,nm), -eye(n)];
       
       

       if nb
         dbdz=[sparsity.dbdtf, sparsity.dbdp, sparsity.dbdx0, sparsity.dbdu0, zeros(nb,nx+nu-2*(n+m)), sparsity.dbduf, sparsity.dbdxf];
         bzz=spones(dbdz'*dbdz);
        else 
         bzz=0;
        end
        
     
        jS=[dcdz;dgdz;dbdz];
        data.jacStruct=spones(jS);
    
 end


% Format initial guess/reference for NLP
%---------------------------------------
infoNLP.z0=zeros(nz,1);
x_guess=zeros(M,n);u_guess=zeros(N,m);


if isempty(guess.states)
   
  x0_lg=x0_l;x0_lg(x0_l==-inf)=-100;
  x0_ug=x0_u;x0_ug(x0_u==inf)=100;
  xf_lg=xf_l;xf_lg(xf_l==-inf)=-100;
  xf_ug=xf_u;xf_ug(xf_u==inf)=100;
  guess.states=[(x0_lg+(x0_ug-x0_lg).*rand(n,1))';(xf_lg+(xf_ug-xf_lg).*rand(n,1))'];    % Prune bounds for initial conditions
end

if isempty(guess.inputs)
   ul_g=ul;ul_g(ul==-inf)=-100;
   uu_g=uu;uu_g(uu==inf)=100;
   guess.inputs=[(ul_g+(uu_g-ul_g).*rand(m,1))';(ul_g+(uu_g-ul_g).*rand(m,1))'];
end    

for i=1:n
  x_guess(:,i)=linspace(guess.states(1,i),guess.states(2,i),M);
end


for i=1:m
 if N>1  
   u_guess(:,i)=linspace(guess.inputs(1,i),guess.inputs(2,i),N);
  else
   u_guess(:,i)=guess.inputs(1,i);
  end
 end

if strcmp(options.transcription,'multiple_shooting')
    
    u_guess=[u_guess;zeros(1,m)];
    infoNLP.z0=reshape([x_guess u_guess]',M*(n+m),1);
    infoNLP.z0(end-m+1:end)=[];
    
else
    u_guess=u_guess';u_guess=u_guess(:);
    x_guess=x_guess';x_guess=x_guess(:);
    infoNLP.z0=data.map.xV*x_guess+data.map.uV*u_guess;
end





if np;infoNLP.z0(nt+1:nt+np)=guess.parameters(:);end
if nt;infoNLP.z0(1)=guess.tf;end





%  Formatting options of the ipopt solver

  data.funcs.objective         = @costFunction;
  data.funcs.gradient          = @costGradient;
  data.funcs.constraints       = @constraintFunction;
  data.funcs.jacobian          = @constraintJacobian;
  %data.funcs.iterfunc          = @callback; 
  data.funcs.jacobianstructure = @jacobianstructure;
  

  data.multipliers=[];


%------------- END OF CODE --------------



function checkErrors(sizes,fcn,data)
%CHECKERRORS - Do some error checking on the ODE and constraint functions
%
% Syntax:  checkErrors(sizes,fcn,method)
%
% Inputs: Defined in main file
%
%------------- BEGIN CODE --------------


% Get functions and dimensions
[f,g,b]=deal(fcn{3:5});[np,n,m,ng,nb]=deal(sizes{2:6});

% Generate some test vectors
x=rand(1,n);u=rand(1,m);p=rand(1,np);


% Check function dimensions
if length(f(x,u,p,rand,data))~=n
    error('Number of states incorrect')
end

if length(g(x,u,p,rand,data))~=ng
    error('Number of path constraints incorrect')
end

if length(b(x(:),x(:),u(:),u(:),p(:),rand,data))~=nb
    error('Number of terminal boundary constraints incorrect')
end
%------------- END OF CODE --------------


function data=transcriptionMatrix(method,nx,nu,tau,data)
%TRANSCRIPTIONMATRIX -  Format matrices for transcription method
%c(z)=A.Vx.z+B.F(z) and generate the quadrature vector w(tau)
%
% Syntax:  data=transcriptionMatrix(method,nx,nu,tau,data)
%
% Inputs: described in main file
%
% Outputs:
%    data - Structure constaining matrices and vectors
%      vx - Matrix that extracts x from (x,u)
%      vu - Matrix that extracts u from (x,u)
%
%------------- BEGIN CODE --------------

% Gets sizes
[nt,np,n,m,ng,nb,M,N]=deal(data.sizes{1:8});

% Preallocate some memory
A=spalloc(nx-n,nx,3*M*n);B=spalloc(nx-n,nx,3*M*n);w=zeros(M,1);
Mi=(M-1)/N;

% Create matrix to extract x vector from z 
Vx=[zeros(nx,nt+np) [kron(speye(N),[speye(n*Mi) zeros(n*Mi,m)]); zeros(n,(n*Mi+m)*N)] [zeros((M-1)*n,n);speye(n)]];
xV=Vx\speye(n*M);

% Create matrix to extract u vector from z
Vu=[zeros(nu,nt+np) kron(speye(N),[zeros(m,n*Mi) speye(m)])  zeros(nu,n)];   
uV=Vu\speye(m*N);

% Create matrices to extract x and u from (x,u) 
vx=Vx;vu=Vu;
if nt; vx(:,1)=[];vu(:,1)=[];           end
if np; vx(:,1:np)=[];vu(:,1:np)=[];     end



% Generate A,B matrices and quadrature vector w
switch method

    case{'discrete'}

        disp('Formatting matrices for the discrete-time system')

        A=spdiags(ones(M*n,1),n,(M-1)*n,M*n);

        B=-spdiags(ones(n*M,1),0,n*(M-1),n*M);

        w=[ones(M-1,1);0];
        W=sparse(1:M,1:M,w);

    case{'euler'}

        disp('Formatting matrices for the euler approximation')
        d1=spdiags(-ones(M*n,1),0,(M-1)*n,M*n);
        d2=spdiags(ones(M*n,1),n,(M-1)*n,M*n);
        A=d1+d2;

       % h=repmat(tau,n,1);
        
        h=kron(tau,ones(n,1));   %PF
        B=-spdiags(h,0,(M-1)*n,n*M);
        
        w=[tau;0];
        W=sparse(1:M,1:M,w);

    case{'trapezoidal'}

        disp('Formatting matrices for the trapezoidal approximation')
        d1=spdiags(-ones(M*n,1),0,(M-1)*n,M*n);
        d2=spdiags(ones(M*n,1),n,(M-1)*n,M*n);
        A=d1+d2;

     

       h1=-0.5*kron(tau,ones(n,1));
       B=spdiags(h1,0,(M-1)*n,n*M)+spdiags(h1,n,(M-1)*n,n*M);


        w=0.5*[tau(1);tau(1:M-2)+tau(2:M-1);tau(M-1)];   % PF
        W=sparse(1:M,1:M,w);                     
        
     
      case{'hermite'}

        disp('Formatting matrices for the hermite-simpson approximation')

        A=spdiags(-0.5*kron(ones(M,1),[ones(n,1);zeros(n,1)]),2*n,A);
        A=spdiags(ones((2*M)*n,1),n,A);
        A=spdiags(-0.5*kron(ones(M,1),[ones(n,1);zeros(n,1)]),0,A);
        A=spdiags(-kron(ones(M-1,1),[zeros(n,1);ones(n,1)]),-n,A);

      
        t0 =-tau(:,ones(1,n)).';ta = t0(:)/8;tb=t0(:)/6;tc=2*t0(:)/3;
        B=spdiags(ta.*repmat([ones(n,1);zeros(n,1)],(M-1)/2,1),0,B);
        B=spdiags(ta.*repmat([ones(n,1);zeros(n,1)],(M-1)/2,1),2*n,B);
        B=spdiags(tb.*repmat([zeros(n,1);ones(n,1)],(M-1)/2,1),-n,B);
        B=spdiags(tb.*repmat([zeros(n,1);ones(n,1)],(M-1)/2,1),n,B);
        B=B+spdiags(tc.*repmat([zeros(n,1);ones(n,1)],(M-1)/2,1),0,nx-n,nx);

        vtau=tau(1:2:M-1);
        w0(1:2:M-2,1)=4*vtau;
        w0(2:2:M-2,1)=vtau(1:end-1)+vtau(2:end);
        
        w=[tau(1);w0;tau(end)]/6;
        W=sparse(1:length(w),1:length(w),w);

        
    otherwise
        disp('Method unknown. Please check spelling')
end

% Store the matrices and return to formatNLP.m
map.A=A;map.B=B;map.w=w;map.W=W;
map.Vx=Vx;map.xV=xV;map.Vu=Vu;map.uV=uV;
data.map=map;

%------------- END OF CODE --------------


