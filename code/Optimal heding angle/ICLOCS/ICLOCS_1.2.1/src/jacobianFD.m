% -------------------------------------------------------------------------
%  2. Evaluate the Jacobian of the constraints
% -------------------------------------------------------------------------

function jac=jacobianFD(f,g,X,U,P,T,b,x0,xf,u0,uf,p,tf,data)

% jacobianFD - It evaluates numerically the Jacobian of the constraints
%
% Copyright (C) 2010 Paola Falugi, Eric Kerrigan and Eugene van Wyk. All Rights Reserved.
% This code is published under the BSD License.
% Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) 5 May 2010
% iclocs@imperial.ac.uk



e=data.options.perturbation.J;                                 % pertubation size
[nt,np,n,m,ng,nb,M,N]=deal(data.sizes{1:8});
nz=nt+np+M*n+N*m;
vdat=data.data;
t0=data.t0;
% Compute fz
%------------

fz=spalloc(n*M,nz,n*M*(n+m)+np+nt);

idx=data.FD.index.f;nfd=size(idx,2);
etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;
ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;


for i=1:nfd

    fp=(tf+etf(i)-t0)*f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0).*T+data.k0,vdat);
    fm=(tf-etf(i)-t0)*f(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0).*T+data.k0,vdat);
    fz=fz+sparse(1:M*n,idx(:,i),reshape((fp-fm)'/(2*e),M*n,1),M*n,nz);
end


% Compute gz
%------------

gz=spalloc(ng*M,nz,ng*M*(n+m+np+nt));
if ng

    idx=data.FD.index.g;nfd=size(idx,2);
    etf=e*data.FD.vector.g.etf;ep=data.FD.vector.g.ep;
    ex=data.FD.vector.g.ex;eu=data.FD.vector.g.eu;

   
for i=1:nfd

    gp=g(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(tf+etf(i)-t0).*T+data.k0,vdat);
    gm=g(X-ex{i}*e,U-eu{i}*e,P-ep{i}*e,(tf-etf(i)-t0).*T+data.k0,vdat);
    gz=gz+sparse(1:M*ng,idx(:,i),reshape((gp-gm)'/(2*e),M*ng,1),M*ng,nz);
    
end
end


% Compute bz
%------------
bz=spalloc(nb,nz,(2*m+2*n+nt+np)*nb);
if nb

idx=data.FD.index.b;nfd=size(idx,2);
etf=e*data.FD.vector.b.etf;ep=e*data.FD.vector.b.ep;
ex0=e*data.FD.vector.b.ex0;eu0=e*data.FD.vector.b.eu0;
exf=e*data.FD.vector.b.exf;euf=e*data.FD.vector.b.euf;

for i=1:nfd
    bp=b(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),tf+etf,vdat);
    bm=b(x0-ex0(:,i),xf-exf(:,i),u0-eu0(:,i),uf-euf(:,i),p-ep(:,i),tf-etf,vdat);
    bz=bz+sparse(1:nb,idx(:,i),(bp-bm)/(2*e),nb,nz);
end
end



% Map derivatives to the jacobian
%---------------------------------

jac=[[zeros(n,nt) zeros(n,np) eye(n), zeros(n,(M-1)*n+N*m)]*data.cx0;data.map.A*data.map.Vx+data.map.B*fz;gz;bz];

