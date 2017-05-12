
% -------------------------------------------------------------------------
%  3. Evaluate the Hessian of the Lagrangian
% -------------------------------------------------------------------------

function hessian=hessianFD(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,tf,data)

% hessianFD - It evaluates the Hessian of the Lagrangian with finite differences
%  considering the forward difference formula
%
%
% Syntax:  [Lzz,Ezz,fzz,gzz,bzz]=hessianFD(L,f,g,X,U,P,T,E,b,x0,xf,u0,uf,p,tf,data)
%
% Inputs:
%    see directCollocation.m
%
% Outputs:
%    Lzz - hessian of wL wrt z
%    Ezz - hessian of E wrt z
%    fzz - hessian of f wrt z
%    gzz - hessian of g wrt z
%    bzz - hessian of b wrt z
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Copyright (C) 2010 Paola Falugi, Eric Kerrigan and Eugene van Wyk. All Rights Reserved.
% This code is published under the BSD License.
% Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) 5 May 2010
% iclocs@imperial.ac.uk

%------------- BEGIN CODE ---------------


% Define some useful variables
[nt,np,n,m,ng,nb,M,N]=deal(data.sizes{1:8});
nz=nt+np+M*n+N*m;                           % Length of the primal variable
Xr=data.references.xr;Ur=data.references.ur;
lambda=data.lambda(:);
adjoint_f=reshape(lambda(n+1:n*M)'*data.map.B,n,M)';
adjoint_g=reshape(lambda(n*M+1:n*M+ng*M)',ng,M)';
vdat=data.data;
t0=data.t0;
DT=tf-t0;



e=data.options.perturbation.H;
e2=e*e;                     % Pertubation size


% Compute fzz
% ------------

fzz=spalloc(nz,nz,M*((m+n)*(m+n)+nt+np));


idx=data.FD.index.f;nfd=size(idx,2);
etf=e*data.FD.vector.f.etf;ep=data.FD.vector.f.ep;
ex=data.FD.vector.f.ex;eu=data.FD.vector.f.eu;



 fo=f(X,U,P,DT*T+data.k0,vdat);
for i=1:nfd
    fp1=f(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(DT+etf(i)).*T+data.k0,vdat);
    for j=1:i
      if j==i;fp2=fp1;else
        fp2=f(X+ex{j}*e,U+eu{j}*e,P+ep{j}*e,(DT+etf(j)).*T+data.k0,vdat);
      end

     fpp=f(X+(ex{i}+ex{j})*e,U+(eu{i}+eu{j})*e,P+(ep{i}+ep{j})*e,(DT+etf(i)+etf(j)).*T+data.k0,vdat);
    fte=DT*(fpp-fp2+fo-fp1)+etf(i)*(fpp-fp1)+etf(j)*(fpp-fp2);
    ft=fte.*adjoint_f/e2;
    fzz=fzz+sparse(idx(:,i),idx(:,j),reshape(ft',M*n,1),nz,nz);
   end
end



% Compute gzz
% ------------

gzz=spalloc(nz,nz,M*((m+n)*(m+n)+nt+np));

if ng
idx=data.FD.index.g;nfd=size(idx,2);
etf=e*data.FD.vector.g.etf;ep=data.FD.vector.g.ep;
ex=data.FD.vector.g.ex;eu=data.FD.vector.g.eu;

go=g(X,U,P,DT*T+data.k0,vdat);

for i=1:nfd
    gp1=g(X+ex{i}*e,U+eu{i}*e,P+ep{i}*e,(DT+etf(i)).*T+data.k0,vdat);

    for j=1:i
    
    if j==i;gp2=gp1;else
    gp2=g(X+ex{j}*e,U+eu{j}*e,P+ep{j}*e,(DT+etf(j)).*T+data.k0,vdat);
    end

    gpp=g(X+(ex{i}+ex{j})*e,U+(eu{i}+eu{j})*e,P+(ep{i}+ep{j})*e,...
        (DT+etf(i)+etf(j)).*T+data.k0,vdat);
    
    gt=(gpp-gp2+go-gp1).*adjoint_g/e2; 
    gzz=gzz+sparse(idx(:,i),idx(:,j),reshape(gt',M*ng,1),nz,nz);
    end
end


end







% Compute (w'L)zz
% ----------------

idx=data.FD.index.Ly;
nfd=size(idx,2);                               

etf=data.FD.vector.Ly.etf;ep=data.FD.vector.Ly.ep;
ex=data.FD.vector.Ly.ex;eu=data.FD.vector.Ly.eu;

Lzz=spalloc(nz,nz,M*((m+n)*(m+n+1)/2)+nt+np*np);
Lo=DT*L(X,Xr,U,Ur,P,DT*T+data.k0,vdat);


for i=1:nfd
dt1=e*etf{i};dp1=e*ep{i};dx1=e*ex{i}; du1=e*eu{i};

Lp1=(DT+dt1)*L(X+dx1,Xr,U+du1,Ur,P+dp1,(DT+dt1)*T+data.k0,vdat);

    for j=1:i
    dt2=e*etf{j};dp2=e*ep{j};dx2=e*ex{j}; du2=e*eu{j};

    if j==i;Lp2=Lp1;else
        Lp2=(DT+dt2)*L(X+dx2,Xr,U+du2,Ur,P+dp2,(DT+dt2)*T+data.k0,vdat);
    end

    Lpp=(DT+dt1+dt2)*L(X+dx1+dx2,Xr,U+du1+du2,Ur,P+dp1+dp2,(DT+dt1+dt2)*T+data.k0,vdat);
    Lt=(Lpp-Lp2+Lo-Lp1).*data.map.w/e2; 
    Lzz=Lzz+sparse(idx(:,i),idx(:,j),reshape(Lt',M,1),nz,nz);

    end
end



% Compute Ezz
% ------------

Ezz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));
idx=data.FD.index.Ey;nfd=size(idx,2);                               
etf=e*data.FD.vector.Ey.etf;ep=e*data.FD.vector.Ey.ep;
ex0=e*data.FD.vector.Ey.ex0;eu0=e*data.FD.vector.Ey.eu0;
exf=e*data.FD.vector.Ey.exf;euf=e*data.FD.vector.Ey.euf;

Eo=E(x0,xf,u0,uf,p,tf,vdat);

for i=1:nfd
    Ep1=E(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),tf+etf(:,i),vdat);

    for j=1:i
    
    if j==i;Ep2=Ep1;else
    Ep2=E(x0+ex0(:,j),xf+exf(:,j),u0+eu0(:,j),uf+euf(:,j),p+ep(:,j),tf+etf(:,j),vdat);
    end

    Epp=E(x0+ex0(:,i)+ex0(:,j),xf+exf(:,i)+exf(:,j),u0+eu0(:,i)+eu0(:,j),...
          uf+euf(:,i)+euf(:,j),p+ep(:,i)+ep(:,j),tf+etf(:,i)+etf(:,j),vdat);
    
    Ezz(idx(i),idx(j))=(Epp+Eo-Ep1-Ep2)/e2; 
     
    end
end


% Compute bzz
% ------------
bzz=spalloc(nz,nz,(2*m+2*n+nt+np)*(2*m+2*n+nt+np));
if nb


idx=data.FD.index.b;nfd=size(idx,2);
etf=e*data.FD.vector.b.etf;ep=e*data.FD.vector.b.ep;
ex0=e*data.FD.vector.b.ex0;eu0=e*data.FD.vector.b.eu0;
exf=e*data.FD.vector.b.exf;euf=e*data.FD.vector.b.euf;

adjoint=data.lambda(n*M+M*ng+(~~nb):n*M+M*ng+nb).';
bo=b(x0,xf,u0,uf,p,tf,vdat);

for i=1:nfd
    bp1=b(x0+ex0(:,i),xf+exf(:,i),u0+eu0(:,i),uf+euf(:,i),p+ep(:,i),tf+etf(:,i),vdat);

    for j=1:i
    
    if j==i;bp2=bp1;else
    bp2=b(x0+ex0(:,j),xf+exf(:,j),u0+eu0(:,j),uf+euf(:,j),p+ep(:,j),tf+etf(:,j),vdat);
    end

    bpp=b(x0+ex0(:,i)+ex0(:,j),xf+exf(:,i)+exf(:,j),u0+eu0(:,i)+eu0(:,j),...
          uf+euf(:,i)+euf(:,j),p+ep(:,i)+ep(:,j),tf+etf(:,i)+etf(:,j),vdat);
    bt=(bpp-bp2+bo-bp1).*adjoint'/e2; 
    bzz=bzz+sparse(idx(:,i),idx(:,j),bt,nz,nz);
    end
end

end





% Return the Hessian of the Lagrangian
% -------------------------------------

hessc=data.sigma*(Lzz+Ezz)+fzz+gzz+bzz;
hessc(end-n+1:end,end-n-m+1:end-n)=hessc(end-n+1:end,end-n-m+1:end-n)+hessc(end-n-m+1:end-n,end-n+1:end)';
hessian=tril(hessc);

%------------- END OF CODE --------------

