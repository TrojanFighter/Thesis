function [a] =  trapacc(~,t1,t2,t3,amax,T)

%%Parameters of the trapezoid function
%amax=max accelartion
%t0=beginning of the ramp
%t1=end of the ramp & beginning of the const. acc.
%t2=end of the const. acc. & beginning of the decent
%t3=end of the decent

T=mod(T,2*t3);

if T<t1
    a= (amax/t1)*T;
elseif T<t2
    a= amax;
elseif T<2*t3-t2
    a= (amax/(t2-t3))*(T-t3);
elseif T<2*t3-t1
    a= -amax;
elseif T<2*t3
    a= amax *(T-2*t3);
    
elseif 2*t3<T
    error('time value out of range');
    
end


