function phi = vec_Phi(v,k,Gamma,sigma)
%VEC_PHI Summary of this function goes here
%   Detailed explanation goes here
m=size(k,1);
s=size(v,1);
Phi_star=zeros(s,m);
for i=1:s
    Phi_star(i,:)=vec__phi(v(i,:),k,Gamma,sigma)';
end
phi=Phi_star;
end

