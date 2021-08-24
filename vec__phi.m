function phi = vec__phi(v,k,Gamma,sigma)
%VEC__PHI Summary of this function goes here
%   Detailed explanation goes here
v_k=v-k;
v_k2=vecnorm(v_k').^2;
phi=1./Gamma.*exp(-v_k2'./(sigma.^2)/2);
end

