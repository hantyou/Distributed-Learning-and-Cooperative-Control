function value = func_phi(v,k,Gamma,sigma)
%FUNC_PHI Summary of this function goes here
%   Detailed explanation goes here
v_k=v-k;
v_k2=vecnorm(v_k,2,3).^2;
value=1./Gamma.*exp(-v_k2./(sigma.^2)/2);
end

