function rho = func_rho(z)

h=0.8;

rho=zeros(size(z,1),size(z,2));
rho(z>=0&z<h)=rho((z>=0&z<h))+1;
rho(z>=h&z<=1)=0.5.*(1+cos(pi.*(z(z>=h&z<=1)-h)./(1-h)));
end

