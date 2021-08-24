function mesh_value = generate_region(reg_x,reg_y,m,Gamma,sigma,Kernels,theta)
x_min=reg_x(1);
x_max=reg_x(2);
y_min=reg_y(1);
y_max=reg_y(2);
M=512;
x=linspace(x_min,x_max,M);
y=linspace(y_min,y_max,M);
y=flip(y);
[mesh_x,mesh_y]=meshgrid(x,y);
Mesh=zeros(M,M,2);
Mesh(:,:,1)=mesh_x;
Mesh(:,:,2)=mesh_y;
mesh_value=zeros(M,M);
for k=1:m
    K=ones(M,M,2);
    K(:,:,1)=Kernels(k,1)*K(:,:,1);
    K(:,:,2)=Kernels(k,2)*K(:,:,2);
    mesh_value = mesh_value+theta(k).*func_phi(Mesh,K,Gamma(k),sigma(k));
end

end

