clc,clear,close all;
%% Generate the Scalar Field
m=25; % The number of radial basis function used to generate the field is m

reg_x=[-5,5]; % These are x and y limite for the region
reg_y=[-5,5];


% Kernels=[(reg_x(2)-reg_x(1))*rand(m,1)+reg_x(1) (reg_y(2)-reg_y(1))*rand(m,1)+reg_y(1)]; % This contains the positions of the RBFs
Kernels=[kron(linspace(reg_x(1),reg_x(2),5)/1.1,[1 1 1 1 1])',repmat(linspace(reg_y(1),reg_y(2),5)/1.1,1,5)'];
% Gamma=2*abs(randn(m,1)); % Gamma in the RBFs
Gamma=0.1*ones(m,1);
% sigma=2*abs(randn(m,1)); % sigma in the RBFs
sigma=2*ones(m,1);
% Theta=2*abs(randn(m,1)); % scaling factors of RBFs
Theta=zeros(m,1);
Theta(19)=5; % This set a peak for the field
% This returns the true scalar value of the whole field
mesh_value = generate_region(reg_x,reg_y,m,Gamma,sigma,Kernels,Theta);

M=size(mesh_value,1);
[mesh_x,mesh_y]=meshgrid(linspace(reg_x(1),reg_x(2),M),flip(linspace(reg_y(1),reg_y(2),M)));
image(linspace(reg_x(1),reg_x(2),M),linspace(reg_y(2),reg_y(1),M),mesh_value,'CDataMapping','scaled');
title("Initialization of agents and fields")
% view(2);
% plot3(mesh_x,mesh_y,255.*mesh_value./max(max(mesh_value)))
% hold on;
% surf(mesh_x,mesh_y,mesh_value,'edgecolor', 'none');
% hold off;
%% Global parameters definition
N=25; % Number of agents
d=0.6; % a parameter
d0=1.5*d; % a parameter
d1=3.5*d; % a parameter
r=4*d; % The communication range of an agent
s=1; % Sampling rate of agents
sigma_w=1; % variance of the sampling noise
%% Agents initialization
t=0;
InitialPosi=-2+randn(N,2); % Initialize the agent positions
InitialPosi(InitialPosi>5)=4; % correct positions outside the border
InitialPosi(InitialPosi<-5)=-4;
Agents=agent.empty(N,0);
for n=1:N
    Agents(n).Code=n;
    Agents(n).Position=InitialPosi(n,:);
    Agents(n).CommuDist=r;
    Agents(n).Speed=0;
    Agents(n).d=d;
    Agents(n).d0=d0;
    Agents(n).d1=d1;
    Agents(n).Theta_est=randn(m,1);
    Agents(n).P=eye(m);
    Agents(n).Gamma=Gamma;
    Agents(n).sigma=sigma;
    Agents(n).Kernels=Kernels;
    Agents(n).v=[0;0];
    Agents(n).gamma=0.2;
    Agents(n).k3=0.1;
    Agents(n).k4=0.1;
    Agents(n).delta_t=0.1;
    Agents(n).BorderLimitX=reg_x;
    Agents(n).BorderLimitY=reg_y;
end

hold on;
scatter(InitialPosi(:,1),InitialPosi(:,2),'g*');
scatter(Kernels(:,1),Kernels(:,2),'ro');
Agents=Agents.UpdateNeighbour;

%% Distributed Learning
IterationFlag=1;
IterNumMax=1000;
IterCount=0;
figure(2),
errorRecord=[];
while IterationFlag
    IterCount=IterCount+1;
    
    Agents=Agents.Measure(Theta,Kernels,Gamma,sigma,sigma_w); % The agents measure their local scalar field values
    Agents=Agents.ReceiveNeighbourMeasurements; % The agents exchange their measurement with neighbours
    Agents=Agents.Learn; % The agents estimate their local field with the measurement
    Agents=Agents.GenerateControl; % The control signal is generated, which will be used to move the agents
    Agents=Agents.Move; % The agents move to new positions according to the control signal
    Agents=Agents.UpdateNeighbour; % agents update their new list of neighbours
    
    for i=1:N % Hard decision that restrict the agents stay inside the region M of interest
        Agents(i).Position(Agents(i).Position>=5)=4.99;
        Agents(i).Position(Agents(i).Position<=-5)=-4.99;
    end
    
    Posi=reshape([Agents.Position],[2,N])';
    
    % plot the current agent positions and the estimated field of agent 1
    if ~mod(IterCount,100)
        figure(2),
        mesh_value_est_1 = generate_region(reg_x,reg_y,m,Gamma,sigma,Kernels,Agents(1).Theta_est);
        image(linspace(reg_x(1),reg_x(2),M),linspace(reg_y(2),reg_y(1),M),mesh_value_est_1,'CDataMapping','scaled');
        hold on;
        scatter(Posi(:,1),Posi(:,2),'g*');
        scatter(Agents(1).Position(1),Agents(1).Position(2),'r*');
        scatter(Kernels(:,1),Kernels(:,2),'ro');
        hold off;
    end
    
    pause(0.001)
    if IterCount>IterNumMax
        IterationFlag=0;
    end
    errorRecord=[errorRecord ReportMSEError(Agents)];
end

figure(2),
mesh_value_est_1 = generate_region(reg_x,reg_y,m,Gamma,sigma,Kernels,Agents(1).Theta_est);
image(linspace(reg_x(1),reg_x(2),M),linspace(reg_y(2),reg_y(1),M),mesh_value_est_1,'CDataMapping','scaled');
hold on;
scatter(Posi(:,1),Posi(:,2),'g*');
scatter(Agents(1).Position(1),Agents(1).Position(2),'r*');
scatter(Kernels(:,1),Kernels(:,2),'ro');
title("Final positions of agents and estimated field of agent-1")
hold off;

figure,semilogy(errorRecord);
title("MSE - Iteration")






