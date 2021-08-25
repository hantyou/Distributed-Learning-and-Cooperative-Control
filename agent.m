classdef agent
    
    properties
        Code % agent serial number
        BorderLimitX % X limit of the whole region
        BorderLimitY % Y limit of the whole region
        CommuDist % The farthest distance of communication of between two agents
        Position % The current position of an agent
        Neighbour % record the neighbors of an agents (not including itself, agent code from small to large)
        Speed % This variable is useless, replaced by variable v below 
        SampleNumber % number of samples from self and neighbors
        RealLocalValue % The real value of local scalar field
        Measurement_self % The measured value of an agent itself
        Measurement_around % The collected measurements from the neighbor of an agent
        Ys % The measured values from itself and its neighbors
        K % a variable in the distributed learning stage (field estimation) 
        P % a variable in the distributed learning stage (field estimation) 
        d % a parameter related to the inner potential of the agent flock
        d0 % a parameter related to the inner potential of the agent flock
        d1 % a parameter related to the inner potential of the agent flock
        Phi_star % a variable used in the distributed learning stage (field estimation) 
        Theta_est % the estimated parameters of the field 
        Mu_est % the estimated value of current agent position 
        Gamma % scalar values of the radial basis functions
        sigma % variance of the radial basis functions
        Kernels % kernel positions of the radial basis functions
        v % the control of the agent, can be regarded as its speed
        gamma % a standard adaptive gain sequence, reduces as iteration number goes up
        delta_t % the time between two update of estimation and movement 
        k_di % a scalar of the importance of speed feedback in control stage
        k3 % a scalar of the importance of flock speed in control stage 
        k4 % a scalar of the importance of field gradient in control stage
    end
    
    methods
        function obj = agent(CommuDist,Position)
            if nargin>0
                obj.CommuDist=CommuDist;
                obj.Position=Position;
            end
        end
        
        % function MEASURE gives the current measurement of the agent
        function obj=Measure(obj,theta,Kernel,Gamma,sigma,sigma_w)
            N=length(obj);
            for i=1:N
                %                 RealValue=0;
                Phi_star=vec_Phi(obj(i).Position,Kernel,Gamma,sigma);
                RealValue=Phi_star*theta;
                obj(i).RealLocalValue=RealValue;
                MeasureValue=RealValue+sigma_w*randn;
                obj(i).Measurement_self=MeasureValue;
            end
        end
        
        % function REPORTMSEERROR calculate the current MSE of the agents
        % estimation of the scalar field
        function error=ReportMSEError(obj)
            N=length(obj);
            error=0;
            for i=1:N
                error=error+(obj(i).Mu_est-obj(i).RealLocalValue)^2;
            end
            error =error/N;
        end
        
        % function UPDATENEIGHBOUR update the neighbours of agents
        function obj=UpdateNeighbour(obj)
            N=length(obj);
            for i=1:N
                obj(i).Neighbour=[];
            end
            for i=1:N-1
                for j=i+1:N
                    Dist=norm(obj(i).Position-obj(j).Position);
                    if obj(i).CommuDist>Dist
                        obj(i).Neighbour=[obj(i).Neighbour j];
                        obj(j).Neighbour=[obj(j).Neighbour i];
                    end
                end
            end
            for i=1:N
                obj(i).SampleNumber=1+length(obj(i).Neighbour);
            end
        end
        
        % function RECEIVENEIGHBOURMEASUREMENTS let agents collect
        % measurement values from their neighbours
        function obj=ReceiveNeighbourMeasurements(obj)
            N=length(obj);
            for i=1:N
                ThisAgent=obj(i);
                L=length(ThisAgent.Neighbour);
                s=L+1;
                ThisAgent.Measurement_around=[];
                for j=ThisAgent.Neighbour
                    ThisAgent.Measurement_around=[ThisAgent.Measurement_around obj(j).Measurement_self];
                end
                ThisAgent.Ys=[ThisAgent.Measurement_self ThisAgent.Measurement_around];
                B=sortrows([[ThisAgent.Code;ThisAgent.Neighbour'],ThisAgent.Ys']);
                ThisAgent.Ys=B(:,2);
                obj(i)=ThisAgent;
            end
        end
        
        % function LEARN is the recursive update the estimation of the
        % scalar field
        function obj=Learn(obj)
            N=length(obj);
            for i=1:N
                Y_star=obj(i).Ys;
                kernel=obj(i).Kernels;
                G=obj(i).Gamma;
                s=obj(i).sigma;
                m=size(kernel,1);
                Positions=obj(i).Position;
                for j=obj(i).Neighbour
                    Positions=[Positions;obj(j).Position];
                end
                B=sortrows([[obj(i).Code;obj(i).Neighbour'],Positions]);
                Positions=B(:,2:3);
                Phi_star=vec_Phi(Positions,kernel,G,s);
                % Update
                obj(i).K=obj(i).P*Phi_star'*inv(eye(obj(i).SampleNumber)+Phi_star*obj(i).P*Phi_star');
                obj(i).P=(eye(m)-obj(i).K*Phi_star)*obj(i).P;
                obj(i).Theta_est=obj(i).Theta_est+obj(i).K*(Y_star-Phi_star*obj(i).Theta_est);
                obj(i).Mu_est=vec_Phi(obj(i).Position,kernel,G,s)*obj(i).Theta_est;
            end
            
        end
        
        % function  CALCDIFFPOTENTIALU calculate the potential U used in
        % control stage
        function DiffPotential=CalcDiffPotentialU(obj,obj_self)
            %             obj=obj_all;
            alpha=1;
            beta=0.1;
            DiffPotential=[0,0];
            k1=9;
            k2=0;
            k1=k1/(k1+k2);
            k2=1-k1;
            for j=obj_self.Neighbour
                r_ij=norm(obj_self.Position-obj(j).Position)^2;
                if r_ij<(obj_self.d0^2)
                    DiffPotential=DiffPotential+(r_ij-obj_self.d^2)*(obj_self.Position-obj(j).Position)/(alpha+r_ij)^2;
                else
                    DiffPotential=DiffPotential+func_rho((sqrt(r_ij)-obj_self.d0)/(abs(obj_self.d1-obj_self.d0)))*norm(obj_self.d0^2-obj_self.d^2)/(alpha+obj_self.d0^2)^2*(obj_self.Position-obj(j).Position);
                end
            end
            PotentialInner=DiffPotential;
            order=2;
            PotentialBorder1=-beta/(obj_self.Position(1)-obj_self.BorderLimitX(1))^order+beta/(obj_self.Position(1)-obj_self.BorderLimitX(2))^order;
            PotentialBorder2=-beta/(obj_self.Position(2)-obj_self.BorderLimitY(1))^order+beta/(obj_self.Position(2)-obj_self.BorderLimitY(2))^order;
            PotentialBorder=[PotentialBorder1,PotentialBorder2];
            PotentialBorder(abs(PotentialBorder)>50)=PotentialBorder(abs(PotentialBorder)>10)*0.1;
            DiffPotential=k1*PotentialInner+k2*PotentialBorder;
        end
        
        % function GENRATECONTROL generate the control signal v for the
        % agents
        function obj=GenerateControl(obj)
            N=length(obj);
            k_di=0.5*eye(2);
            delta_t=obj(1).delta_t;
            for i=1:N
                DiffPotential=1*obj.CalcDiffPotentialU(obj(i));
                u_i=-1*DiffPotential'-k_di*(delta_t/obj(i).gamma)*obj(i).v;
                for j=obj(i).Neighbour
                    u_i=u_i+obj(i).k3*(delta_t*(obj(j).v-obj(i).v)/obj(i).gamma);
                end
                u_i=u_i+obj(i).k4*(vec__phi(obj(i).Position,obj(i).Kernels,obj(i).Gamma,obj(i).sigma)./obj(i).sigma.^2.*(-obj(i).Position+obj(i).Kernels))'*obj(i).Theta_est;
                gamma_old=obj(i).gamma;
                obj(i).gamma=obj(i).gamma*1/1.00001;
                obj(i).v=obj(i).gamma/delta_t*(delta_t/gamma_old*obj(i).v+gamma_old*u_i);
            end
        end
        
        % function MOVE let agents move to a new position according to
        % control signal v
        function obj=Move(obj)
            N=length(obj);
            for i=1:N
                obj(i).Position=obj(i).Position+obj(i).delta_t*obj(i).v';
            end
        end
        
    end
end

