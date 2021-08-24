classdef agent
    %AGENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Code
        CommuDist
        Position
        Neighbour
        Speed
        SampleNumber
        RealLocalValue
        Measurement_self
        Measurement_around
        K
        P
        d
        d0
        d1
        Phi_star
        Theta_est
        Mu_est
        Ys
        Gamma
        sigma
        Kernels
        v
        gamma
        delta_t
        k_di
        k3
        k4
    end
    
    methods
        function obj = agent(CommuDist,Position)
            if nargin>0
                obj.CommuDist=CommuDist;
                obj.Position=Position;
            end
        end
        
        function obj=Measure(obj,theta,Kernel,Gamma,sigma,sigma_w)
            N=length(obj);
            for i=1:N
                %                 RealValue=0;
                Phi_star=vec_Phi(obj(i).Position,Kernel,Gamma,sigma);
                RealValue=Phi_star*theta;
                MeasureValue=RealValue+sigma_w*randn;
                obj(i).Measurement_self=MeasureValue;
            end
        end
        
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
        
        function DiffPotential=CalcDiffPotentialU(obj,obj_self)
            %             obj=obj_all;
            alpha=1;
            DiffPotential=[0,0];
            for j=obj_self.Neighbour
                r_ij=norm(obj_self.Position-obj(j).Position)^2;
                if r_ij<(obj_self.d0^2)
                    DiffPotential=DiffPotential+(r_ij-obj_self.d^2)*(obj_self.Position-obj(j).Position)/(alpha+r_ij)^2;
                else
                    DiffPotential=DiffPotential+func_rho((sqrt(r_ij)-obj_self.d0)/(abs(obj_self.d1-obj_self.d0)))*norm(obj_self.d0^2-obj_self.d^2)/(alpha+obj_self.d0^2)^2*(obj_self.Position-obj(j).Position);
                end
            end
        end
        
        function obj=GenerateControl(obj)
            N=length(obj);
            k_di=1*eye(2);
            delta_t=obj(1).delta_t;
            for i=1:N
                DiffPotential=1*obj.CalcDiffPotentialU(obj(i));
                u_i=-1*DiffPotential'-k_di*(delta_t/obj(i).gamma)*obj(i).v;
                for j=obj(i).Neighbour
                    u_i=u_i+obj(i).k3*(delta_t*(obj(j).v-obj(i).v)/obj(i).gamma);
                end
                u_i=u_i+obj(i).k4*(vec__phi(obj(i).Position,obj(i).Kernels,obj(i).Gamma,obj(i).sigma)./obj(i).sigma.^2.*(-obj(i).Position+obj(i).Kernels))'*obj(i).Theta_est;
                gamma_old=obj(i).gamma;
                obj(i).gamma=obj(i).gamma*1/1.0001;
                obj(i).v=obj(i).gamma/delta_t*(delta_t/gamma_old*obj(i).v+gamma_old*u_i);
            end
        end
        
        function obj=Move(obj)
            N=length(obj);
            for i=1:N
                obj(i).Position=obj(i).Position+obj(i).delta_t*obj(i).v';
            end
        end
        
    end
end

