%% simulator.m
%%% APRIL 10, 2021

classdef simulator
    
    properties (SetAccess = protected)
        t;
        x;
    end
    
    properties (SetAccess = public)
        opt = odeset('RelTol', 1e-08, 'MaxStep', 1e-01, 'Refine', 10);
        tf;
        x0;
        N;
        ts;
    end
    
    methods (Access = public)
        
        function obj = simulator % Constructor
            obj.tf = 2000;
            obj.N = 20;
            obj.ts = 1;
            obj = obj.set_init_con;
        end
        
        function obj = simulate_model(obj)
            tspan = 0:(obj.ts):(obj.tf);
            [obj.t, obj.x] = ode23t(@obj.model_odes, tspan, ...
                obj.x0, obj.opt);
        end
        
        function plot_simulation(obj)
            Theta = mod(obj.x(:,1:2:2*obj.N), 2*pi);
            figure('Position', [1 1 720 720], 'DefaultAxesFontSize', ...
                24, 'DefaultAxesLineWidth', 3,'Renderer', 'Painters');
            subplot(2,1,1)
            plot(obj.t, Theta, 'LineWidth', 2)
            set(gca, 'XTickLabel', {''})
            set(gca, 'YLim', [0 6.5], 'YTick', [0 pi 2*pi], ...
                'YTickLabel', {'0' '\pi' '2\pi'})
            ylabel('$\vartheta_i$', 'Interpreter', 'latex')
            subplot(2,1,2)
            plot(obj.t, obj.x(:,2:2:2*obj.N), 'LineWidth', 2)
            xlabel('Time (min)')
            ylabel('$I_i$', 'Interpreter', 'latex')
        end
        
    end
    
    methods (Access = protected)
        
        function dxdt = model_odes(obj,~, x)
            %% Set model parameters
            omega = 2 * pi / 80;
            theta_c = .25 * 2 * pi;
            delta = .01;
            alpha = 1;
            eta = 2;
            Q = .8;
            nu = .5;
            h = 4; 
            kappa = alpha ./ (delta + eta) .* (delta + eta .*(1 + Q .* ...
                (nu - 1))) ./ (delta + eta .* (1 - Q));
            
            %% Preallocate array dxdt
            dxdt = nan(2*obj.N,1);
            
            %% Retrieve state
            Theta(:,1) = mod(x(1:2:2*obj.N), 2*pi);
            I(:,1) = x(2:2:2*obj.N);
            E = Q * mean(I);
            
            %% Compute switching function
            sigma = Theta < theta_c;
            u = (I .^ h) ./ (I .^ h + kappa .^ h);
            
            %% ODEs
            dxdt(1:2:(2*obj.N),1) = omega .* (1 + sigma .* (u - 1));
            dxdt(2:2:(2*obj.N),1) = - delta .* I + alpha .* .5 .* ...
                (sin(Theta)+1) + eta .* (E - I);
        end
        
        function obj = set_init_con(obj)
            %% Set initial conditions of oscillators' phase.
            Theta_0 = linspace(0, 2*pi, obj.N+1);
            %% Set state initial conditions
            obj.x0 = zeros(2*obj.N,1);
            obj.x0(1:2:2*obj.N,:) = Theta_0(1:end-1);
        end
        
    end
    
end