% PSO Code

clear all
clc
close all 
%% Test Functions

% You can validate the optmizer with test functions.

% Easom Function: -cos(x)*cos(y)*exp(-((x-pi)^2 + (y-pi)^2))
% Eggholder Function: -(y+47)*sin(sqrt(abs((x/2)+(y+47))) -
% x*sin(sqrt(abs(x-(y+47))))
% McCormick Function: sin(x+y)+(x-y)^2 - 1.5*x + 2.5*y + 1
% Holder Table Function: -abs(sin(x)*cos(y)*exp(abs(1-(sqrt(x^2 + y^2))/pi)) 



%% Design Variables Bounds

% Based off of your design variables and the bounds you have, changes the
% bounds below.
% Also if you have more design varibles then add based on the syntax below
% Only change line 19 and 20
x_B= [-10 10]
y_B=[-10 10]

x_max_p = x_B(2);
x_min_p = x_B(1);
y_max_p = y_B(2);
y_min_p = y_B(1);



%% Problem Definiton

% Change nVar based on the number of variables you have
nVar = 2;                                                                   % Number of Unknown (Design) Variables

%% Parameters of PSO

% Determine the number of iterations you want the population to have.
MaxIt = 100;                                                                % Maximum Number of Iterations
% Determine the size of the population you want.
nPop = 500;                                                                 % Population Size (Swarm Size)
VarSize = [1 nVar];                                                         % Matrix Size of Design Variables
% Do not touch the next three lines
alpha = .7;
A = 2;
B = 2;

%% LHS

% Based on the number of design variables you have, modify the line below
data = lhsdesign_modified(nPop,[x_B(1) y_B(1)],[x_B(2) y_B(2)])             % Spreads out the population based off of LHS (Latin HyperSquare)

%% velocity

% Now we have to determine the speed of each particle.
% You can mess around with the four lines below but I find that 2 and 1
% usually works as the max and min.
% Add on the velocity based off of the number of design variables you have.
x_max_v = 2;
x_min_v = 1;
y_max_v = 2;
y_min_v = 1;
vmax_x = (x_max_v - x_min_v)*.2;
vmin_x = -vmax_x;
vmax_y = .2*(y_max_v - y_min_v);
vmin_y = -vmax_y;



%% Initialization

empty_particle.Position = [];
empty_particle.Velocity = [];
empty_particle.Cost = [];
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];

% Create Population Array
particle = repmat(empty_particle, nPop, 1);

% Initialize Global Best
GlobalBest.Cost = 100000;

% Initialize Population Members
for i=1:nPop
    
    % Modifiy lines 88-96 based on the number of desing variables you have.
    % Change line 97 to your objective function or whatever you are
    % trying to minimize.
    
    % Generate Random Solution
    particle(i).Position(1) = data(i,1);
    particle(i).Position(2) = data(i,2);
    
    % Initialize Velocity
    particle(i).Velocity = zeros(VarSize);
    
    % Evaluation
    x = particle(i).Position(1);                                          
    y = particle(i).Position(2);
    
    % Objective Function
    particle(i).Cost = -cos(x)*cos(y)*exp(-((x-pi)^2 + (y-pi)^2));
    
    % Update the Personal Best
    particle(i).Best.Position = particle(i).Position;
    particle(i).Best.Cost = particle(i).Cost;
    % Update Global Best
    if particle(i).Best.Cost < GlobalBest.Cost
        GlobalBest = particle(i).Best;
    end
end

% Array to Hold Best Cost Value on Each Iteration
BestCosts = zeros(MaxIt, 1);
%% Main Loop of PSO

% Notice that the loop below is almost the same as the one you modified
% above. 
% Make the same changes as above but also there are other modifications
% that need to be added on to it.
% Follow the comments to make the proper changes

for it=1:MaxIt
    for i=1:nPop
        
        % Update Velocity
        % Modify the lines below based on the number of design variables
        % you have.
        
        x_rand = x_B(1) + (x_B(2)-x_B(1))*rand(1);
        y_rand = y_B(1) + (y_B(2)-y_B(1))*rand(1);
        particle(i).Velocity = alpha*particle(i).Velocity + A*[x_rand y_rand ].*(particle(i).Best.Position - particle(i).Position) + B*[x_rand y_rand].*(GlobalBest.Position - particle(i).Position);
        
        % Apply velocity limits
        % Add the velocity limits for each design variable.
        
        if particle(i).Velocity(1) > vmax_x
            particle(i).Velocity(1) = vmax_x;
        elseif particle(i).Velocity(1) < vmin_x
            particle(i).Velocity(1) = vmin_x;
        end
        if particle(i).Velocity(2) > vmax_y
            particle(i).Velocity(2) = vmax_y;
        elseif particle(i).Velocity(2) < vmin_y
            particle(i).Velocity(2) = vmin_y;
        end
        
        % Update Position
        particle(i).Position = particle(i).Position + particle(i).Velocity;
        
        % Apply Position Limits
        % This part checks to makes sure no design varible is outside the
        % bounds that were given.
        % Modify based on the number of design varibles you have.
        
        if particle(i).Position(1) > x_max_p
            particle(i).Position(1) = x_max_p;
        elseif particle(i).Position(1) < x_min_p
            particle(i).Position(1) = x_min_p;
        end
        if particle(i).Position(2) > y_max_p
            particle(i).Position(2) = y_max_p;
        elseif particle(i).Position(2) < y_min_p
            particle(i).Position(2) = y_min_p;
        end
        
        
        % Evaluation
        % Set each part of the particle equal to their respective design
        % variables and evaluate the objective function.
        
        x = particle(i).Position(1);
        y = particle(i).Position(2); 
        particle(i).Cost = -cos(x)*cos(y)*exp(-((x-pi)^2 + (y-pi)^2));
        
        % Update Personal Best
        if particle(i).Cost < particle(i).Best.Cost
            particle(i).Best.Position = particle(i).Position;
            particle(i).Best.Cost = particle(i).Cost;
            
            % Update Global Best
            if particle(i).Best.Cost < GlobalBest.Cost
                GlobalBest = particle(i).Best;
            end
        end
    end
    
    % Store the Best Cost Value
    BestCosts(it) = GlobalBest.Cost;
    
    % The best value is stored in GobalBest.Cost and GlobalBest.Position
    
    % Display Iteration Information
    % This part to just create a visual to show the convergence of the
    % particles.
    % It is pretty basic so I will not comment on how it works.
    
    figure(1)
    parts = zeros(3,nPop);
    for part = 1:nPop
        disp(particle(part))
        parts(1,part) = particle(part).Position(1);
        parts(2,part) = particle(part).Position(2);
        parts(3,part) = particle(part).Cost;
    end
    clf
    plot(parts(1,:),parts(2,:),'bx')
    xlim([-5 5])
    ylim([-5 5])
    zlim([-1 0])
    hold on
    figure (2)
    plot3(parts(1,:),parts(2,:),parts(3,:),'bx')
    xlim([-5 5])
    ylim([-5 5])
    zlim([-1 0])
    pause(.1)
    disp(['Iteration' num2str(it) ':Best Cost =' num2str(BestCosts(it))]);
    
end