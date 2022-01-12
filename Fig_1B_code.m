%% Interventions targeting nonsymptomatic cases can be important to prevent local outbreaks: COVID-19 as a case-study
% Code for generating Fig 1B
% Francesca Lovell-Read (francesca.lovell-read@merton.ox.ac.uk)
% Version of: Friday 5th March 2021
clear

%% 

% Define vector of rho values
rho_vec = 0:0.1:20;
% Define vector of delta values
delta_vec = [0.6 0.7 0.8 0.9];
% Compute corresponding f values
f_of_rho = zeros(length(delta_vec),length(rho_vec));
for i=1:length(delta_vec)
    f_of_rho(i,:) = delta_vec(i).*rho_vec./(1+rho_vec);
end

%% PLOTTING
figure(1); hold on; box off; set(gca,'fontsize',16,'linewidth',2);

color1 = [0.8 0.1 0.3]; color2 = [0.5 0.3 1]; color3 = [0.3 0.6 1]; color4 = [0.4 0.7 0];
colours = [color1; color2; color3; color4;];

% Plot graphs of f(rho) and asymptotes
for i = 1:length(delta_vec)
    myplot(i) = plot(rho_vec,f_of_rho(i,:),'linewidth',2,'color',colours(i,:));
    myplot(i+length(delta_vec)) = yline(delta_vec(i),'linewidth',2,'linestyle','--','color',colours(i,:));
end

% Legend
Legend=cell(length(delta_vec),1);
for i=1:length(delta_vec)
    Legend{i}=strcat('\delta =',32,num2str(delta_vec(i)));
end
leg = legend(myplot(1:length(delta_vec)),Legend);
leg.Position = [0.75 0.3 0 0]; leg.Box = 'on'; leg.FontSize = 16;

xlabel('Surveillance intensification effort (\rho)');
ylabel({'Reduction in expected time to isolation (f(\rho,\delta))'});
xlim([0 rho_vec(end)]); ylim([0 1]);
yticks = [0 0.2 0.4 0.6 0.8 1];
set(gca,'YTick',yticks);

