%% Interventions targeting nonsymptomatic cases can be important to prevent local outbreaks: SARS-CoV-2 as a case-study
% Code for generating Fig 2A
% Francesca Lovell-Read (francesca.lovell-read@merton.ox.ac.uk)
% Version of: Wednesday 21st April 2021
clear

%% SET EPIDEMIOLOGICAL PARAMETERS

% BASELINE PARAMETERS
R0_bl = 3;          % R0 (basic reproduction number)
xi = 0.2;           % Proportion of infections which are asymptomatic
gamma = 0.0924;     % Isolation rate of symptomatic hosts without intensified surveillance
epsilon = 0.1;      % Relative isolation rate of nonsymptomatic hosts
lambda = 1/2;       % Rate at which presymptomatic hosts develop symptoms
mu = 1/8;           % Recovery rate of symptomatic hosts
nu = 0.1;           % Recovery rate of asymptomatic hosts
Kp = 0.4894;        % Proportion of transmissions coming from presymptomatic hosts
Ka = 0.1064;        % Proportion of transmissions coming from asymptomatic hosts

% Calculate alpha and eta (chosen to give specified proportions of presymptomatic and asymptomatic transmissions)
alpha = (lambda/(gamma+mu))*(Kp/(1-Ka-Kp));
eta = ((1-xi)/xi)*((nu+epsilon*gamma)/(lambda+epsilon*gamma))*(Ka/Kp)*alpha;

% Set vector of R0 values to consider
R0_vec = 0.1:0.1:5; 
% Set vector of lambda values to consider
lambda_vec = [1 1/2 1/4];
% Compute corresponding nu values that maintain a constant ratio of
% presymptomatic:asymptomatic transmissions
C = Ka/Kp;
nu_vec = (1/C)*(xi/(1-xi))*(eta/alpha)*(lambda_vec+epsilon*gamma)-epsilon*gamma;

%% COMPUTE OUTBREAK PROBABILITIES

% Create empty matrices for storing values
M_asymp = zeros(length(lambda_vec),length(R0_vec));
M_presymp = zeros(length(lambda_vec),length(R0_vec));
M_symp = zeros(length(lambda_vec),length(R0_vec));

for j = 1:length(lambda_vec) % Loop over lambda/nu values
    lambda = lambda_vec(j);
    nu = nu_vec(j);
    
    % The corresponding new values of Kp, Ka can be calculated with the two
    % lines below, if desired
    Kp = ((1-xi)*(alpha/(lambda+epsilon*gamma)))/(xi*eta/(nu+epsilon*gamma)+(1-xi)*(alpha/(lambda+epsilon*gamma)+(lambda/(lambda+epsilon*gamma))*(1/(gamma+mu))));
    Ka = (xi*eta/(nu+epsilon*gamma))/(xi*eta/(nu+epsilon*gamma)+(1-xi)*(alpha/(lambda+epsilon*gamma)+(lambda/(lambda+epsilon*gamma))*(1/(gamma+mu))));

    for i = 1:length(R0_vec) % Loop over R0 values
        R0 = R0_vec(i);

        % Calculate beta (chosen to give specified R0)
        g2 = gamma;         % Baseline isolation rate for symptomatic hosts
        g1 = epsilon*gamma; % Baseline isolation rate for presymptomatic hosts (scaled by epsilon)
        g2_inv = 1/g2;
        g1_sum = g1 + lambda;
        g1_sum_inv = 1/g1_sum;
        g = (1-xi)*g1_sum_inv*(alpha+lambda/(mu+g2)) + xi*eta/(nu+g1);

        beta = R0/g;
        beta2 = beta;       % Transmission coefficient from symptomatic hosts
        beta1 = alpha*beta; % Transmission coefficient from presymptomatic hosts (scaled by alpha)
        beta3 = eta*beta;   % Transmission coefficient from asymptomatic hosts (scaled by eta)
    
        % Determine minimal non-negative solution to eqns 3,4,5:

        % Define a,b,c,d   
        a = (beta1)./(g1+beta1+lambda);
        b = (lambda)./(g1+beta1+lambda);
        c = (beta2)./(g2+beta2+mu);
        d = (beta3)./(g1+beta3+nu);

        % Define polynomial coefficients    
        Z4 = d*(a-d)*xi*(d-c);
        Z3 = d*(a-d)*xi*c*(1-d) + (d-a-a*d*xi+(d^2)*(a-1+b*(1-xi)+xi))*(d-c) - b*(1-c)*(d^3)*(1-xi);
        Z2 = (d-a-a*d*xi+(d^2)*(a-1+b*(1-xi)+xi))*(c*(1-d)) + (d*(d-2*a-1)+2*a)*(d-c);
        Z1 = (d*(d-2*a-1)+2*a)*c*(1-d) -a*(1-d)^2*((d-c));
        Z0 = -a*(1-d)^2*(c*(1-d));

        % Define polynomial to solve and determine smallest non-negative root
        poly = [Z4 Z3 Z2 Z1 Z0];
        poly_roots = roots(poly);
        poly_roots = real(poly_roots(abs(imag(poly_roots))<1e-5));
        poly_roots = poly_roots(poly_roots>=0);
        q_minus = min(poly_roots);

        % If this root is in [0,1] then take it as solution; otherwise take 1. 
        if(q_minus>=0 && q_minus<1)
            z = q_minus;
            x = (1/(d*(1-xi)))*(1-d*xi*q_minus-(1-d)/q_minus);
            y = (1-c)*d*z/(d*z-c*(z-1+d));
            PMO_asymp = 1-z;
            PMO_presymp = 1-x;
            PMO_symp = 1-y;
        else 
            PMO_asymp = 0;
            PMO_presymp = 0;
            PMO_symp = 0;
        end
        
        M_asymp(j,i) = PMO_asymp;
        M_presymp(j,i) = PMO_presymp;
        M_symp(j,i) = PMO_symp;
    end

end

% Compute matrix of local outbreak probabilities starting from a single nonsymptomatic
% host (weighted average of asymptomatic and presymptomatic matrices)
M_nonsymp = xi*M_asymp+(1-xi)*M_presymp;

%% PLOTTING

figure(1); hold on; box off;

color1 = [0.5 0.3 1]; color2 = [0.3 0.6 1]; color3 = [0.4 0.7 0];
colours = [color1; color2; color3];

% For each lambda, plot computed outbreak probabilities
for j = 1:length(lambda_vec)
    myplot(j+1) = plot(R0_vec,M_nonsymp(j,:),'linewidth',2,'color',colours(j,:));
end

% Plot '1-1/R0' 
R0_vec_geq_1 = R0_vec(R0_vec>=1);
myplot(1) = plot(R0_vec_geq_1,1-1./R0_vec_geq_1,'r-.','linewidth',2);
line([0 1],[0 0],'color','r','linestyle','-.','linewidth',2)

% Plot baseline R0
line([R0_bl R0_bl],[0 1],'color',[0.7 0.7 0.7],'linestyle',':','linewidth',2)

% Legend
Legend=cell(length(lambda_vec),1);
Legend{1}='1-1/R_0';
for i=1:length(lambda_vec)
if lambda_vec(i)==1;
    Legend{i+1}=strcat('1/\lambda =',32,num2str(1/lambda_vec(i)),32,'day');
else
    Legend{i+1}=strcat('1/\lambda =',32,num2str(1/lambda_vec(i)),32,'days');
end
end
leg = legend(myplot(1:1+length(lambda_vec)),Legend);
leg.Position = [0.75 0.285 0 0]; leg.Box = 'on'; leg.FontSize = 16;

xlim([0 R0_vec(end)]); ylim([0 1]);
xlabel('Basic reproduction number (R_0)');
ylabel('Probability of a local outbreak (p)');

set(gca,'Fontsize', 16,'linewidth',1); 

