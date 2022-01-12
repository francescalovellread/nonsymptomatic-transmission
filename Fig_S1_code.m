%% Interventions targeting nonsymptomatic cases can be important to prevent local outbreaks: SARS-CoV-2 as a case-study
% Code for generating Supplementary Figure S1
% Francesca Lovell-Read (francesca.lovell-read@merton.ox.ac.uk)
% Version of: Wednesday 21st April 2021
clear

%% SPECIFY PARAMETERS

P = 1000;       % Total population size                     
I0 = 1;         % Initial number infected
S0 = P-I0;      % Initial number susceptible
R0 = 0;         % Initial number removed

beta = 1/2500;  % Infection rate
mu = 1/5;       % Removal rate
tFinal = 100;   % Final time for simulations
numSims = 100000; % Number of stochastic runs

%% RUN STOCHASTIC SIMULATIONS: COMPARTMENTAL

% Create cell array for storing results
simData = cell(numSims,2);

for i=1:numSims
    % Initialize
    t = 0; I = I0; S = S0; R = R0;
    tvec = []; Ivec = []; Svec = []; Rvec  = [];
    index = 1; tvec(index)=t; Ivec(index)=I0; Svec(index)=S0; Rvec(index) = R0;
    infCount = 0;
    while (t<tFinal)
        a1 = beta*S*I; a2 = mu*I; a0 = a1+a2; % Compute reaction propensities
        r1 = rand(1); tau = (1/a0)*log(1/r1); % Compute time to next event
        r2 = rand(1); % Decide which even happens next
        if r2<a1/a0
            I = I+1; S = S-1; infCount = infCount+1;
        else
            I = I-1; R = R+1;
        end
        t = t+tau; 
        index = index + 1; tvec(index) = t; Ivec(index) = I; Svec(index) = S;
    end
simData{i,1}=tvec; simData{i,2}=Ivec; % Store simulated t and I vectors
Imax(i) = max(Ivec); % Compute and store maximum number simultaneously infected
everInf(i) = infCount; % Store maximum number ever infected
end

%% PLOTTING: COMPARTMENTAL

figure(); hold on; box on; set(gca,'Fontsize',22); grid on; grid minor;

% Sort infection counts into histogram
edges = 0:P/50:P;
midpoints = P/100:P/50:P-P/100;
h = histcounts(everInf,edges); h = h/sum(h);
hbar = bar(midpoints,h);
hbar.FaceColor = [0.5 0.8 1];

xlim([0 P]); ylim([0 0.55]);
xlabel('Number ever infected')
ylabel('Probability')

%% RUN STOCHASTIC SIMULATIONS: BRANCHING PROCESS

% Create cell array for storing results
simData = cell(numSims,2);

for i=1:numSims
    if mod(i,1000)==0
        i
    end
    % Initialize
    t = 0; I = I0; S = S0; R = R0;
    tvec = []; Ivec = []; Svec = []; Rvec  = [];
    index = 1; tvec(index)=t; Ivec(index)=I0; Svec(index)=S0; Rvec(index) = R0;
    infCount = 0;
    while (t<tFinal && I<=1000)
        a1 = beta*S*I; a2 = mu*I; a0 = a1+a2; % Compute reaction propensities
        r1 = rand(1); tau = (1/a0)*log(1/r1); % Compute time to next event
        r2 = rand(1); % Decide which even happens next
        if r2<a1/a0
            I = I+1; infCount = infCount+1;
        else
            I = I-1; R = R+1;
        end
        t = t+tau; 
        index = index + 1; tvec(index) = t; Ivec(index) = I; Svec(index) = S;
    end
simData{i,1}=tvec; simData{i,2}=Ivec; % Store simulated t and I vectors
Imax(i) = max(Ivec); % Compute and store maximum number simultaneously infected
everInf(i) = infCount; % Store maximum number ever infected
end

%% PLOTTING: BRANCHING PROCESS

figure(); hold on; box on; set(gca,'Fontsize',22); grid on; grid minor;

% Sort infection counts into histogram
edges = [0:P/50:P-P/50 max(everInf)];
midpoints = [P/100:P/50:P-P/100];
h = histcounts(everInf,edges); h = h/sum(h);
hbar = bar(midpoints,h);
hbar.FaceColor = [0.5 0.8 1];

xlim([0 P]); ylim([0 0.55]);
xticklabels = ({'0'; '200'; '400'; '600'; '800'; '1000+'});
xticks = linspace(0,1000, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);
xlabel('Number ever infected')
ylabel('Probability')
