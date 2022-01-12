%% Interventions targeting nonsymptomatic cases can be important to prevent local outbreaks: SARS-CoV-2 as a case-study
% Code for generating Figs 2B,C,D, Fig 3, Figs S2-S11
% Francesca Lovell-Read (francesca.lovell-read@merton.ox.ac.uk)
% Version of: Wednesday 21st April 2021
clear

%% SET EPIDEMIOLOGICAL PARAMETERS

R0 = 3; % Basic reproduction number
xi = 0.2; % Proportion of infections which are asymptomatic
gamma = 0.0924; % Isolation rate of symptomatic hosts without intensified surveillance
epsilon = 0.1; % Relative isolation rate of nonsymptomatic hosts
lambda = 1/2; % Rate at which presymptomatic hosts develop symptoms
mu = 1/8; % Recovery rate of symptomatic hosts
nu = 0.1; % Recovery rate of asymptomatic hosts
delta = 0.8; % Upper bound on fractional reduction in time to isolation
rho_max = 20; % Upper bound on surveillance intensification effort (rho)
Kp = 0.4894; % Proportion of transmissions coming from presymptomatic hosts
Ka = 0.1064; % Proportion of transmissions coming from asymptomatic hosts

%% SET PROGRAM PARAMETERS

M_size = 100; % Model resolution (size of outbreak probability matrix)

% Define search angles and radius for computing smooth and 'two-way' contours
theta_vec_1 = 0:pi/64:pi/2; % Search angles for smooth radial search contours
theta_vec_2 = 0:pi/2:pi/2;  % Search angles for 'two-way' contours
search_radius = 0.01;       % Search radius/step size (0.01 good step size for smoothness)

% Define search points for computing linear search contours
rho1_step = 0:0.1*search_radius:search_radius;
rho2_step = search_radius-rho1_step;

buffer = 20000; % A padding buffer to ensure vectors have consistent dimensions for concatenation
% NB If you see the error 'Dimensions of arrays being concatenated are not
% consistent', try increasing the buffer size.

total_effort_vec = [2 6 10 14 18 22 26 30 34 38]; % Vector of 'total effort' values considered.
% This determines which 'fixed effort' contours are plotted in Fig 3C, 
clims = [0.555 0.731]; % Limits for plotting heat map

%% PRELIMINARY CALCULATIONS

% Calculate alpha and eta (chosen to give specified proportions of presymptomatic and asymptomatic transmissions)
alpha = (lambda/(gamma+mu))*(Kp/(1-Ka-Kp));
eta = ((1-xi)/xi)*((nu+epsilon*gamma)/(lambda+epsilon*gamma))*(Ka/Kp)*alpha;

% Calculate beta (chosen to give specified R0)
g2 = gamma;         % Baseline isolation rate for symptomatic hosts
g1 = epsilon*gamma; % Baseline isolation rate for nonsymptomatic hosts (scaled by epsilon)
g2_inv = 1/g2;
g1_sum = g1 + lambda;
g1_sum_inv = 1/g1_sum;
g = (1-xi)*g1_sum_inv*(alpha+lambda/(mu+g2)) + xi*eta/(nu+g1);

beta = R0/g;
beta2 = beta;       % Transmission coefficient from symptomatic hosts
beta1 = alpha*beta; % Transmission coefficient from presymptomatic hosts (scaled by alpha)
beta3 = eta*beta;   % Transmission coefficient from asymptomatic hosts (scaled by eta)

%% COMPUTE LOCAL OUTBREAK PROBABILITY MATRIX
% This is a matrix containing the local outbreak probability values
% over the whole range of (rho_1,rho_2) values
fprintf('Computing probability matrix...');

M_symp = zeros(M_size+1);     % Empty matrix for storing PLO values starting from symptomatic individual
M_presymp = zeros(M_size+1);  % Empty matrix for storing PLO values starting from presymptomatic individual
M_asymp = zeros(M_size+1);    % Empty matrix for storing PLO values starting from asymptomatic individual

for j=1:M_size+1
    rho1 = (rho_max/M_size)*(j-1);  % loop over rho_1 (nonsymptomatic surveillance)
for k=1:M_size+1
    rho2 = (rho_max/M_size)*(k-1);  % loop over rho_2 (symptomatic surveillance)
    
    % For each (rho_1,rho_2) determine minimal non-negative solution to eqns 3,4,5:
    
    % Define a,b,c,d (event probabilities)
    a = ((1-(delta*rho1)./(1+rho1))*beta1)./(g1+(1-(delta*rho1)./(1+rho1))*(beta1+lambda));
    b = ((1-(delta*rho1)./(1+rho1))*lambda)./(g1+(1-(delta*rho1)./(1+rho1))*(beta1+lambda));
    c = ((1-(delta*rho2)./(1+rho2))*beta2)./(g2+(1-(delta*rho2)./(1+rho2))*(beta2+mu));
    d = ((1-(delta*rho1)./(1+rho1))*beta3)./(g1+(1-(delta*rho1)./(1+rho1))*(beta3+nu));

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
        PLO_asymp = 1-z;
        PLO_presymp = 1-x;
        PLO_symp = 1-y;
    else 
        PLO_asymp = 0;
        PLO_presymp = 0;
        PLO_symp = 0;
    end
    
    % Fill PLO matrices!
    M_asymp(j,k) = PLO_asymp;
    M_presymp(j,k) = PLO_presymp;
    M_symp(j,k) = PLO_symp;

end
end
fprintf(' Done!\n');

% Compute matrix of local outbreak probabilities starting from a single nonsymptomatic
% host (weighted average of asymptomatic and presymptomatic matrices)
M_nonsymp = xi*M_asymp+(1-xi)*M_presymp;

%% CREATE MATRIX CONTAINING INFORMATION ON WHICH INTERVENTION HAS GREATER EFFECT ('GRADS MATRIX')
fprintf('Computing grads matrix...');

grads = zeros(M_size,M_size); % Matrix of zeros to fill up with values

% Scan over nonsymptomatic PLO matrix (M_nonsymp) and at every point
% determine whether increasing surveillance for symptomatic hosts or for
% nonsymptomatic hosts makes more difference

for j=1:M_size
    for k=1:M_size

        PLO_value = M_nonsymp(j,k);         % PLO at given point
        PLO_value_rho1 = M_nonsymp(j+1,k);  % New PLO when nonsymptomatic surveillance is increased
        PLO_value_rho2 = M_nonsymp(j,k+1);  % New PLO when symptomatic surveillance is increased

        del_rho1 = PLO_value-PLO_value_rho1;  % Difference in PLO from increasing nonsymptomatic surveillance (rho1)
        del_rho2 = PLO_value-PLO_value_rho2;  % Difference in PLO from increasing symptomatic surveillance (rho2)

        if abs(del_rho1-del_rho2)<0.000001  % Negligible difference between effect of changing rho1 and rho2
            grads(j,k)=1; 
        elseif del_rho1>del_rho2            % Nonsymptomatic surveillance has more effect
            grads(j,k)=2;
        end
    end
end
fprintf(' Done!\n');

% Grads matrix is 0 when increasing symptomatic surveillance has a greater
% effect, 1 when there is a negligible differene between increasing
% symptomatic surveillance and nonsymptomatic surveillance and 2 when
% increasing nonsymptomtatic surveillance has a greater effect. 

%% COMPUTE STEEPEST DESCENT CONTOURS USING THETA_VEC_1: 'SMOOTH' RADIAL SEARCH CONTOURS
fprintf('Computing smooth radial search contours...');

% Compute contours starting from rho2=0 (i.e. along rho1 axis)
rho1_radial_search_1 = [];  
rho2_radial_search_1 = []; % Empty vectors ready to fill with values

for j=0:4:1*rho_max-1 % j indexes the contours: consider each contour one at a time
    
    rho1 = 1*j;
    rho2 = 0;
    rho1_plot = [rho1];
    rho2_plot = [rho2];

    while (rho1<rho_max+1)
        % Create vectors of rho1 and rho2 values achievable using given search angles and search radius
        rho1 = rho1 + search_radius*cos(theta_vec_1);
        rho2 = rho2 + search_radius*sin(theta_vec_1);

        % For each search angle compute PLO for associated rho1 and rho2:

        PLO_asymp = zeros(1,length(theta_vec_1));
        PLO_presymp = zeros(1,length(theta_vec_1));
        PLO_symp = zeros(1,length(theta_vec_1));

        % Define a,b,c,d   
        a = ((1-(delta*rho1)./(1+rho1))*beta1)./(g1+(1-(delta*rho1)./(1+rho1))*(beta1+lambda));
        b = ((1-(delta*rho1)./(1+rho1))*lambda)./(g1+(1-(delta*rho1)./(1+rho1))*(beta1+lambda));
        c = ((1-(delta*rho2)./(1+rho2))*beta2)./(g2+(1-(delta*rho2)./(1+rho2))*(beta2+mu));
        d = ((1-(delta*rho1)./(1+rho1))*beta3)./(g1+(1-(delta*rho1)./(1+rho1))*(beta3+nu));

        % Define polynomial coefficients    
        Z4 = d.*(a-d).*xi.*(d-c);
        Z3 = d.*(a-d).*xi.*c.*(1-d) + (d-a-a.*d.*xi+(d.^2).*(a-1+b.*(1-xi)+xi)).*(d-c) - b.*(1-c).*(d.^3).*(1-xi);
        Z2 = (d-a-a.*d.*xi+(d.^2).*(a-1+b.*(1-xi)+xi)).*(c.*(1-d)) + (d.*(d-2.*a-1)+2.*a).*(d-c);
        Z1 = (d.*(d-2.*a-1)+2.*a).*c.*(1-d) -a.*(1-d).^2.*((d-c));
        Z0 = -a.*(1-d).^2.*(c.*(1-d));

        for i=1:length(theta_vec_1)
            poly = [Z4(i) Z3(i) Z2(i) Z1(i) Z0(i)];
            poly_roots = roots(poly); % Find roots of polynomial

            % Obtain the minimal non-negative real root:
            poly_roots = real(poly_roots(abs(imag(poly_roots))<1e-5));
            poly_roots = poly_roots(poly_roots>=0);
            q_minus = min(poly_roots);

                if(q_minus>=0 && q_minus<1)
                    z = q_minus;
                    x = (1/(d(i)*(1-xi)))*(1-d(i)*xi*q_minus-(1-d(i))/q_minus);
                    y = (1-c(i))*d(i)*z/(d(i)*z-c(i)*(z-1+d(i)));
                    PLO_asymp(i) = 1-z;
                    PLO_presymp(i) = 1-x;
                    PLO_symp(i) = 1-y;
                else 
                    PLO_asymp(i) = 0;
                    PLO_presymp(i) = 0;
                    PLO_symp(i) = 0;
                end        
        end

        PLO_nonsymp = xi*PLO_asymp+(1-xi)*PLO_presymp; % Vector of PLO values corresponding to each search angle

        % Find the search angle which minimises the PLO value: 
        if min(PLO_nonsymp)>0   % i.e. if contour has not hit zero
            min_id = min(find(PLO_nonsymp==min(PLO_nonsymp)));  % Find the index of the point for which the PLO is minimised
            rho1 = rho1(min_id);    % Take the rho1 corresponding to this index
            rho2 = rho2(min_id);    % Take the rho2 corresponding to this index
            rho1_plot = [rho1_plot rho1];   % Append new rho1 value to vector of previous rho1 values
            rho2_plot = [rho2_plot rho2];   % Append new rho2 value to vector of previous rho2 values
        elseif((min(PLO_nonsymp)==0)&&(max(PLO_nonsymp)>0)) % i.e. if the contour has hit zero
            min_id = floor(mean(find(PLO_nonsymp==0))); % Find indices for which the PLO is zero and go down the middle of this region
            rho1 = rho1(min_id);
            rho2 = rho2(min_id);
            rho1_plot = [rho1_plot rho1];
            rho2_plot = [rho2_plot rho2];
        end
    end

    % Store rho_1 and rho_2 values for this contour
    rho1_plot(buffer) = 0;  % Pad with buffer to make sure dimensions are consistent (since contours will be different lengths)
    rho2_plot(buffer) = 0;  % -"-
    rho1_radial_search_1 = [rho1_radial_search_1;rho1_plot];    % Append new rho1 vector to matrix of vectors for previous contours
    rho2_radial_search_1 = [rho2_radial_search_1;rho2_plot];    % Append new rho2 vector to matrix of vectors for previous contours
end

% Compute contours starting from rho1=0 (i.e. along rho2 axis)
rho1_radial_search_2 = [];
rho2_radial_search_2 = []; % Empty vectors ready to fill with values

for j=0:4:1*rho_max-1 % j indexes the contours: consider each contour one at a time
        
    rho1 = 0;
    rho2 = 1*j;
    rho1_plot = [rho1];
    rho2_plot = [rho2];

    while (rho1<rho_max+1)
        % Create vectors of rho1 and rho2 values achievable using given search angles and search radius
        rho1 = rho1 + search_radius*cos(theta_vec_1);
        rho2 = rho2 + search_radius*sin(theta_vec_1);
        
        % For each search angle compute PLO for associated rho1 and rho2:

        PLO_asymp = zeros(1,length(theta_vec_1));
        PLO_presymp = zeros(1,length(theta_vec_1));
        PLO_symp = zeros(1,length(theta_vec_1));
        
        % Define a,b,c,d   
        a = ((1-(delta*rho1)./(1+rho1))*beta1)./(g1+(1-(delta*rho1)./(1+rho1))*(beta1+lambda));
        b = ((1-(delta*rho1)./(1+rho1))*lambda)./(g1+(1-(delta*rho1)./(1+rho1))*(beta1+lambda));
        c = ((1-(delta*rho2)./(1+rho2))*beta2)./(g2+(1-(delta*rho2)./(1+rho2))*(beta2+mu));
        d = ((1-(delta*rho1)./(1+rho1))*beta3)./(g1+(1-(delta*rho1)./(1+rho1))*(beta3+nu));

        % Define polynomial coefficients    
        Z4 = d.*(a-d).*xi.*(d-c);
        Z3 = d.*(a-d).*xi.*c.*(1-d) + (d-a-a.*d.*xi+(d.^2).*(a-1+b.*(1-xi)+xi)).*(d-c) - b.*(1-c).*(d.^3).*(1-xi);
        Z2 = (d-a-a.*d.*xi+(d.^2).*(a-1+b.*(1-xi)+xi)).*(c.*(1-d)) + (d.*(d-2.*a-1)+2.*a).*(d-c);
        Z1 = (d.*(d-2.*a-1)+2.*a).*c.*(1-d) -a.*(1-d).^2.*((d-c));
        Z0 = -a.*(1-d).^2.*(c.*(1-d));

        for i=1:length(theta_vec_1)
            poly = [Z4(i) Z3(i) Z2(i) Z1(i) Z0(i)];
            poly_roots = roots(poly); % Find roots of polynomial
            
            % Obtain the minimal non-negative real root:
            poly_roots = real(poly_roots(abs(imag(poly_roots))<1e-5));
            poly_roots = poly_roots(poly_roots>=0);
            q_minus = min(poly_roots);

                if(q_minus>=0 && q_minus<1)
                    z = q_minus;
                    x = (1/(d(i)*(1-xi)))*(1-d(i)*xi*q_minus-(1-d(i))/q_minus);
                    y = (1-c(i))*d(i)*z/(d(i)*z-c(i)*(z-1+d(i)));
                    PLO_asymp(i) = 1-z;
                    PLO_presymp(i) = 1-x;
                    PLO_symp(i) = 1-y;
                else 
                    PLO_asymp(i) = 0;
                    PLO_presymp(i) = 0;
                    PLO_symp(i) = 0;
                end           
        end

        PLO_nonsymp = xi*PLO_asymp+(1-xi)*PLO_presymp; % Vector of PLO values corresponding to each search angle

        % Find the search angle which minimises the PLO value: 
        if min(PLO_nonsymp)>0 % i.e. if contour has not hit zero
            min_id = min(find(PLO_nonsymp==min(PLO_nonsymp))); % Find the index of the point for which the PLO is minimised
            rho1 = rho1(min_id); % Take the rho1 corresponding to this index
            rho2 = rho2(min_id); % Take the rho2 corresponding to this index
            rho1_plot = [rho1_plot rho1]; % Append new rho1 value to vector of previous rho1 values
            rho2_plot = [rho2_plot rho2]; % Append new rho2 value to vector of previous rho2 values
         elseif((min(PLO_nonsymp)==0)&&(max(PLO_nonsymp)>0)) % i.e. if the contour has hit zero
            min_id = floor(mean(find(PLO_nonsymp==0))); % Find indices for which the PLO is zero and go down the middle of this region
            rho1 = rho1(min_id);
            rho2 = rho2(min_id);
            rho1_plot = [rho1_plot rho1];
            rho2_plot = [rho2_plot rho2];
        end
    end
    
    % Store rho_1 and rho_2 values for this contour
    rho1_plot(buffer) = 0; % Pad with buffer to make sure dimensions are consistent (since contours will be different lengths)
    rho2_plot(buffer) = 0; % -"-
    rho1_radial_search_2 = [rho1_radial_search_2;rho1_plot]; % Append new rho1 vector to matrix of vectors for previous contours
    rho2_radial_search_2 = [rho2_radial_search_2;rho2_plot]; % Append new rho2 vector to matrix of vectors for previous contours
end

fprintf(' Done!\n');

%% COMPUTE STEEPEST DESCENT CONTOURS USING THETA_VEC_2: TWO WAY CONTOURS
fprintf('Computing two way contours...');

% Compute contours starting from rho2=0 (i.e. along rho1 axis)
rho1_two_way_search_1 = [];
rho2_two_way_search_1 = []; % Empty vectors ready to fill with values

for j=0:4:1*rho_max-1 % j indexes the contours: consider each contour one at a time
    
    rho1 = 1*j;
    rho2 = 0;
    rho1_plot = [rho1];
    rho2_plot = [rho2];

    while (rho1<rho_max+1)
        % Create vectors of rho1 and rho2 values achievable
        rho1 = rho1 + search_radius*[1 0];
        rho2 = rho2 + search_radius*[0 1];
        
        % For each search angle compute PLO for associated rho1 and rho2:

        PLO_asymp = zeros(1,2);
        PLO_presymp = zeros(1,2);
        PLO_symp = zeros(1,2);
        
        % Define a,b,c,d   
        a = ((1-(delta*rho1)./(1+rho1))*beta1)./(g1+(1-(delta*rho1)./(1+rho1))*(beta1+lambda));
        b = ((1-(delta*rho1)./(1+rho1))*lambda)./(g1+(1-(delta*rho1)./(1+rho1))*(beta1+lambda));
        c = ((1-(delta*rho2)./(1+rho2))*beta2)./(g2+(1-(delta*rho2)./(1+rho2))*(beta2+mu));
        d = ((1-(delta*rho1)./(1+rho1))*beta3)./(g1+(1-(delta*rho1)./(1+rho1))*(beta3+nu));

        % Define polynomial coefficients    
        Z4 = d.*(a-d).*xi.*(d-c);
        Z3 = d.*(a-d).*xi.*c.*(1-d) + (d-a-a.*d.*xi+(d.^2).*(a-1+b.*(1-xi)+xi)).*(d-c) - b.*(1-c).*(d.^3).*(1-xi);
        Z2 = (d-a-a.*d.*xi+(d.^2).*(a-1+b.*(1-xi)+xi)).*(c.*(1-d)) + (d.*(d-2.*a-1)+2.*a).*(d-c);
        Z1 = (d.*(d-2.*a-1)+2.*a).*c.*(1-d) -a.*(1-d).^2.*((d-c));
        Z0 = -a.*(1-d).^2.*(c.*(1-d));

        for i=1:2
            poly = [Z4(i) Z3(i) Z2(i) Z1(i) Z0(i)];
            poly_roots = roots(poly); % Find roots of polynomial
            
            % Obtain the minimal non-negative real root:
            poly_roots = real(poly_roots(abs(imag(poly_roots))<1e-5));
            poly_roots = poly_roots(poly_roots>=0);
            q_minus = min(poly_roots);

            if(q_minus>=0 && q_minus<1)
                z = q_minus;
                x = (1/(d(i)*(1-xi)))*(1-d(i)*xi*q_minus-(1-d(i))/q_minus);
                y = (1-c(i))*d(i)*z/(d(i)*z-c(i)*(z-1+d(i)));
                PLO_asymp(i) = 1-z;
                PLO_presymp(i) = 1-x;
                PLO_symp(i) = 1-y;
            else 
                PLO_asymp(i) = 0;
                PLO_presymp(i) = 0;
                PLO_symp(i) = 0;
            end           
        end

        PLO_nonsymp = xi*PLO_asymp+(1-xi)*PLO_presymp; % Vector of PLO values corresponding to each search angle

        % Find the search angle which minimises the PLO value:
        if min(PLO_nonsymp)>0 % i.e. if contour has not hit zero
            min_id = min(find(PLO_nonsymp==min(PLO_nonsymp))); % Find the index of the point for which the PLO is minimised
            rho1 = rho1(min_id); % Take the rho1 corresponding to this index
            rho2 = rho2(min_id); % Take the rho2 corresponding to this index
            rho1_plot = [rho1_plot rho1]; % Append new rho1 value to vector of previous rho1 values
            rho2_plot = [rho2_plot rho2]; % Append new rho2 value to vector of previous rho2 values
        elseif((min(PLO_nonsymp)==0)&&(max(PLO_nonsymp)>0))
            min_id = floor(mean(find(PLO_nonsymp==0))); % Find indices for which the PLO is zero and go down the middle of this region
            rho1 = rho1(min_id); 
            rho2 = rho2(min_id);
            rho1_plot = [rho1_plot rho1];
            rho2_plot = [rho2_plot rho2];
        end
    end
    
    % Store rho_1 and rho_2 values for this contour
    rho1_plot(buffer) = 0; % Pad with buffer to make sure dimensions are consistent (since contours will be different lengths)
    rho2_plot(buffer) = 0; % -"-
    rho1_two_way_search_1 = [rho1_two_way_search_1;rho1_plot]; % Append new rho1 vector to matrix of vectors for previous contours
    rho2_two_way_search_1 = [rho2_two_way_search_1;rho2_plot]; % Append new rho2 vector to matrix of vectors for previous contours
end

% Compute contours starting from rho1=0 (i.e. along rho2 axis)
rho1_two_way_search_2 = [];
rho2_two_way_search_2 = []; % Empty vectors ready to fill with values

for j=0:4:1*rho_max-1 % j indexes the contours: consider each contour one at a time

    rho1 = 0;
    rho2 = 1*j;
    rho1_plot = [rho1];
    rho2_plot = [rho2];

    while (rho1<rho_max+1)
        % Create vectors of rho1 and rho2 values achievable
        rho1 = rho1 + search_radius*[1 0];
        rho2 = rho2 + search_radius*[0 1];
        
        % For each search angle compute PLO for associated rho1 and rho2:
        
        PLO_asymp = zeros(1,2);
        PLO_presymp = zeros(1,2);
        PLO_symp = zeros(1,2);

        % Define a,b,c,d   
        a = ((1-(delta*rho1)./(1+rho1))*beta1)./(g1+(1-(delta*rho1)./(1+rho1))*(beta1+lambda));
        b = ((1-(delta*rho1)./(1+rho1))*lambda)./(g1+(1-(delta*rho1)./(1+rho1))*(beta1+lambda));
        c = ((1-(delta*rho2)./(1+rho2))*beta2)./(g2+(1-(delta*rho2)./(1+rho2))*(beta2+mu));
        d = ((1-(delta*rho1)./(1+rho1))*beta3)./(g1+(1-(delta*rho1)./(1+rho1))*(beta3+nu));

        % Define polynomial coefficients    
        Z4 = d.*(a-d).*xi.*(d-c);
        Z3 = d.*(a-d).*xi.*c.*(1-d) + (d-a-a.*d.*xi+(d.^2).*(a-1+b.*(1-xi)+xi)).*(d-c) - b.*(1-c).*(d.^3).*(1-xi);
        Z2 = (d-a-a.*d.*xi+(d.^2).*(a-1+b.*(1-xi)+xi)).*(c.*(1-d)) + (d.*(d-2.*a-1)+2.*a).*(d-c);
        Z1 = (d.*(d-2.*a-1)+2.*a).*c.*(1-d) -a.*(1-d).^2.*((d-c));
        Z0 = -a.*(1-d).^2.*(c.*(1-d));

        for i=1:2
            poly = [Z4(i) Z3(i) Z2(i) Z1(i) Z0(i)];
            poly_roots = roots(poly); % Find roots of polynomial
            
            % Obtain the minimal non-negative real root:
            poly_roots = real(poly_roots(abs(imag(poly_roots))<1e-5));
            poly_roots = poly_roots(poly_roots>=0);
            q_minus = min(poly_roots);

            if(q_minus>=0 && q_minus<1)
                z = q_minus;
                x = (1/(d(i)*(1-xi)))*(1-d(i)*xi*q_minus-(1-d(i))/q_minus);
                y = (1-c(i))*d(i)*z/(d(i)*z-c(i)*(z-1+d(i)));
                PLO_asymp(i) = 1-z;
                PLO_presymp(i) = 1-x;
                PLO_symp(i) = 1-y;
            else 
                PLO_asymp(i) = 0;
                PLO_presymp(i) = 0;
                PLO_symp(i) = 0;
            end           
        end

        PLO_nonsymp = xi*PLO_asymp+(1-xi)*PLO_presymp; % Vector of PLO values corresponding to each search angle

        % Find the search angle which minimises the PLO value:
        if min(PLO_nonsymp)>0 % i.e. if contour has not hit zero
            min_id = min(find(PLO_nonsymp==min(PLO_nonsymp))); % Find the index of the point for which the PLO is minimised
            rho1 = rho1(min_id); % Take the rho1 corresponding to this index
            rho2 = rho2(min_id); % Take the rho2 corresponding to this index
            rho1_plot = [rho1_plot rho1]; % Append new rho1 value to vector of previous rho1 values
            rho2_plot = [rho2_plot rho2]; % Append new rho2 value to vector of previous rho2 values
        elseif((min(PLO_nonsymp)==0)&&(max(PLO_nonsymp)>0))
            min_id = floor(mean(find(PLO_nonsymp==0))); % Find indices for which the PLO is zero and go down the middle of this region
            rho1 = rho1(min_id);
            rho2 = rho2(min_id);
            rho1_plot = [rho1_plot rho1];
            rho2_plot = [rho2_plot rho2];
        end
    end
    
    % Store rho_1 and rho_2 values for this contour
    rho1_plot(buffer) = 0; % Pad with buffer to make sure dimensions are consistent (since contours will be different lengths)
    rho2_plot(buffer) = 0; % -"-
    rho1_two_way_search_2 = [rho1_two_way_search_2;rho1_plot]; % Append new rho2 vector to matrix of vectors for previous contours
    rho2_two_way_search_2 = [rho2_two_way_search_2;rho2_plot]; % Append new rho2 vector to matrix of vectors for previous contours

end

fprintf(' Done!\n');

%% COMPUTE STEEPEST DESCENT CONTOURS USING RHO1_STEP and RHO2_STEP: LINEAR SEARCH CONTOURS
fprintf('Computing linear search contours...');

% Compute contours starting from rho2=0 (i.e. along rho1 axis)
rho1_linear_search_1 = [];  
rho2_linear_search_1 = []; % Empty vectors ready to fill with values

for j=0:4:1*rho_max-1 % j indexes the contours: consider each contour one at a time
    
    rho1 = 1*j;
    rho2 = 0;
    rho1_plot = [rho1];
    rho2_plot = [rho2];

    while (rho1<rho_max+1)
        % Create vectors of rho1 and rho2 values achievable using given search angles and search radius
        rho1 = rho1 + rho1_step;
        rho2 = rho2 + rho2_step;
        
        % For each search angle compute PLO for associated rho1 and rho2:

        PLO_asymp = zeros(1,length(rho1_step));
        PLO_presymp = zeros(1,length(rho1_step));
        PLO_symp = zeros(1,length(rho1_step));
        
        % Define a,b,c,d   
        a = ((1-(delta*rho1)./(1+rho1))*beta1)./(g1+(1-(delta*rho1)./(1+rho1))*(beta1+lambda));
        b = ((1-(delta*rho1)./(1+rho1))*lambda)./(g1+(1-(delta*rho1)./(1+rho1))*(beta1+lambda));
        c = ((1-(delta*rho2)./(1+rho2))*beta2)./(g2+(1-(delta*rho2)./(1+rho2))*(beta2+mu));
        d = ((1-(delta*rho1)./(1+rho1))*beta3)./(g1+(1-(delta*rho1)./(1+rho1))*(beta3+nu));

        % Define polynomial coefficients    
        Z4 = d.*(a-d).*xi.*(d-c);
        Z3 = d.*(a-d).*xi.*c.*(1-d) + (d-a-a.*d.*xi+(d.^2).*(a-1+b.*(1-xi)+xi)).*(d-c) - b.*(1-c).*(d.^3).*(1-xi);
        Z2 = (d-a-a.*d.*xi+(d.^2).*(a-1+b.*(1-xi)+xi)).*(c.*(1-d)) + (d.*(d-2.*a-1)+2.*a).*(d-c);
        Z1 = (d.*(d-2.*a-1)+2.*a).*c.*(1-d) -a.*(1-d).^2.*((d-c));
        Z0 = -a.*(1-d).^2.*(c.*(1-d));

        for i=1:length(rho1_step)
            poly = [Z4(i) Z3(i) Z2(i) Z1(i) Z0(i)];
            poly_roots = roots(poly); % Find roots of polynomial
            
            % Obtain the minimal non-negative real root:
            poly_roots = real(poly_roots(abs(imag(poly_roots))<1e-5));
            poly_roots = poly_roots(poly_roots>=0);
            q_minus = min(poly_roots);

            if(q_minus>=0 && q_minus<1)
                z = q_minus;
                x = (1/(d(i)*(1-xi)))*(1-d(i)*xi*q_minus-(1-d(i))/q_minus);
                y = (1-c(i))*d(i)*z/(d(i)*z-c(i)*(z-1+d(i)));
                PLO_asymp(i) = 1-z;
                PLO_presymp(i) = 1-x;
                PLO_symp(i) = 1-y;
            else 
                PLO_asymp(i) = 0;
                PLO_presymp(i) = 0;
                PLO_symp(i) = 0;
            end           
        end

        PLO_nonsymp = xi*PLO_asymp+(1-xi)*PLO_presymp; % Vector of PLO values corresponding to each search angle

        if min(PLO_nonsymp)>0   % i.e. if contour has not hit zero and should therefore keep going
            min_id = min(find(PLO_nonsymp==min(PLO_nonsymp)));  % Find the index of the point for which the PLO is minimised
            rho1 = rho1(min_id);    % Take the rho1 corresponding to this index
            rho2 = rho2(min_id);    % Take the rho2 corresponding to this index
            rho1_plot = [rho1_plot rho1];   % Append new rho1 value to vector of previous rho1 values
            rho2_plot = [rho2_plot rho2];   % Append new rho2 value to vector of previous rho2 values
        elseif((min(PLO_nonsymp)==0)&&(max(PLO_nonsymp)>0)) 
            min_id = floor(mean(find(PLO_nonsymp==0))); % Find indices for which the PLO is zero and go down the middle of this region
            rho1 = rho1(min_id);
            rho2 = rho2(min_id);
            rho1_plot = [rho1_plot rho1];
            rho2_plot = [rho2_plot rho2];
        end
    end
    
    % Store rho_1 and rho_2 values for this contour
    rho1_plot(buffer) = 0;  % Pad with buffer to make sure dimensions are consistent
    rho2_plot(buffer) = 0;  % -"-
    rho1_linear_search_1 = [rho1_linear_search_1;rho1_plot]; % Append new rho1 vector to matrix of vectors for previous contours
    rho2_linear_search_1 = [rho2_linear_search_1;rho2_plot]; % Append new rho2 vector to matrix of vectors for previous contours

end

% Compute contours starting from rho1=0 (i.e. along rho2 axis): comments as above
rho1_linear_search_2 = [];
rho2_linear_search_2 = []; % Empty vectors ready to fill with values

for j=0:4:1*rho_max-1 % j indexes the contours: consider each contour one at a time
   
    rho1 = 0;
    rho2 = 1*j;
    rho1_plot = [rho1];
    rho2_plot = [rho2];

    while (rho1<rho_max+1)
        % Create vectors of rho1 and rho2 values achievable using given search angles and search radius
        rho1 = rho1 + rho1_step;
        rho2 = rho2 + rho2_step;
        
        % For each search angle compute PLO for associated rho1 and rho2:

        PLO_asymp = zeros(1,length(rho1_step));
        PLO_presymp = zeros(1,length(rho1_step));
        PLO_symp = zeros(1,length(rho1_step));
        
        % Define a,b,c,d   
        a = ((1-(delta*rho1)./(1+rho1))*beta1)./(g1+(1-(delta*rho1)./(1+rho1))*(beta1+lambda));
        b = ((1-(delta*rho1)./(1+rho1))*lambda)./(g1+(1-(delta*rho1)./(1+rho1))*(beta1+lambda));
        c = ((1-(delta*rho2)./(1+rho2))*beta2)./(g2+(1-(delta*rho2)./(1+rho2))*(beta2+mu));
        d = ((1-(delta*rho1)./(1+rho1))*beta3)./(g1+(1-(delta*rho1)./(1+rho1))*(beta3+nu));

        % Define polynomial coefficients    
        Z4 = d.*(a-d).*xi.*(d-c);
        Z3 = d.*(a-d).*xi.*c.*(1-d) + (d-a-a.*d.*xi+(d.^2).*(a-1+b.*(1-xi)+xi)).*(d-c) - b.*(1-c).*(d.^3).*(1-xi);
        Z2 = (d-a-a.*d.*xi+(d.^2).*(a-1+b.*(1-xi)+xi)).*(c.*(1-d)) + (d.*(d-2.*a-1)+2.*a).*(d-c);
        Z1 = (d.*(d-2.*a-1)+2.*a).*c.*(1-d) -a.*(1-d).^2.*((d-c));
        Z0 = -a.*(1-d).^2.*(c.*(1-d));

        for i=1:length(rho1_step)
            poly = [Z4(i) Z3(i) Z2(i) Z1(i) Z0(i)];
            poly_roots = roots(poly); % Find roots of polynomial
            
            % Obtain the minimal non-negative real root:
            poly_roots = real(poly_roots(abs(imag(poly_roots))<1e-5));
            poly_roots = poly_roots(poly_roots>=0);
            q_minus = min(poly_roots);

            if(q_minus>=0 && q_minus<1)
                z = q_minus;
                x = (1/(d(i)*(1-xi)))*(1-d(i)*xi*q_minus-(1-d(i))/q_minus);
                y = (1-c(i))*d(i)*z/(d(i)*z-c(i)*(z-1+d(i)));
                PLO_asymp(i) = 1-z;
                PLO_presymp(i) = 1-x;
                PLO_symp(i) = 1-y;
            else 
                PLO_asymp(i) = 0;
                PLO_presymp(i) = 0;
                PLO_symp(i) = 0;
            end           
        end

        PLO_nonsymp = xi*PLO_asymp+(1-xi)*PLO_presymp; % Vector of PLO values corresponding to each search angle

        if min(PLO_nonsymp)>0 % i.e. if contour has not hit zero and should therefore keep going
            min_id = min(find(PLO_nonsymp==min(PLO_nonsymp))); % Find the index of the point for which the PLO is minimised
            rho1 = rho1(min_id); % Take the rho1 corresponding to this index
            rho2 = rho2(min_id); % Take the rho2 corresponding to this index
            rho1_plot = [rho1_plot rho1]; % Append new rho1 value to vector of previous rho1 values
            rho2_plot = [rho2_plot rho2]; % Append new rho2 value to vector of previous rho2 values
        elseif((min(PLO_nonsymp)==0)&&(max(PLO_nonsymp)>0))
            min_id = floor(mean(find(PLO_nonsymp==0))); % Find indices for which the PLO is zero and go down the middle of this region
            rho1 = rho1(min_id);
            rho2 = rho2(min_id);
            rho1_plot = [rho1_plot rho1];
            rho2_plot = [rho2_plot rho2];

        end
    end
    
    % Store rho_1 and rho_2 values for this contour
    rho1_plot(buffer) = 0; % Pad with buffer to make sure dimensions are consistent
    rho2_plot(buffer) = 0; % -"-
    rho1_linear_search_2 = [rho1_linear_search_2;rho1_plot]; % Append new rho1 vector to matrix of vectors for previous contours
    rho2_linear_search_2 = [rho2_linear_search_2;rho2_plot]; % Append new rho2 vector to matrix of vectors for previous contours

end

fprintf(' Done!\n');

fprintf('Plotting...\n');

%% Plot PLO with probability isolines (e.g. Figs 2B,C,D)
figure(1); 

% Plot heat map of PLO matrix
clims = [min(min(M_nonsymp)),max(max(M_nonsymp))];
imagesc(flip(M_nonsymp'),clims); hold on;

% Add colorbar
colbar = colorbar;
ylabel(colbar, 'Probability of a local outbreak (p)');
set(colbar,'Fontsize',16,'linewidth',1);

% Plot contours of constant probability
levels = [0.6 0.625 0.65 0.7];
[C,h] = imcontour(flip(M_nonsymp'),levels); 
h.LineColor = [1 0 0]; h.LineStyle = ':'; h.LineWidth = 2;
mylabels = clabel(C,h,'manual','color','k','Fontsize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end;

xlabel({'Surveillance targeting nonsymptomatic hosts (\rho_1)'});
xticklabels = ({'0'; '4'; '8'; '12'; '16'; '20'});
xticks = linspace(1, M_size+1, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels);

ylabel({'Surveillance targeting symptomatic hosts (\rho_2)'})
yticklabels = ({'0'; '4'; '8'; '12'; '16'; '20'});
yticks = linspace(1, M_size+1, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels));

set(gca,'Fontsize', 16,'linewidth',1); box off; pbaspect([1.2 1 1]); 

%% Plot PLO matrix with smooth radial search contours (e.g. Fig 3A)
figure(2); 

% Plot heat map of PLO matrix
clims = [min(min(M_nonsymp)),max(max(M_nonsymp))];
imagesc(flip(M_nonsymp'),clims); hold on;

% Add colorbar
colbar = colorbar;
ylabel(colbar, 'Probability of a local outbreak (p)');
set(colbar,'Fontsize',16,'linewidth',1);

% Plot radial search contours starting from rho2=0
for k = 1:length(rho1_radial_search_1(:,1)) % k indexes the contours themselves (so in this loop we consider each contour in turn)
    xplot = rho1_radial_search_1(k,:);      % Vector of rho1 values for this contour
    i1 = find(xplot,1,'last');              % Finds the last index which corresponds to a non-zero element in xplot
    if length(i1)==0; i1=1; end             % If all entries are 0, set i1=1
    xplot = xplot(1:i1);                    % Removes trailing zeros from rho1 vector
    yplot = rho2_radial_search_1(k,:);      % Vector of rho2 values for this contour
    i2 = find(yplot,1,'last');              % Finds the last index which corresponds to a non-zero element in yplot
    if length(i2)==0; i2=1; end             % If all entries are 0, set i2=1;
    yplot = yplot(1:i2);                    % Removes trailing zeros from rho2 vector
    plot(((M_size+1)/rho_max)*xplot+1,M_size+1-((M_size+1)/rho_max)*yplot,'color',[1 1 1],'linewidth',2)
end

% Plot radial search contours starting from rho1=0
for k = 1:length(rho1_radial_search_2(:,1)) % k indexes the contours themselves (so in this loop we consider each contour in turn)
    xplot = rho1_radial_search_2(k,:);      % Vector of rho1 values for this contour
    i1 = find(xplot,1,'last');              % Finds the last index which corresponds to a non-zero element in xplot
    if length(i1)==0; i1=1; end             % If all entries are 0, set i1=1
    xplot = xplot(1:i1);                    % Removes trailing zeros from rho1 vector
    yplot = rho2_radial_search_2(k,:);      % Vector of rho2 values for this contour
    i2 = find(yplot,1,'last');              % Finds the last index which corresponds to a non-zero element in yplot
    if length(i2)==0; i2=1; end             % If all entries are 0, set i2=1;
    yplot = yplot(1:i2);                    % Removes trailing zeros from rho2 vector
    plot(((M_size+1)/rho_max)*xplot+1,M_size+1-((M_size+1)/rho_max)*yplot,'color',[1 1 1],'linewidth',2)
end

xlabel({'Surveillance targeting nonsymptomatic hosts (\rho_1)'}, 'Fontsize', 16);
xticklabels = ({'0'; '4'; '8'; '12'; '16'; '20'});
xticks = linspace(1, M_size+1, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

ylabel({'Surveillance targeting symptomatic hosts (\rho_2)'}, 'Fontsize', 16)
yticklabels = ({'0'; '4'; '8'; '12'; '16'; '20'});
yticks = linspace(1, M_size+1, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels), 'Fontsize', 16);

set(gca,'Fontsize', 16,'linewidth',1); box off; pbaspect([1.2 1 1]); 

%% Plot 'grads' matrix with two way search contour starting from rho1=rho2=0 (e.g. Fig 3B)
figure(3); 

% Plot heat map of 'grads' matrix
colormap([0.5 0.7 1;0.9 0.9 0.9;0.5 0.9 0.5]);
clims2=[0 3];
imagesc(flip(grads'),clims2); hold on

% Plot two way search contour starting from rho1=rho2=0
xplot = rho1_two_way_search_1(1,:);  % Vector of rho1 values for this contour      
i1 = find(xplot,1,'last');           % Finds the last index which corresponds to a non-zero element in xplot    
if length(i1)==0; i1=1; end          % If all entries are 0, set i1=1                                
xplot = xplot(1:i1);                 % Removes trailing zeros from rho1 vector
yplot = rho2_two_way_search_1(1,:);  % Vector of rho2 values for this contour
i2 = find(yplot,1,'last');           % Finds the last index which corresponds to a non-zero element in yplot
if length(i2)==0; i2=1; end          % If all entries are 0, set i2=1; 
yplot = yplot(1:i2);                 % Removes trailing zeros from rho2 vector
thisplot = plot((M_size/rho_max)*xplot+1,M_size-(M_size/rho_max)*yplot,'color',[1 1 1],'linewidth',3)

xlabel({'Surveillance targeting nonsymptomatic hosts (\rho_1)'}, 'Fontsize', 16);
xticklabels = ({'0'; '4'; '8'; '12'; '16'; '20'});
xticks = linspace(1, M_size, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

ylabel({'Surveillance targeting symptomatic hosts (\rho_2)'}, 'Fontsize', 16)
yticklabels = ({'0'; '4'; '8'; '12'; '16'; '20'});
yticks = linspace(1, M_size, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels), 'Fontsize', 16);

text1 = text(22,22,{'Targeting','nonsymptomatic','hosts makes','more','difference'});
text2 = text(65,60,{'Targeting symptomatic hosts', 'makes more difference'});
text1.FontSize = 16; text1.HorizontalAlignment = 'center';
text2.FontSize = 16; text2.HorizontalAlignment = 'center';

set(gca,'Fontsize', 16,'linewidth',1); box off; pbaspect([1.2 1 1]); 

%% Plot PLO matrix with fixed total effort contours and optimal strategies (e.g. Fig 3C)
figure(4);

% Plot heat map of PLO matrix
clims = [min(min(M_nonsymp)),max(max(M_nonsymp))];
imagesc(flip(M_nonsymp'),clims); hold on;

% Add colorbar
colbar = colorbar;
ylabel(colbar, 'Probability of a local outbreak (p)');
set(colbar,'Fontsize',16,'linewidth',1);

% Plot linear search contour starting from rho1=rho2=0
xplot = rho1_linear_search_1(1,:);      % Vector of rho1 values for this contour
i1 = find(xplot,1,'last');              % Finds the last index which corresponds to a non-zero element in xplot
if length(i1)==0; i1=1; end             % If all entries are 0, set i1=1
xplot = xplot(1:i1);                    % Removes trailing zeros from rho1 vector
yplot = rho2_linear_search_1(1,:);      % Vector of rho2 values for this contour
i2 = find(yplot,1,'last');              % Finds the last index which corresponds to a non-zero element in yplot
if length(i2)==0; i2=1; end             % If all entries are 0, set i2=1;
yplot = yplot(1:i2);                    % Removes trailing zeros from rho2 vector
plot(((M_size+1)/rho_max)*xplot+1,M_size+1-((M_size+1)/rho_max)*yplot,'color',[1 1 1],'linewidth',3)

% For each 'total effort contour' considered, find the point along the
% contour at which the PLO is minimised

for j=1:length(total_effort_vec)
    n=total_effort_vec(j); % Define 'total effort'
    rho1_vec = 0:0.01:n;
    rho2_vec = n-rho1_vec; % These values of rho_1, rho_2 make up the contour
    
    PLO_fixed_effort = zeros(1,length(rho1_vec));
    PLO_asymp = zeros(1,length(rho1_vec));
    PLO_presymp = zeros(1,length(rho1_vec));
    PLO_symp = zeros(1,length(rho1_vec));
    
    % Define a,b,c,d   
    a = ((1-(delta*rho1_vec)./(1+rho1_vec))*beta1)./(g1+(1-(delta*rho1_vec)./(1+rho1_vec))*(beta1+lambda));
    b = ((1-(delta*rho1_vec)./(1+rho1_vec))*lambda)./(g1+(1-(delta*rho1_vec)./(1+rho1_vec))*(beta1+lambda));
    c = ((1-(delta*rho2_vec)./(1+rho2_vec))*beta2)./(g2+(1-(delta*rho2_vec)./(1+rho2_vec))*(beta2+mu));
    d = ((1-(delta*rho1_vec)./(1+rho1_vec))*beta3)./(g1+(1-(delta*rho1_vec)./(1+rho1_vec))*(beta3+nu));

    % Define polynomial coefficients    
    Z4 = d.*(a-d).*xi.*(d-c);
    Z3 = d.*(a-d).*xi.*c.*(1-d) + (d-a-a.*d.*xi+(d.^2).*(a-1+b.*(1-xi)+xi)).*(d-c) - b.*(1-c).*(d.^3).*(1-xi);
    Z2 = (d-a-a.*d.*xi+(d.^2).*(a-1+b.*(1-xi)+xi)).*(c.*(1-d)) + (d.*(d-2.*a-1)+2.*a).*(d-c);
    Z1 = (d.*(d-2.*a-1)+2.*a).*c.*(1-d) -a.*(1-d).^2.*((d-c));
    Z0 = -a.*(1-d).^2.*(c.*(1-d));

    for i=1:length(rho1_vec)
        poly = [Z4(i) Z3(i) Z2(i) Z1(i) Z0(i)];
        poly_roots = roots(poly); % Find polynomial roots
        
        % Find minimal non-negative real solution
        poly_roots = real(poly_roots(abs(imag(poly_roots))<1e-5));
        poly_roots = poly_roots(poly_roots>=0);
        q_minus = min(poly_roots);

        if(q_minus>=0 && q_minus<1)
            z = q_minus;
            x = (1/(d(i)*(1-xi)))*(1-d(i)*xi*q_minus-(1-d(i))/q_minus);
            y = (1-c(i))*d(i)*z/(d(i)*z-c(i)*(z-1+d(i)));
            PLO_asymp(i) = 1-z;
            PLO_presymp(i) = 1-x;
            PLO_symp(i) = 1-z;
        else 
            PLO_asymp(i) = 0;
            PLO_presymp(i) = 0;
            PLO_symp(i) = 0;
        end           
    end

    PLO_fixed_effort = xi*PLO_asymp+(1-xi)*PLO_presymp;
    
    [M I]=min(PLO_fixed_effort); % Find the point at which the PLO is minimised
    rho1_index = rho1_vec(I);    % Extract the corresponding rho_1 value
    rho2_index = rho2_vec(I);    % Extract the corresponding rho_1 value
    
    % Plot fixed effort contour
    xplot = rho1_vec; yplot = rho2_vec;        
    plot((M_size/rho_max)*xplot+1,M_size-(M_size/rho_max)*yplot+1,'color',[1 0 0],'linewidth',2,'linestyle',':')
   
    % Plot minimum point
    rho1_plot = 1 + (M_size/rho_max)*rho1_index;
    rho2_plot = 1 + M_size-(M_size/rho_max)*rho2_index;
    redplot = plot(rho1_plot,rho2_plot,'r.','markersize',30);

end

xlabel({'Surveillance targeting nonsymptomatic hosts (\rho_1)'}, 'Fontsize', 16);
xticklabels = ({'0'; '4'; '8'; '12'; '16'; '20'});
xticks = linspace(1, M_size+1, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

ylabel({'Surveillance targeting symptomatic hosts (\rho_2)'}, 'Fontsize', 16)
yticklabels = ({'0'; '4'; '8'; '12'; '16'; '20'});
yticks = linspace(1, M_size+1, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels), 'Fontsize', 16);

set(gca,'Fontsize', 16,'linewidth',1); box off; pbaspect([1.2 1 1]); 

%% Plot PLO matrix with fixed probability contours and optimal strategies (e.g. Fig 3D)
figure(5); 

% Plot heat map of PLO matrix
clims = [min(min(M_nonsymp)),max(max(M_nonsymp))];
imagesc(flip(M_nonsymp'),clims); hold on;

% Add colorbar
colbar = colorbar;
ylabel(colbar, 'Probability of a local outbreak (p)');
set(colbar,'Fontsize',16,'linewidth',1);

% Plot contours of constant probability
levels = [0.6 0.625 0.65 0.68 0.69 0.7];
[C,h] = imcontour(flip(M_nonsymp'),levels); 
h.LineColor = [1 0 0]; h.LineStyle = ':'; h.LineWidth = 2;
mylabels = clabel(C,h,'manual','color','k','FontSize',16);
for i=1:length(mylabels); mylabels(i).Color = [1 1 1]; end;

% Plot linear search contour starting from rho1=rho2=0
xplot = rho1_linear_search_1(1,:);      % Vector of rho1 values for this contour
i1 = find(xplot,1,'last');              % Finds the last index which corresponds to a non-zero element in xplot
if length(i1)==0; i1=1; end             % If all entries are 0, set i1=1
xplot = xplot(1:i1);                    % Removes trailing zeros from rho1 vector
yplot = rho2_linear_search_1(1,:);      % Vector of rho2 values for this contour
i2 = find(yplot,1,'last');              % Finds the last index which corresponds to a non-zero element in yplot
if length(i2)==0; i2=1; end             % If all entries are 0, set i2=1;
yplot = yplot(1:i2);                    % Removes trailing zeros from rho2 vector
plot(((M_size+1)/rho_max)*xplot+1,M_size+1-((M_size+1)/rho_max)*yplot,'color',[1 1 1],'linewidth',3)

% Identify and plot points on probability contours which minimise PLO
i=1; contour_start = 1; 
while i<=length(h.LevelList) % For each level in list...
    contour_length = C(2,contour_start); % The number of entries which make up this contour
    mycontour = C(:,[contour_start+1:contour_start+contour_length]); % Extract the elements corresponding to this contour
    mycontours{i} = mycontour; % Store the elements corresponding to this contour
    mycontour2 = [mycontour(1,:);M_size+2-mycontour(2,:)]; % Convert into (rho_1,rho_2) coordinates
    total_effort = sum(mycontour2); % Compute total effort at each point along the contour
    [M,I] = min(total_effort); % Find the minimum total effort and its index
    rho1_plot = mycontour(1,I); % Find corresponding value of rho_1
    rho2_plot = mycontour(2,I); % Find corresponding value of rho_2
    if rho1_plot<M_size+1 && rho2_plot>1
        redplot = plot(rho1_plot,rho2_plot,'r.','markersize',30); % Plot minimum point
    end
    contour_start = contour_start+contour_length+1;
    i=i+1;
end

xlabel({'Surveillance targeting nonsymptomatic hosts (\rho_1)'}, 'Fontsize', 16);
xticklabels = ({'0'; '4'; '8'; '12'; '16'; '20'});
xticks = linspace(1, M_size+1, numel(xticklabels));
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels, 'Fontsize', 16);

ylabel({'Surveillance targeting symptomatic hosts (\rho_2)'}, 'Fontsize', 16)
yticklabels = ({'0'; '4'; '8'; '12'; '16'; '20'});
yticks = linspace(1, M_size+1, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels), 'Fontsize', 16);

set(gca,'Fontsize', 16,'linewidth',1); box off; pbaspect([1.2 1 1]); 

%%
fprintf('FINISHED\n\n');

