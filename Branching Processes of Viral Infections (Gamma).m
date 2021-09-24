
% Branching Processes of Viral Infections (Gamma)


% SARS-COV-1 SECTION

 % Q1:
 X = sim_branching_poisson(1.88, 10);
 figure;
 plot(0:1:10, X);
 xlabel('Generation')
 ylabel('Outbreak Size')
 
 % Q2:
 Xs = NaN(11, 100);
 for i = 1:100
     Xs(:, i) = sim_branching_poisson(1.88, 10);
 end
 figure;
 plot(0:1:10, Xs);
 xlabel('Generation')
 ylabel('Outbreak Size')
 
 % Q3:
 poisson_mean_10 = mean(Xs(11, :));
 
 % Q4:
 poisson_prop_ended_10 = sum(Xs(11,:)==0)/100;
 
 % Q5:
 expected_poisson_total_size_10 = sum(1.88.^(0:1:10));
 
 % Q6: The prediction compares well to the distribution of total infections.
 Xs_total = cumsum(Xs);
 mean_poisson_total_size_10 = mean(Xs_total(11, :));
 
 % Q7a:
 load('secondary_infections.mat');
 nbin_params = nbinfit(counts);
 r = nbin_params(1);
 p = nbin_params(2);
 
 % Q7b:
 figure;
 hold on
 histogram(counts, 'Normalization', 'probability');
 plot(0:1:35, nbinpdf(0:1:35, r, p), 'r');
 xlabel('Secondary Infections')
 ylabel('Probability')
 
 % Q7c: Negative binomial better captures the offspring distribution.
 plot(0:1:35, poisspdf(0:1:35, 1.88), 'b');
 legend('Data', 'NB fit', 'Poisson fit');
 
 % Q8a:
 Xs_NB = NaN(11, 100);
 for i = 1:100
     Xs_NB(:,i) = sim_branching_nbin(r, p, 10);
 end
 nbin_mean_10 = mean(Xs_NB(11, :));
 
 % Q8b:
 nbin_prop_ended_10 = sum(Xs_NB(11, :)==0)/100;
 
 % Q8c:
 Xs_NB_total = cumsum(Xs_NB);
 mean_Xs_NB_total_10 = mean(Xs_NB_total(11,:));

 Q8d: The total outbreak sizes tend to be comparable, but the negative
 binomial outbreaks tend to end at a much higher frequency.


% COVID-19 SECTION

% Q9:
load('hubei_cases.mat');
fun_exp = @(x, days) x(1)*exp(x(2)*days);
params_exp = lsqcurvefit(fun_exp, [1, .5], days_hubei, cases_hubei, [0,0], [2,5]);
figure;
hold on
scatter(days_hubei, cases_hubei, 'ok');
days = 15:0.1:25;
plot(days, fun_exp(params_exp, days), 'r')
% transmissibility = .4539


% Q10:
mu = 7.6;
sigma = 3.4;
alpha_sars = (mu/sigma)^2;
beta_sars = mu/sigma^2;

% Q11: R0 = 13.7785
R0 = 1/laplace_gamma(params_exp(2), alpha_sars, beta_sars);

% Q12: Very high risk for pandemic given its high reproduction number.


% Q13:
mu_si_covid = 4.8;
sigma_si_covid = 2.3;
alpha_covid = (mu_si_covid/sigma_si_covid)^2;
beta_covid = mu_si_covid/sigma_si_covid^2;
R0_updated = 1/laplace_gamma(params_exp(2), alpha_covid, beta_covid);

% Q14: A shorter serial interval leads to faster growth of cases in the
% time domain, so by decreasing the interval the R0 likewise decreased.


% Q15:
load('ca_cases.mat');
params_exp_ca = lsqcurvefit(fun_exp, [1, .5], days_ca, cases_ca, [0, 0], [2, 5]);

% Q16: The reproduction number in California is 2.3411 whereas the
% reproduction number in Hubei is 5.8515
R = 1/laplace_gamma(params_exp_ca(2), alpha_covid, beta_covid);

% Q17: With a much lower reproduction number than that of Hubei, efforts in
% California to attenuate the spread of the virus likely did work.
% California also has a lower population density than Hubei, which likely
% played a part in lowering the reproudction number among other factors.



function X = sim_branching_poisson(lambda, n_max)
 
X = zeros(1, n_max+1); % Initialize infection counts to zero
X(1) = 1; % Single infection at X0
 
for n = 1:n_max
    
    counts = 0;
 
    for i = 1:X(n)
        counts = counts + poissrnd(lambda);
    end
 
    X(n+1) = counts; % Total number of secondary infections
 
end
end 

function X = sim_branching_nbin(r, p, n_max)
 
X = zeros(1, n_max+1); % Initialize infection counts to zero
X(1) = 1; % Single infection at X0
 
for n = 1:n_max
    
    counts = 0;
 
    for i = 1:X(n)
        counts = counts + nbinrnd(r, p);
    end
 
    X(n+1) = counts; % Total number of secondary infections
 
end
end 

function M = laplace_gamma(s, alpha, beta)
M = (beta^alpha)/((s + beta)^alpha);
end