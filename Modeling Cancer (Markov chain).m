

% Modeling Cancer with Markov Chains

% Q1: Columns should be all zeros.


% Q2: 
P = xlsread('cancer_model.xls');
answers_Q2 = (P(:, 2)); %all zeros -> agrees with Q1 answer


% Q3:
figure;
bar(P(23, :));
xlabel('Site')
ylabel('Transition Probability')
xticks(1:50)
title('Transition Probabilities from the Lung')
% Lung cancer is most likely to migrate to site 24, the regional lymph nodes.

% Q4:
pi0 = zeros(1, 50);
pi0(23) = 1;
pi100 = pi0*P^1000;
display(pi100')

% Q5:
chain = sim_MC(P, 23, 100);
counts = histc(chain, 1:50);
display(counts)
% The regional lymph nodes do tend to appear at the highest frequency.
% ~20 times per 100 time steps

% Q6:
N_sim = 10000;
fp_heart = NaN(N_sim, 1);
fp_pancreas = NaN(N_sim, 1);
fp_lymph = NaN(N_sim, 1);

for i = 1:N_sim
    chain = sim_MC(P, 23, 500);
    fp_heart(i) = find(chain==17, 1)-1;
    fp_pancreas(i) = find(chain==28, 1)-1;
    fp_lymph(i) = find(chain==24, 1)-1;
end

mfpt_heart = mean(fp_heart);
mfpt_pancreas = mean(fp_pancreas);
mfpt_lymph = mean(fp_lymph);
% The mean first-passage time for the lymph nodes is much lower than that
% of the heart and pancreas.

% Q7:
Q28 = P;
Q28(28, :) = 0;
Q28(:, 28) = 0;
T28 = (eye(50)-Q28)\ones(50,1);
answers.pancreas = T28(23);

Q24 = P;
Q24(24, :) = 0;
Q24(:, 24) = 0;
T24 = (eye(50)-Q24)\ones(50,1);
answers.lymph = T24(23);

display(answers)
