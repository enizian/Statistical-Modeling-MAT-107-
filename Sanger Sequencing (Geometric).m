
% Sanger Sequencing Simulation with Geometric Distribution

k = 1:10;
pdf_1 = geopdf(k-1, 1-.95);
pdf_2 = geopdf(k-1, 1-.7);
pdf_3 = geopdf(k-1, 1-.5);

figure;
hold on
plot(k, pdf_1, 'bs')
plot(k, pdf_2, 'rs')
plot(k, pdf_3, 'gs')

ylim([0, 1])
ylabel('Probability')
xlim([0.5 10.5])
xticks(1:10)
xlabel('k')
legend('p=0.9')
title('Probability of Ending at k-th T in S')

% Q2: p = 0.95


% Q3: As k becomes larger the probability of sequencing long fragments
% decreases substantially. This inhibits the abiliy to sequence by Sanger very
% large DNA fragments as there will be too few longer fragments (of the larger DNA fragment)
% to give a reliable read on the sequence of the DNA fragment.


% Q4:
figure;
hold on
p = 0:0.001:1;

n_exp_1 = 1./geopdf(1-1, 1-p);
n_exp_2 = 1./geopdf(2-1, 1-p);
n_exp_3 = 1./geopdf(3-1, 1-p);

plot(p, n_exp_1)
plot(p, n_exp_2)
plot(p, n_exp_3)

xlim([0 1])
xlabel('Proportion of Good Nucleotides')
ylabel('Expected Number of Experiments')
ylim([0 10])
legend('k=1', 'k=2', 'k=3')


% Q5: When k=1, the expected number of experiments increases in an
% exponential manner as the proportion of good nucleotides increases.
% When k=2, the expected number of experiments follows an almost parabolic
% curve where the number of experiments is high when the proportion of good
% nucleotides is low and high.
% When k=3, the expected number of experiments follows a similar
% "parabolic" curve as k=2, but more narrow and localized around a minimum
% of large p value ~7.


% Q6: These points represent the minimum number of expected experiments for
% respective values of k = 1, 2, 3
k = 1:3;
p_opt = (k-1)./k;

n_opt = [1./geopdf(0, 1-p_opt(1)),...
         1./geopdf(1, 1-p_opt(2)),...
         1./geopdf(2, 1-p_opt(3))];

plot(p_opt, n_opt, 'k*')

% Q7: p=.9 appears to work better for this experiment because there are far
% more fragments of longer lengths which allows the latter half of the DNA
% molecule sequence to be reliably read in gel electrophoresis (or otherwise).
p = 0.8;
Nk = 30;
N_experiments = 10000;
outcomes = NaN(N_experiments, 1);
for i = 1:N_experiments
 outcomes(i) = min([geornd(1-p)+1, Nk+1]);
end
%display(outcomes)
figure;
histogram(outcomes)
xlabel('k')
ylabel('Number of Fragments')
title('p=.8')

p = 0.9;
Nk = 30;
N_experiments = 10000;
outcomes = NaN(N_experiments, 1);
for i = 1:N_experiments
 outcomes(i) = min([geornd(1-p)+1, Nk+1]);
end
%display(outcomes)
figure;
histogram(outcomes)
xlabel('k')
ylabel('Number of Fragments')
title('p=.9')

% Q8: The distribution of fragments changes to a relatively uniform amount
% across every fragment except the whole sequence (k=31) which is far
% larger than rest, indicating in an experiment using the optimal p value
% the majority of fragments sequenced are the whole sequence. However, this
% uniform amount across the rest of the fragments allows for more accurate
% sequencing in gel electrophoresis (or otherwise).
Nk = 30;
p = (Nk-1)/Nk;
N_experiments = 10000;
outcomes = NaN(N_experiments, 1);
for i = 1:N_experiments
 outcomes(i) = min([geornd(1-p)+1, Nk+1]);
end
%display(outcomes)
figure;
histogram(outcomes)
hold on
xlabel('k')
ylabel('Number of Fragments')
title(['p=p_opt (' num2str(p) ')'])

% Q9: All of the fragments in the distribution from Question 8 are detected
% at D=.01. A challenge that emerges from longer lengths of DNA is ensuring
% that all fragments are able to reach the detection limit because as the
% number of fragments grow due to the increased size of the DNA, the
% percentage each fragment makes up decreases with the same number of
% experiments (so the number of experiments would necessarily have to
% increase to account for this increase in length).
plot([0,32], [100,100], '--r')
