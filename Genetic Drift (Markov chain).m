% Simulating Genetic Drift with a Markov model

% Q1:
Q1 = 8/20;

% Q2:
Q2 = (8/20)^2;

% Q3:
Q3 = binopdf(8, 20, .4);

% Q4:
Q4 = 1 - Q3;

% Q5:
Q5 = 1 - binocdf(8, 20, 8/20);
% Probability of decreasing = 4.16 > probabiliy of increasing = 4.044

% Q6:
chains = sim_wright_fisher(20, 10, 15, 100);
plot(chains)
xticks(1:16);
xticklabels(cellstr(num2str((0:15)')))
ylabel('Generation')
ylabel('Copies of Allele A')

% Q7:
figure;
counts_generation_1 = chains(2, :);
histogram(counts_generation_1);
xlabel('Copies of Allele A')
ylabel('Counts')

% Q8:
figure;
counts_generation_15 = chains(16, :);
histogram(counts_generation_15);
xlabel('Copies of Allele A')
ylabel('Counts')
% 14 reached 0 copies, 21 reached 20 copies for a total of 35 fixed copies.

% Q9a:
N = 10;
P = zeros(2*N+1, 2*N+1);
for i = 0:2*N
    P(i+1, :) = binopdf(0:2*N, 2*N, i/(2*N));
end

% Q9b:
pi0 = zeros(1, 2*N+1);
pi0(N+1) = 1;

% Q9c:
figure;
pi1 = pi0*P;
bar(pi1)
xticks(1:21);
xticklabels(cellstr(num2str(0:20)'))
xlabel('Copies of Allele A')
ylabel('Counts')

% Q9d:
figure;
pi15 = pi0*P^15;
bar(pi15)
xticks(1:21);
xticklabels(cellstr(num2str(0:20)'))
xlabel('Copies of Allele A')
ylabel('Counts')

% Q9e:
% The simulations from Q7 and Q8 match up well to the theoretical
% distribution displayed in Q9c and Q9d.

% Q10:
figure;
buri_data = csvread('buri_data.csv');
bar3(buri_data(2:end,:), 'w')
xlabel('Copies of Allele')
yticks(1:length(buri_data(:,1))-1)
xticks(1:33)
xticklabels(cellstr(num2str((0:32)')))
ylabel('Generation')
zlabel('Number of Observations')
% The distribution is consistent with the theory of genetic drift as most
% populations converged towards having just a single allele uniform throughout
% the population which is representative of populations diverging
% genetically through random chance.

% Q11:
[chains_buri, counts_buri] = sim_wright_fisher(32, 16, 19, 107);
figure;
bar3(counts_buri(2:end,:), 'w')
xlabel('Copies of Allele')
yticks(1:length(counts_buri(:,1))-1)
xticks(1:33)
xticklabels(cellstr(num2str((0:32)')))
ylabel('Generation')
zlabel('Number of Observations')
% The simulations are much less extreme in the growth of the recursive
% states (single allele populations), but they still do generally converge
% towards the recursive states overall there is just a more even
% distribution between.

