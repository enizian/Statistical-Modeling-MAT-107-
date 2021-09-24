

% Simulating RNA-Seq with a Poisson Distribution


% Q1:
samples = poissrnd(10, 1000, 1);
mean_samples = mean(samples);
display(mean_samples);
var_samples = var(samples);
display(var_samples);

% Q2:
k = 0:1:20;
n = 10000;
p = .001;

probs_binomial = binopdf(k, n, p);
probs_poisson = poisspdf(k, n*p);

figure;
bar(k, [probs_binomial' probs_poisson'])
ylabel('Probability')
legend('Binomial', 'Poisson')

% Q3:
ti = [200 300 1000 500 450];
ni = [5 3 5 7 10];
np = 200;

counts = rna_seq_sim(ti, ni, np);

% Q4:
p_genes = (counts./(ti - 5))/(sum(counts./(ti-5)));


% Q5: The values calculated in Q4 are relatively close in value to the
% underlying relative abundance.


% Q6:
Ns = 1000;
counts_A = NaN(Ns, 1);
for i = 1:Ns
    counts = rna_seq_sim(ti, ni, np);
    counts_A(i) = counts(1);
end
figure;
histogram(counts_A)

% Q7: The data from Question 6 agrees quite well with the poisson
% distribution overlay generated.
Ns = 1000;
counts_A = NaN(Ns, 1);
for i = 1:Ns
    counts = rna_seq_sim(ti, ni, np);
    counts_A(i) = counts(1);
end
figure;
hold on
histogram(counts_A, 'Normalization', 'pdf');

mc = mean(counts_A);
k = 0:1:10;
plot(k, poisspdf(k, mc));


% Q8: 
Ns = 1000;
counts_A = NaN(Ns, 1);
for i = 1:Ns
    counts = rna_seq_sim(ti, poissrnd(ni), np);
    counts_A(i) = counts(1);
end
figure;
hold on
histogram(counts_A)


% Q9: The Poisson distribution does not fit the data well and appears to
% overestimate outcomes for experiments with counts 2-5 and underestimates 
% outcomes 0 and above 9.
Ns = 1000;
counts_A = NaN(Ns, 1);
for i = 1:Ns
    counts = rna_seq_sim(ti, poissrnd(ni), np);
    counts_A(i) = counts(1);
end
figure;
hold on
histogram(counts_A, 'Normalization', 'pdf');

mc = mean(counts_A);
k = 0:1:20;
plot(k, poisspdf(k, mc));

% Q10: Negative binomial appears to still overesimate at certain count
% values but overall it is a better fit to the data than Poisson.
parameters = nbinfit(counts_A);
plot(k, nbinpdf(k, parameters(1), parameters(2)));


% Q11: A negative binomial is a more suitable choice for modeling RNA seq
% counts. Poisson and negative binomial distributions are both count-based,
% but negative binomial is better for RNA seq because it does not assume
% variance and expectation to be equal like Poisson. An RNA-seq experiment
% will involve large variation in the counts per gene (and cell) better 
% suited for modeling by negative binomial.

