
% Lab 7: Genomic composition of E.coli modeled with Markov Chains

% Q1:
P = [.90 .10; .48 .52];

% Q2:
P200 = P^200;

% Q3:
pi0_dry = [1, 0];
pi200_dry = pi0_dry*P200;

pi0_wet = [0, 1];
pi200_wet = pi0_wet*P200;

% Distribution converges to the limiting matrix, regardless of initial
% being wet or dry.

% Q4:
[V, D] = eig(P');

% Q5:
display(V(:,1));

% Q6:
display(sum(V(:,1)));
% No, the eigenvector does not sum to 1.

% Q7:
pi_s = V(:,1)/sum(V(:,1));

% Q8: Over the course of a full year with about 1/6 of all days being wet
% and about 5/6 being dry it is probably about accurate. Seasonality would
% not be represented well in this model as the distribution of transitions
% would be even throughout when there should be more or less depending on
% the season. The model may be improved by splitting it into four separate 
% markov chains for each season.


% Q9:
P = [0.9, 0.10; 0.48, 0.52]; % Transitions
T = 365;
X = NaN(1,T);
X(1) = 1;
for i = 2:T
    probs = P(X(i-1),:);
    X(i) = binornd(1, probs(2))+1;
end
%display(X)

% Q10:
figure;
hold on
p_wet = cumsum(X-1)./(1:T);
plot(1:T, p_wet)
plot([1, T], [.1724, .1724], 'color', 'r')  
xlabel('Time (days)')
ylabel('Proportion of Days which were wet')

% Q11:
% a) [.5, .5]
% b)
P2 = [.5, .5; .5, .5];
X2 = NaN(1,T);
X2(1) = 1;
for i = 2:T
    probs = P2(X2(i-1),:);
    X2(i) = binornd(1, probs(2))+1;
end

figure;
hold on
p2_wet = cumsum(X2-1)./(1:T);
plot(1:T, p2_wet, 'r--')
plot([1, T], [.5, .5], 'color', 'r')
plot(1:T, p_wet, 'b--')
plot([1, T], [.1724, .1724], 'color', 'b')  
xlabel('Time (days)')
ylabel('Proportion of Days which were wet')
% In most simulations Davis appears to converge faster than "Random Town."

%Q12:
file = fastaread('Ecoli-k12-genome.fasta');
seq = file.Sequence;
di_nuc_counts = zeros(4);
for i = 1:length(seq)-1
    di_nuc = seq(i:i+1); % Read dinucleotide
    row = nuc_to_index(di_nuc(1)); % Convert letters to index for matrix
    col = nuc_to_index(di_nuc(2)); % Convert letters to index for matrix
    di_nuc_counts(row, col) = di_nuc_counts(row, col) + 1;
end
transition_matrix = di_nuc_counts./sum(di_nuc_counts, 2);
display(transition_matrix);


% Q13:
[V, D] = eig(transition_matrix');
vec = V(:,1)/sum(V(:,1));
display(vec)

% Q14:
file = fastaread('Organism_A.fasta');
seq = file.Sequence;
di_nuc_counts = zeros(4);
for i = 1:length(seq)-1
    di_nuc = seq(i:i+1); % Read dinucleotide
    row = nuc_to_index(di_nuc(1)); % Convert letters to index for matrix
    col = nuc_to_index(di_nuc(2)); % Convert letters to index for matrix
    di_nuc_counts(row, col) = di_nuc_counts(row, col) + 1;
end
transition_matrix_A = di_nuc_counts./sum(di_nuc_counts, 2);
display(transition_matrix_A);

file = fastaread('Organism_B.fasta');
seq = file.Sequence;
di_nuc_counts = zeros(4);
for i = 1:length(seq)-1
    di_nuc = seq(i:i+1); % Read dinucleotide
    row = nuc_to_index(di_nuc(1)); % Convert letters to index for matrix
    col = nuc_to_index(di_nuc(2)); % Convert letters to index for matrix
    di_nuc_counts(row, col) = di_nuc_counts(row, col) + 1;
end
transition_matrix_B = di_nuc_counts./sum(di_nuc_counts, 2);
display(transition_matrix_B);

% Q15:
[VA, ~] = eig(transition_matrix_A);
comp_A = VA(:,1)/sum(VA(:,1));
display(comp_A);

[VB, ~] = eig(transition_matrix_B);
comp_B = VB(:,1)/sum(VB(:,1));
display(comp_B);

% Q16:
% Analysis of the dinucleotide counts for the genomes does allow us to
% distinguish which genome belongs to which organism. transition_matrix_A
% likely belongs to E. atrepeatium as the probability of self-transitioning
% for A and T for Organism_A is much higher than that of Organism_B. 

function index = nuc_to_index(nuc)
switch nuc
 case {'A', 'a'}
 index = 1;
 return
 case {'T', 't'}
 index = 2;
 return
 case {'G', 'g'}
 index = 3;
 return
 case {'C', 'c'}
 index = 4;
 return
end
end
