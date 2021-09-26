% problem: find best x & y within [0,31] where 557917^2 - (x-557917)^2 reaches the greatest
% individual coded as 20-bit binary
% fitness function is the same as objective function: 557917^2 - (x-557917)^2
clear; close all
M = 1000; % amount of individuals
nCode = 20; % bits of code
pc = 0.7; % probability of gene cross
pm = 0.1; % probability of gene variation
MAX_ITER = 100000; % max iteration
CODE_TRANS = [2^19; 2^18; 2^17; 2^16; 2^15; 2^14; 2^13; 2^12; 2^11; 2^10; 2^9; 2^8; 2^7; 2^6; 2^5; 2^4; 2^3; 2^2; 2; 1]; % transform matrix, from binary to decimal
fit_rec = zeros(MAX_ITER, 1);

%% initialization

ppl = zeros(M, nCode); % population of 1000 individuals, with 5 codes each

for iter = 1:100 % 100 iterations for initialzation, picking up 10 best-fit-individuals every time
    t_indv = rand(100, nCode) > 0.5; % code of 100 indivuals
    t_fit = 557917^2 - (t_indv * CODE_TRANS - 557917).^2; % fitness value of 100 individuals
    [~, i] = sort(t_fit); % sort from min to max
    ppl(10*(iter-1)+1:10*(iter-1)+10, :) = t_indv(i(end-9:end), :); % select top 10 individuals as initial
end

fit = 557917^2 - (ppl * CODE_TRANS - 557917).^2; % fitness value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ITERATION = 1:MAX_ITER
%% reproduction
% 轮盘赌的问题：适应度都很高时，概率差不多
prob = fit / sum(fit); % probability
sum_prob = zeros(M, 1); % accumulated probability
sum_prob(1) = prob(1);
t_ppl = zeros(M, nCode);
for i = 2:M
    sum_prob(i) = sum_prob(i-1) + prob(i);
end
r = rand(M,1);
for n = 1:M % number of random choice
    i = 1;
    while(r(n)>sum_prob(i))
        i = i+1;
    end
    t_ppl(n,:) = ppl(i,:);
end
ppl = t_ppl;

%% gene cross

r = randi(M, [floor(M*pc), 1]); % individuals that go gene cross, in pairs
rPos = randi(nCode-1, [floor(M*pc/2),1]); % random cross point
for i = 1:floor(M*pc/2) % every pair of cross: r(2i-1) & r(2i), cross point at rPos(i)
    [ppl(r(2*i-1), (rPos(i)+1):end), ppl(r(2*i), (rPos(i)+1):end)] = deal(ppl(r(2*i), (rPos(i)+1):end), ppl(r(2*i-1), (rPos(i)+1):end)); % swap
end

%% gene variation

r = randi(M, [floor(M*pm), 1]); % individuals that go gene variation
rPos = randi(nCode, [floor(M*pm),1]); % random variation point
for i = 1:floor(M*pm)
    ppl(r(i),rPos(i)) = ~ppl(r(i),rPos(i));
end

%% fitness evaluation

fit = 557917^2 - (ppl * CODE_TRANS - 557917).^2; % fitness value
fit_rec(ITERATION) = mean(fit);
end

figure, plot(fit_rec);