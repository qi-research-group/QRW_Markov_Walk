%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mean conf across time only seven states
clc
clear
%%% initial parameters
ns = 7;
% odd no. evidence states
ws = 3;
% ws = start width
tv = 0:.1:20;
% no time steps
nt = size(tv,2);
Mid = (ns+1)/2;
mv = -(Mid-1):(Mid-1);

%%%% Quantum model
mu = 1;
% drift rate
ap = 1;
% diﬀusion
% build start state
S0 = zeros(ns,1);
S0((Mid-ws):(Mid+ws)) = 1;
S0 = [0 0 1 10 1  0 0]';
% S0((Mid-ws):(Mid)) = 1;
S0 = S0./sqrt(S0'*S0);
% build Hamiltonian
b = mu*mv;
a = ap*ones(ns,1);
H = buildH(a,b,a);
% function given below
% time loop
PM1 = [];
for n=1:nt
t = tv(n);
U = expm(-1i*t*H);
St = U*S0;
Mc = mv*(abs(St).^2);
PM1 = [PM1 ; [n+1 Mc]];
end

%%%% Markov model
mu = .5; % drift rate
var = 2; % diﬀusion
% build start state
S0 = zeros(ns,1);
S0((Mid-ws):(Mid+ws)) = 1;
S0 = S0./sum(S0);
% build intensity matrix
mk = ones(ns,1);
b = -var*mk;
a1 = .5*(var-mu)*mk;
a2 = .5*(var+mu)*mk;
K = buildK(a1,b,a2); % function given below
% time loop
PM2=[];
for n=1:nt
t = tv(n);
T = expm(t*K);
Pt = T*S0;
Mc = mv*Pt;
PM2 = [PM2 ; [n+1 Mc]];
end
%%%%% Plot results
figure;
plot(tv,PM2(:,2),'-o',tv,PM1(:,2),'.-')
xlabel('Time')
ylabel('Mean Conﬁdence')
legend('Markov','quantum')
