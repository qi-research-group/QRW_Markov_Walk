
nstates = 9; % (an uneven numbet of states)

% a completely random walk means things should just stay IN THE MIDDLE.
% so the way to do that is:
diag_H_base = -((nstates-1)/2):1:((nstates-1)/2);

% BUT, if we are doing some "moving to the direction that we want"; it may
% be reasonable that a hamiltonian that conforms/increases reliability
% has a slightly other diag:

% let's just use the framework of the random walk QModel itself
% and change the H according to the ehm; events and responses that
% happened.
% increment and decrement are model parameters
% as well as the drift

% and when there is a disagreement, we can model that by pushing the
% reliability down:


% increment_trust and increment_distrust are our MODEL PARAMETERS.


% base level of trust: this will just be the offset; We will take the
% rating % AFTER PRACTICE BLOCK - and just subtract all ratings with that
% so the results are around 0.
% then the H might do the actual trick.

dhist = dir('../../history_files/history_*.txt');
drating = dir('../../history_files/ratings_*.txt');

h = load([dhist(1).folder filesep dhist(1).name]);
h = h(h<100);

r = load([drating(1).folder filesep drating(1).name]);
begin_rating = r(2);


% let's make S0 the start rating; and let's put that in the middle pls:
% #sum(np.square(S0)) # sum of all elements of S0 squared has to be equal 1
% well, let's set the S0 to:
% S0 = zeros(1, nstates);
S0 = [0 0 0 0 0 0 0 1 0];
% S0((nstates+1)/2) = 1;
S0 = S0 / (S0*S0'); % a but superfluous, but for formality, needs to be == 1.
S0 = reshape(S0, numel(S0), 1);

% let's begin the simulation.
St = S0;
R = [0];
t = [0]; % time = 0 at the beginning;
diag_H = diag_H_base;
% diag_H = diag_H*0;
timestep = 4;


sigmasquared = 1/10;
mv = 1/10;
delta_increment = 1; %  0.1;
delta_decrement = 1;  %0.2;

for i=1:10 % numel(h) % i.e. the history:
    

    % get the hamiltonian based on match ve mismatch
    switch h(i)
        case {4, 6} % which is HRAIR, and HFAIF
            % disp('match');
            if max(diag_H) < (nstates - delta_increment)
                diag_H = diag_H * 1;
            end
        case {5, 7} % which is HRAIF, and HFAIR
            % disp('mismatch');
            if min(diag_H) > (-nstates + delta_increment)
                diag_H = diag_H * -1; %- delta_decrement;
            end

    end
    disp(diag_H);
         

            
    H = make_H_v2(mv * diag_H, sigmasquared, 0);

 
    t(end+1) = t(end) + timestep; % for artuments sake; make it 4 ; but can be sibjected to change.
    disp(H);
    
    
    % evolve things:
    U = expm(-1j * timestep * H);
    St = U*St;
    
    
    reliability_as_per_qmodel = diag_H_base * (abs(St).^2);
    
    R(end+1) = reliability_as_per_qmodel;
end


figure;plot(t, R)


