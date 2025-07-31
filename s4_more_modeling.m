% so the H+, with positive gradient, 'explores' cyclically in higher state
% regimes. The H-, with negative gradient, 'explores state space cyclically
% in lower state regimes.

% the initial state sets the bar of this exploration. Cyclical 'energy'
% likely still is maintained in the imaginary component, so I don't know
% how alternating H+ and H- will play into this.

% lets explore offset and slope of gradients.
% so, the slope of the diagonal - matters (more faster,generally, of the
% exploring stuff; the Offset - absolutely NO effect. hmmmm... why...

% crosse-wise diffusion could be also only come into effect once step has
% occupied that slot.

% next try to model with H+ and H-; with the forces in opposite direction
% (or something). Then fit scale-time, and slope, to the behavioural data.



% lets go:


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
S0 = [0 0 0 0 0 0 0.5 1 0.5];
% S0((nstates+1)/2) = 1;

S0 = S0/sqrt(S0*conj(S0'));
S0 = reshape(S0, numel(S0), 1);

% let's begin the simulation.
St = S0;
R = [0];
t = [0]; % time = 0 at the beginning;
% diag_H = diag_H_base;
% diag_H = diag_H*0;

% delta_t = 4;

timestep = 3.5; % (in our seconds
% parameter to be optimized
time_scaling = 1/40; % general rotational shenanigans; a scaling factor
d_discrepancy_factor = 0; % this is how much an oddball will (extra) rotate with the negative H; how impactful is a miss?
H_SLOPE = 1; % another scaling factor; of the forces that are embedded within the H.

state_vec = [-4:4];



sigmasquared = 1;
% mv = 1;
% delta_increment = 1; %  0.1;
% delta_decrement = 1;  %0.2;

for i=1:numel(h) % i.e. the history:
    

    % get the hamiltonian based on match ve mismatch
    switch h(i)
        case {4, 6} % which is HRAIR, and HFAIF
            
            diag_H = [-4:4];
            discrepancy_factor = 1;
            
        case {5, 7} % which is HRAIF, and HFAIR
            
            diag_H = [4:-1:-4];
            discrepancy_factor = 1 + d_discrepancy_factor;

    end
    
    H = make_H_v2(H_SLOPE * diag_H, sigmasquared, 0);
    
    
    t(end+1) = t(end) + timestep; % for artuments sake; make it 4 ; but can be sibjected to change.
    % disp(H);
    
    
    % evolve things:
    U = expm(-1j * timestep * time_scaling * discrepancy_factor * H);
    St = U*St;
    
    
    reliability_as_per_qmodel = state_vec * (abs(St).^2);
    
    R(end+1) = reliability_as_per_qmodel;
end


figure;plot(t, R)