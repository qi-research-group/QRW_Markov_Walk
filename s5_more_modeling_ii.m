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

dhist = dir('../../history_files/history_*.txt');
drating = dir('../../history_files/ratings_*.txt');

h = load([dhist(1).folder filesep dhist(1).name]);
h = h(h<100);

r = load([drating(1).folder filesep drating(1).name]);
begin_rating = r(2);


this_sub = regexp(dhist(1).name,'[0-9]{3}','match'); this_sub=this_sub{1};



time_scaling = 1/20; % general rotational shenanigans; a scaling factor
% d_discrepancy_factor = 0; % this is how much an oddball will (extra) rotate with the negative H; how impactful is a miss?
Hp_SLOPE = 1; % another scaling factor; of the forces that are embedded within the H.
Hm_SLOPE = 1; 
sigmasquared = 1;



t = [];
fit_params=struct();
fit_params.time_scaling = 1/20;
fit_params.Hp_SLOPE = 1; % the 'potential', the book makes mention of, linear, to POSITIVE states
fit_params.Hm_SLOPE = 1; % the 'potential', the book makes mention of, linear, to NEGATIVE states
fit_params.sigmasquared = 1; % the diffusion off-diagonal

function model_the_behaviour_with_random_walk(h, r, t, fit_params, plot_it)
% h = history (which of the 4 possibles did we get, or none?)
% there are 560 items here. each one, even a missed one, has a certain
% value.
% r = the ratings (22 ratings in total; 1 begin, 1 after prac, and 20 after
% each block.
% t = the times (in seconds) of all the events encoded in h. For now we
% will just leave it empty.
% fit_params (this is a struct with stuff in it that contain the fitting
% parameter of the Hamiltonian (and other parameters) that we will try to
% fit to the data.

% this function assumes:
% - 1 H+ for events in which the AI MATCHES
% - 1 H- for events in which the AI MISMATCHES
% - sigmasquared - a measure of the diffusivity of the system - or how
% likely is someone to change their mind -modeled with the Q R Walk
% - time_scaling; this is a number with which the time is multipled, i.e.
% how fast shall we rotate?

% if you only have 1 unchanging H, then things will always be PERIODIC -
% you can never drift towards anything lower than the initial rating.
% with 2 H, it allows to do non-commutative rotations. Then things become
% more interesting.



time_scaling = fit_params.time_scaling;
Hp_SLOPE = fit_params.Hp_SLOPE; % the 'potential', the book makes mention of, linear, to POSITIVE states
Hm_SLOPE = fit_params.Hm_SLOPE; % the 'potential', the book makes mention of, linear, to NEGATIVE states
sigmasquared = fit_params.sigmasquared;







% base level of trust: this will just be the offset; We will take the
% rating % AFTER PRACTICE BLOCK - and just subtract all ratings with that
% so the results are around 0.
% then the H might do the actual trick.


% let's make S0 the start rating; and let's put that in the middle pls:
% #sum(np.square(S0)) # sum of all elements of S0 squared has to be equal 1
% well, let's set the S0 to:
% S0 = zeros(1, nstates);

% rating after practice block is our initial state vector.
which_one = ceil(r(2)/10);
S0 = zeros(10,1);

S0(which_one) = 1;
if which_one > 1
    S0(which_one-1) = 0.5;
end
if which_one < size(S0,1)
    S0(which_one+1) = 0.5;
end

S0 = S0/sqrt(S0'*conj(S0));
S0 = reshape(S0, numel(S0), 1);


% let's begin the simulation.
St = S0;
R = [0];

timestep = 3.5; % (in our seconds

if numel(t) == 0
    t = [0 (1:numel(h))*timestep]; % time = 0 at the beginning;
end

% delta_t = 4;

% parameter to be optimized


rating_states = [0.5:9.5];




% mv = 1;
% delta_increment = 1; %  0.1;
% delta_decrement = 1;  %0.2;

rating_values = [];
for i=1:numel(h) % i.e. the history:
    

    % get the hamiltonian based on match ve mismatch
    switch h(i)
        case {4, 6} % which is HRAIR, and HFAIF
            
            diag_H = [-4.5:4.5];
            H_SLOPE = Hp_SLOPE;
            
        case {5, 7} % which is HRAIF, and HFAIR
            
            % depending on match or mismatch; changeing the slope(s):
            diag_H = [4.5:-1:-4.5];
            H_SLOPE = Hm_SLOPE;

        otherwise
            
            % just diffuse it for a while more:
            diag_H = [4.5:-1:-4.5] * 0;
            H_SLOPE = 0;
            
    end
    
    H = make_H_v2(H_SLOPE * diag_H, sigmasquared, 0);
    
    timestep = t(i+1) - t(i);
    % t(end+1) = t(end) + timestep; % for artuments sake; make it 4 ; but can be sibjected to change.
    % disp(H);
    
    
    % evolve things:
    U = expm(-1j * timestep * time_scaling * H);
    St = U*St;
    
    
    reliability_as_per_qmodel = rating_states * (abs(St).^2);
    
    R(end+1) = reliability_as_per_qmodel;
    
    if rem(i, 28) ==0
        rating_values(end+1,:) = [i R(end)];
    end
    
end


fh=figure('color','w');
ah=axes;
set(ah,'nextplot','add');
plot(t, R,'k-');
plot(0,  r(2)/10, 'ks');
for i=1:size(rating_values,1)
    plot(t(rating_values(i))*[1 1], rating_values(i,2),'rs');
end

title(sprintf('reliability slider modeling, person: %s\nmodeled with 10 states', this_sub));
xlabel('time (s)');
ylabel('reliability');
set(ah,'yticklabel', [0:10:100],'ytick',[0:10]);