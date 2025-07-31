function [out,varargout] = model_the_behaviour_with_random_walk(h, r, t, fit_params, SD, plot_it, this_sub, COLLAPSE)
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





time_scaling = fit_params(1);
Hp_SLOPE = fit_params(2); % the 'potential', the book makes mention of, linear, to POSITIVE states
Hm_SLOPE = fit_params(3); % the 'potential', the book makes mention of, linear, to NEGATIVE states
sigmasquared = fit_params(4);% .sigmasquared;







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
% S0 = zeros(10,1);
% 
% S0(which_one) = 1;
% if which_one > 1
%     S0(which_one-1) = 0.5;
% end
% if which_one < size(S0,1)
%     S0(which_one+1) = 0.5;
% end
% 
% S0 = S0/sqrt(S0'*conj(S0));
% S0 = reshape(S0, numel(S0), 1);
% my_SD = 0.05;
my_SD = SD;

% well, let's re-do this!
rating_states = [0.5:9.5];

f = @(x, m, s) 1/sqrt(2*pi*s^2) * exp (-1*(x-m)^2 / (2*s^2));
S0_vec = @(m, s) arrayfun(@(x) f(x, m, s), rating_states)';
S0_vec_n = @(m, s) S0_vec(m, s) / sqrt(S0_vec(m, s)' * S0_vec(m, s));

S0 = S0_vec_n(r(2)/10, my_SD);
    






% let's begin the simulation.
St = S0;
R = [rating_states * (abs(St).^2)];

timestep = 3.5; % (in our seconds

if numel(t) == 0
    t = [0 (1:numel(h))*timestep]; % time = 0 at the beginning;
end

% delta_t = 4;

% parameter to be optimized
% keyboard;




% mv = 1;
% delta_increment = 1; %  0.1;
% delta_decrement = 1;  %0.2;

delta_t_values = [];
state_values = [];
rating_values = [];
for_plot = [];
for_plot_t = [];
for_plot_t = 0; % - 3.5;

for i=1:numel(h) % i.e. the history:
    
    for_plot_t(end+1) = t(i) - t(1) + 3.5;
    
    % get the hamiltonian based on match ve mismatch
    switch h(i)
        case {4, 6} % which is HRAIR, and HFAIF
            
            diag_H = [-4.5:4.5];
            H_SLOPE = Hp_SLOPE;
            for_plot(end+1) = 1;
            
        case {5, 7} % which is HRAIF, and HFAIR
            
            % depending on match or mismatch; changeing the slope(s):
            diag_H = [4.5:-1:-4.5];
            H_SLOPE = Hm_SLOPE;
            for_plot(end+1) = 2;
            
        otherwise
            
            % just diffuse it for a while more:
            % diag_H = [4.5:-1:-4.5] * 0;
            diag_H = ones(1,10);
            H_SLOPE = 0;
            for_plot(end+1) = 3;
            
    end
    
    H = make_H_v2(H_SLOPE * diag_H, sigmasquared, 0);
    
    % use standard length of the time pls...
    % we essentially 'pause' the process whenever not in a block.
    % hm, we can make things even more beautiful. but for now this is
    % still pretty good.
    if ismember(i, 1:28:560) 
        if ismember(h(i), 4:7)
            % if press is any of HRAIR, HRAIF, HFAIFHFAIR, then:
            timestep = median(diff(t));
        elseif h(i) == 3
            % so if it's pressed-too-late (no response), 2.5 seconds...
            timestep = 2.5;
        end
    else
        % keyboard;
        timestep = t(i) - t(i-1);
    end
    
    % t(end+1) = t(end) + timestep; % for artuments sake; make it 4 ; but can be subjected to change.
    % disp(H);
    
    delta_t_values(end+1) = timestep;
    
    % evolve things:
    U = expm(-1j * timestep * time_scaling * H);
    St = U*St;
    
    
    reliability_as_per_qmodel = rating_states * (abs(St).^2);
    
    R(end+1) = reliability_as_per_qmodel;
    
    if rem(i, 28) ==0
        rating_values(end+1,:) = [i R(end)];
        
        if COLLAPSE
            
            % keyboard;
            
            % NOW!!! we just need the RESPONDED R; i.e., the collapsed R
            % 
            pick_this_one = floor(i/28);
            measured_value = r(pick_this_one+2)/10;
            
            new_amps_to_use = S0_vec_n(measured_value, my_SD);
            my_angles = angle(St);
            
            new_SC = (cos(my_angles) + 1j*sin(my_angles)).*new_amps_to_use;
            
            % keyboard;
            
            % SC = St .* new_amps_to_use;
            % SC2 = SC/sqrt(SC'*SC);
            SC2 = new_SC;
            
            % keyboard;
            


            %             which_one = ceil(r(pick_this_one+2)/10);
            %             
            %             % which_one = ceil(R(end));
            %             
            %             % COLLAPSE THINGS
            %             % COLLAPSE STRATEGY = LIke setting the initial S0.
            %             % but, we can collapse with different strategies...
            %             
            %             % SC = St;
            % 
            %             SC = zeros(size(St));
            %            
            %             try
            %             SC(which_one) = St(which_one) * 1;
            %             if which_one > 1
            %                 SC(which_one-1) = St(which_one-1) * 0.5;
            %             end
            %             if which_one < size(SC,1)
            %                 SC(which_one+1) = St(which_one+1)* 0.5;
            %             end
            %             catch
            %                 keyboard;
            %             end
            %             
            %             SC2 = SC/sqrt(SC'*SC);
            %             SC=SC2;
            %             % SC = SC/sqrt(SC'*conj(SC));
            %             SC = reshape(SC, numel(SC), 1);
            %             
            %             
            %             % keyboard;
            %             % setting the state:
            %             
            %             % fprintf('Collapsed!');
            
            
            
            % keyboard;
            
            
            St = SC2;
            
            try
                R(end+1) = rating_states * (abs(St).^2);
                for_plot_t = [for_plot_t t(i+1) - t(1)];
            catch
                % well - we're not going to actually use this timestep,
                % as this is the S0 for block 21...
                for_plot_t = [for_plot_t for_plot_t(end) + timestep];
            end
            % disp(St);
            
        end
        
    end
    
    state_values(:, end+1) = St;
    
end


out = rating_values(:,2);

 %  keyboard;

fh=[];
l=@(x) (1:29)+(x-1)*29;
is_to_plot = arrayfun(l, 1:20,'uniformoutput', false);
if plot_it
    fh=figure('color','w');
    ah=axes;
    set(ah,'nextplot','add');
    
    try
    for ibl=1:20
        % keyboard;
        % keyboard;
        plot(for_plot_t(is_to_plot{ibl}), R(is_to_plot{ibl}),'k-','linewidth',2);
        plot(for_plot_t(is_to_plot{ibl}(1)), R(is_to_plot{ibl}(1)),'o','linewidth',1,'markerfacecolor',[0.4 0.4 0.4], 'markeredgecolor',[0.4 0.4 0.4]);
        
    end
    catch
        keyboard;
    end
    % initial state:
    plot(0,  r(2)/10, 'ks','markerfacecolor','k');
    for i=1:size(rating_values,1)
        try
        plot(t(rating_values(i)) - t(1)+3.5, rating_values(i,2),'rs','markerfacecolor','r');
        plot(t(rating_values(i)) - t(1)+3.5, r(i+2)/10,'bo', 'markerfacecolor','b');
        
        plot(t(rating_values(i))*[1 1] - t(1)+3.5, [rating_values(i,2) r(i+2)/10],'linestyle','-','color',[0.7 0.7 0.7]);
        
        catch
            keyboard;
        end
    end
    % keyboard;
    my_x = [0 (t(rating_values(:,1)) - t(1))];
    my_y = r(2:end)/10;
    plot(my_x+3.5, my_y, 'b-');
    % keyboard;
    
    % plot(t(rating_values(:,1)) - t(1), rating_values(:,2),'r-');
    
    
    
    
    
    % keyboard;
    my_colors = {'b','r','k'};
    % figure out the length of t... let's take some reasonablish value...
    reasonable_t_length_value = prctile(diff(t), 25);
    for i=1:numel(for_plot)
        start_t_value = t(i) - t(1)+1.5;
            
        ph=patch( [0 1 1 0] * reasonable_t_length_value + start_t_value, [0 0 1 1]/2, my_colors{for_plot(i)});
        set(ph,'linestyle','none');
        
        
    end
    
    title_to_use = sprintf('reliability slider modeling, person: %s\nmodeled with 10 states', this_sub);
    for i_fit_params = 1:numel(fit_params)
        this_fit_param = fit_params(i_fit_params);
        title_to_use = sprintf('%s\nparam%d: %2.4f', title_to_use, i_fit_params, this_fit_param);
    end
    
    % fit_params
    
    title(title_to_use);
    xlabel('time (s)');
    ylabel('reliability');
    set(ah,'yticklabel', [0:10:100],'ytick',[0:10]);
    set(ah,'xlim', [0 t(end) - t(1) + 3.6]);
end

nout = max(nargout,1) - 1;
if nout > 0
    varargout{1} = fh;
end
if nout > 1
    varargout{2} = delta_t_values;
end
if nout > 2
    varargout{3} = state_values;
end