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
dtimings = dir('../../history_files/times_*.txt');

% location_of_raw_datafiles = '/home/vanderj2/mnt/lyrascratch/faces_experiment/preprocessed';

SDs = [0.75]; %[0.05:0.10:2.0];
saved_model_parameters = {};
for i_d = 1:numel(dhist)
    
    for i_SD = 1:numel(SDs)
        this_SD = SDs(i_SD);
        
        %         % initial state:
        %         plot(0,  r(2)/10, 'ks','markerfacecolor','k');
        %         for i=1:size(rating_values,1)
        %             try
        %                 plot(t(rating_values(i)) - t(1)+3.5, rating_values(i,2),'rs','markerfacecolor','r');
        %                 plot(t(rating_values(i)) - t(1)+3.5, r(i+2)/10,'bo', 'markerfacecolor','b');
        %                 
        %                 plot(t(rating_values(i))*[1 1] - t(1)+3.5, [rating_values(i,2) r(i+2)/10],'linestyle','-','color',[0.7 0.7 0.7]);
        %                 
        %             catch
        %                 keyboard;
        %             end
        %         end
        %         % keyboard;
        %         my_x = [0 (t(rating_values(:,1)) - t(1))];
        %         my_y = r(2:end)/10;
        %         plot(my_x+3.5, my_y, 'b-');
        %         % keyboard;
        
        % plot(t(rating_values(:,1)) - t(1), rating_values(:,2),'r-');
        
        
        
        h = load([dhist(i_d).folder filesep dhist(i_d).name]);
        h = h(h<100); % omit history of practice block
        
        r = load([drating(i_d).folder filesep drating(i_d).name]);
        begin_rating = r(2); % omit rating before practice block
        
        t = load([dtimings(i_d).folder filesep dtimings(i_d).name]);
        t(1:28) = []; % omit practice block
        
        
        this_sub = regexp(dhist(i_d).name,'[0-9]{3}','match'); this_sub=this_sub{1};
        
        % wait a minute; the same routine that harvests the h can also grab the
        % times, no?
        %     EEGS=[];
        %     for i_eeg_f = 1:3
        %         eeg_setfile = sprintf('%s/sub-%s/s1_%s_block_%d_ready_for_inspection.set', location_of_raw_datafiles, this_sub, this_sub, i_eeg_f);
        %         if ~exist(eeg_setfile,'file')
        %             fprintf('file not found!!\n')
        %         else
        %             EEGS = [EEGS pop_loadset('filepath', sprintf('%s/sub-%s', location_of_raw_datafiles, this_sub), 'filename', sprintf('s1_%s_block_%d_ready_for_inspection.set', this_sub, i_eeg_f))];
        %         end
        %     end
        
        
        
        time_scaling = 1/20; % general rotational shenanigans; a scaling factor
        % d_discrepancy_factor = 0; % this is how much an oddball will (extra) rotate with the negative H; how impactful is a miss?
        Hp_SLOPE = 1; % another scaling factor; of the forces that are embedded within the H.
        Hm_SLOPE = 1;
        sigmasquared = 1;
        SD = this_SD; %0.05;
        
        
        
        % let's try to run 1 model, then:
        
        % t = [];
        fit_params=struct();
        fit_params.time_scaling = 1/20;
        fit_params.Hp_SLOPE = 1; % the 'potential', the book makes mention of, linear, to POSITIVE states
        fit_params.Hm_SLOPE = 1; % the 'potential', the book makes mention of, linear, to NEGATIVE states
        fit_params.sigmasquared = 1; % the diffusion off-diagonal
        fit_params.SD = SD;
        
        COLLAPSE = 1;
        
        
    %         
    %         %
    %         %
    %         % let the modeling begin.
    %         %  % initial state:
    %         plot(0,  r(2)/10, 'ks','markerfacecolor','k');
    %         for i=1:size(rating_values,1)
    %             try
    %                 plot(t(rating_values(i)) - t(1)+3.5, rating_values(i,2),'rs','markerfacecolor','r');
    %                 plot(t(rating_values(i)) - t(1)+3.5, r(i+2)/10,'bo', 'markerfacecolor','b');
    %                 
    %                 plot(t(rating_values(i))*[1 1] - t(1)+3.5, [rating_values(i,2) r(i+2)/10],'linestyle','-','color',[0.7 0.7 0.7]);
    %                 
    %             catch
    %                 keyboard;
    %             end
    %         end
    %         % keyboard;
    %         my_x = [0 (t(rating_values(:,1)) - t(1))];
    %         my_y = r(2:end)/10;
    %         plot(my_x+3.5, my_y, 'b-');
    %         % keyboard;
    %         
        % plot(t(rating_values(:,1)) - t(1), rating_values(:,2),'r-');
        
        
        
        %
        %
        
        
        
        
        X0 = [fit_params.time_scaling, fit_params.Hp_SLOPE, fit_params.Hm_SLOPE, fit_params.sigmasquared];
        test_model = model_the_behaviour_with_random_walk(h, r, t, X0 , this_SD, 1, this_sub, COLLAPSE);
        
        
        % out = model_the_behaviour_with_random_walk(h, r, [], X0 , 0, '');
        
        cost_mse = @(x) sum(abs(r(3:end)'/10 - model_the_behaviour_with_random_walk(h, r, t, x , this_SD, 0, '', COLLAPSE))/20); % MSE cost function.
        cost_esq = @(x) sum((r(3:end)'/10 - model_the_behaviour_with_random_walk(h, r, t, x , this_SD, 0, '', COLLAPSE)).^2/20); % MSE cost function.
        
        % keyboard;
        
        cost_fs = {cost_mse cost_esq};
        cost_ls = {'mse','esq'};
        
        for i_cost = 1 %:2
            cost_f = cost_fs{i_cost};
            cost_l = cost_ls{i_cost};
            
            
            do_gradient_descent = 0;
            if do_gradient_descent
                X_optimal_1 = fminsearch(cost_f, X0);
                [out1, fh1] = model_the_behaviour_with_random_walk(h, r, t, X_optimal_1 , this_SD, 1, this_sub, COLLAPSE);
            end
            
            % OK then, let's try a particle swarm then...
            
            % Define the cost function as a sum of squares
            % costFunction = @(x) sum(x.^2);
            
            
            % Set the number of variables
            nvars = 4;  % for a 4-by-1 vector
            
            % Define boundary conditions
            
            lb = [1/100 0.01   0.01   0.01  ];
            ub = [1       10     10     10  ];
            
            do_particle_swarm = 0;
            if do_particle_swarm
                % Set options for particleswarm (optional settings for display and max iterations)
                options = optimoptions('particleswarm', 'Display', 'iter', 'MaxIterations', 2000);
                
                % Run the particleswarm optimization
                [x_optimal_2, fval] = particleswarm(cost_f, nvars, lb, ub, options);
                
                % Display the results
                disp('Optimal vector:');
                disp(x_optimal_2);
                disp('Minimum sum of squares:');
                disp(fval);
                
                [out2, fh2] = model_the_behaviour_with_random_walk(h, r, t, x_optimal_2 , this_SD, 1, this_sub, COLLAPSE);
            end
            
            do_globalsearch = 1;
            if do_globalsearch
                % Set up the problem for fmincon with boundaries
                problem = createOptimProblem('fmincon', 'objective', cost_f, ...
                    'x0', X0, 'lb', lb, 'ub', ub);
                
                % Set up GlobalSearch
                gs = GlobalSearch;
                
                % Run the optimization using GlobalSearch
                try
                [x_optimal_3, fval] = run(gs, problem);
                catch
                    keyboard;
                end
                
                % Display the results
                disp('Optimal vector:');
                disp(x_optimal_3);
                disp('Minimum sum of squares:');
                disp(fval);
                [out3, fh3] = model_the_behaviour_with_random_walk(h, r, t, x_optimal_3 , this_SD, 1, this_sub, COLLAPSE);
                
                saved_model_parameters{end+1} = {x_optimal_3, fval, out3, fh3, drating(i_d).name, this_SD};
            end
            
            
            fh1=0;
            fh2=0;
            fhs = [fh1, fh2, fh3];
            prefixes = {'gdescent','particle','global'};
            for j=3
                fn_to_use = sprintf('figures_lowSD/model_cost_%s_%s_%s.jpg',prefixes{j}, cost_l, this_sub);
                set(fhs(j),'paperunits','centimeters');
                set(fhs(j),'papersize',[35 30]);
                set(fhs(j),'paperposition',[0 0 35 30]);
                pause(2);
                print('-djpeg','-r300',fn_to_use);
                print('-dsvg',[fn_to_use(1:end-3) 'svg']);
                pause(3);
                
            end
            for j=[3];
                close(fhs(j));
            end
            
            
            
            
            
        end
    end
    
end
saved_model_parameters_lowSD = saved_model_parameters;
saved_model_parameters_multiSD = saved_model_parameters;
% save saved_model_parameters_lowSD.mat saved_model_parameters_lowSD
save saved_model_parameters_multiSD.mat saved_model_parameters_multiSD;