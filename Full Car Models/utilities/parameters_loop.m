function [car_cell] = parameters_loop(cP, aP, eP, DTp, Bp, tP, sampling_type, num_samples)
% parameters_loop Generates a cell array of Car objects for DOE.
%
% Inputs:
%   cP, aP, eP, DTp, Bp, tP: Parameter structs from carConfig
%   sampling_type: "FullFactorial", "LHS", or "Random"
%   num_samples: Total runs for LHS/Random (ignored for FullFactorial)
%
% Output:
%   car_cell: Nx2 cell array [Main_Car, Accel_Car]

    % Set defaults
    if nargin < 7, sampling_type = "FullFactorial"; end
    if nargin < 8 && ~strcmpi(sampling_type, "FullFactorial")
        error('num_samples is required for LHS or Random sampling.');
    end

    %% 1. Define the Master Parameter List
    % This is the source of truth. Every variable sampled is listed here.
    % The second column contains the array of levels (FullFactorial) 
    % or the [min, max] bounds (LHS/Random).
    params_list = {
        'mass', cP.mass;                      % 1
        'dr_wt', cP.driver_weight;            % 2
        'acc_dr_wt', cP.accel_driver_weight;  % 3
        'wb', cP.wheelbase;                   % 4
        'wd', cP.weight_dist;                 % 5
        'tw', cP.track_width;                 % 6
        'rad', cP.wheel_radius;               % 7
        'cg', cP.cg_height;                   % 8
        'rchf', cP.roll_center_height_front;  % 9
        'rchr', cP.roll_center_height_rear;   % 10
        'rsf', cP.R_sf;                       % 11
        'izz', cP.I_zz;                       % 12
        'ack', cP.ackermann;                  % 13
        'ccomp', cP.camber_compliance;        % 14
        'cda', aP.cda;                        % 15
        'cla', aP.cla;                        % 16
        'aero_dist', aP.distribution;         % 17
        'acc_cda', aP.accel_cda;              % 18 (Special for Accel)
        'acc_cla', aP.accel_cla;              % 19 (Special for Accel)
        'redline', eP.redline;                % 20
        'shift_pt', eP.shift_point;           % 21
        'shift_t', eP.shift_time;             % 22
        'fd', DTp.final_drive;                % 23
        'eff', DTp.drivetrain_efficiency;     % 24
        'gd1', DTp.G_d1;                      % 25
        'gd2_o', DTp.G_d2_overrun;            % 26
        'gd2_d', DTp.G_d2_driving;            % 27
        'b_dist', Bp.brake_distribution;      % 28
        'b_trq', Bp.max_braking_torque;       % 29
        'gamma', tP.gamma;                    % 30
        'pi', tP.p_i                          % 31
    };

    num_vars = size(params_list, 1);

    %% 2. Sampling Generation Logic
    switch lower(sampling_type)
        case 'fullfactorial'
            grid_args = cell(1, num_vars);
            [grid_args{:}] = ndgrid(params_list{:,2});
            num_runs = numel(grid_args{1});
            final_params = zeros(num_runs, num_vars);
            for v = 1:num_vars
                final_params(:, v) = grid_args{v}(:);
            end
            
        case 'lhs'
            num_runs = num_samples;
            lhs_norm = lhsdesign(num_runs, num_vars, 'criterion', 'maximin');
            final_params = zeros(num_runs, num_vars);
            for v = 1:num_vars
                p_range = params_list{v, 2};
                if isscalar(p_range)
                    final_params(:, v) = p_range;
                else
                    final_params(:, v) = min(p_range) + lhs_norm(:, v) * (max(p_range) - min(p_range));
                end
            end
            
        case 'random'
            num_runs = num_samples;
            rand_norm = rand(num_runs, num_vars);
            final_params = zeros(num_runs, num_vars);
            for v = 1:num_vars
                p_range = params_list{v, 2};
                if isscalar(p_range)
                    final_params(:, v) = p_range;
                else
                    final_params(:, v) = min(p_range) + rand_norm(:, v) * (max(p_range) - min(p_range));
                end
            end
    end

    %% 3. Construction Loop
    car_cell = cell(num_runs, 2);
    
    for i = 1:num_runs
        % Row mapping for clarity
        p = final_params(i, :);

        % Sub-Component Construction
        % Main Aero (Autox/Skid)
        aero = Aero(p(15), p(16), p(17)); 
        % Special Accel Aero (Low Drag)
        accel_aero = Aero(p(18), p(19), p(17)); 

        powertrain = Powertrain(p(20), p(21), eP.gears, eP.primary_reduction, ...
            eP.torque_fn, p(22), p(23), p(7), p(24), p(25), p(26), p(27), p(28), p(29));
            
        tire = Tire2(p(31), tP.Fx_parameters, tP.Fy_parameters, tP.friction_scaling_factor); 

        % Main Car Setup
        car_cell{i, 1} = Car(p(1) + p(2), p(4), p(5), p(6), p(7), p(8), ...
            p(9), p(10), p(11), p(12), p(30), p(14), aero, powertrain, tire, p(13)); 

        % Accel Car Setup (Includes Accel Driver Weight and Accel Aero)
        car_cell{i, 2} = Car(p(1) + p(3), p(4), p(5), p(6), p(7), p(8), ...
            p(9), p(10), p(11), p(12), p(30), p(14), accel_aero, powertrain, tire, p(13)); 
    end
end