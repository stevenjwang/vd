function car = makeGG(paramArr, car)
    
    % Flatten the struct array
    arr = paramArr(:);
    N = numel(arr);
    
    % Pre-allocate cell arrays (This prevents memory reallocation in the loop)
    gg_cell    = cell(N, 1);
    ss_cell    = cell(N, 1);
    accel_cell = cell(N, 1);
    decel_cell = cell(N, 1);
    longA_cell = cell(N, 1);
    longD_cell = cell(N, 1);

    for i = 1:N
        pSet = arr(i);
        
        % Temporary matrices for this iteration
        gg_temp = [];
        
        if pSet.maxLatFlag == 1
            ss_cell{i} = pSet.maxLatx;
            
            % x: long | y: lat | z: velo
            if pSet.maxLongFlag == 1
                p1 = [pSet.maxLongLongAccel, pSet.maxLongLatAccel, pSet.longVel];
                p2 = [pSet.maxLongLongAccel, -pSet.maxLongLatAccel, pSet.longVel];
                
                accel_cell{i} = pSet.maxLongAccelx;
                longA_cell{i} = p1;
                gg_temp = [gg_temp; p1; p2];
            end
            
            if pSet.maxBrakeFlag == 1
                p3 = [pSet.maxBrakeLongDecel, pSet.maxBrakeLatAccel, pSet.longVel];
                p4 = [pSet.maxBrakeLongDecel, -pSet.maxBrakeLatAccel, pSet.longVel];
                
                decel_cell{i} = pSet.maxBrakeDecelx;
                longD_cell{i} = p3;
                gg_temp = [gg_temp; p3; p4];
            end
            
            % Store the combined gg points for this iteration
            gg_cell{i} = gg_temp;
        end
    end

    % vertcat (Vertical Concatenation) is a highly optimized C-backend 
    % function that merges the cells in a single memory operation.
    car.ggPoints        = vertcat(gg_cell{:});
    car.ss_info         = vertcat(ss_cell{:});
    car.accel_info      = vertcat(accel_cell{:});
    car.decel_info      = vertcat(decel_cell{:});
    car.longAccelLookup = vertcat(longA_cell{:});
    car.longDecelLookup = vertcat(longD_cell{:});
end