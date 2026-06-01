% Lapsim2
% Main steady-state lapsim script. This is traditional lapsim. It finds the
% g-g diagram over different velocities (max lateral and longitudinal
% acceleration at a given velocity), then uses those to predict performance
% in the dynamic events.
% HOW TO USE:
% 1) carConfig.m: define the car you want to test
% 2) run this script2
% 3) results are in the comp object, which is stored in the corresponding
% car object. Index into carCell to get the car you want, then open the
% comp object.
% Lapsim2 (Optimized DOE Engine)
% Main steady-state lapsim script optimized for large-scale Latin Hypercube Sampling.

clear
setup_paths

% 1. Generate all cars to sim over using the LHS Parameters Loop
carCell = carConfig(); 
numCars = size(carCell,1);

% 2. Setup Parallel Pool
% Let MATLAB handle the workers dynamically. Hardcoding numWorkers = 16 
% can crash if a teammate runs this on an 8-core laptop.
poolobj = gcp('nocreate'); % Check if pool exists
if isempty(poolobj)
    disp('Starting parallel pool. This takes a moment...');
    parpool(); 
end

fprintf('Starting DOE with %d car setups.\n', numCars);

% 1. Pre-allocate TWO separate 1D cell arrays
carOut_main  = cell(numCars, 1); 
carOut_accel = cell(numCars, 1);

tic;
parfor i = 1:numCars
    % Pull from carCell (broadcast/sliced automatically)
    car = carCell{i, 1};
    accelCar = carCell{i, 2};
    
    % --- Step 1 & 2: Physics & Envelope ---
    raw_gg_params = gg2(car, 0); 
    car = makeGG(raw_gg_params, car); 
    
    % --- Step 3: Integration ---
    comp = EventsDOE(car, accelCar);
    comp.calcTimes();       
    
    % --- Step 4: Storage in separate "buckets" ---
    car.comp = comp;        
    
    % Use single subscripts {i} for both to keep the analyzer happy
    carOut_main{i}  = car;      
    carOut_accel{i} = accelCar; 
    
    fprintf('Completed Setup %d of %d\n', i, numCars);
end

% 2. Zip them back together into the 16x2 structure the RSM needs
carCell = [carOut_main, carOut_accel]; 

% =========================================================================

total_time = toc;
fprintf('\n--- DOE COMPLETE ---\n');
fprintf('Total time elapsed: %.1f seconds (%.2f minutes)\n', total_time, total_time/60);
fprintf('Average time per setup: %.2f seconds\n', total_time/numCars);


