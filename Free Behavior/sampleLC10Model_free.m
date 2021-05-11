%--------------------------------------------------------------------------
%   Spiking model of LC10a circuit for free courtship (v1.0)
%--------------------------------------------------------------------------
%
%   This script simulates spiking of LC10a neurons, modeled as independent
%   integrate-and-fire units, and transforms spiking into behavioral
%   output. The input to the model is the estimated angular position of a 
%   female conspecific fly relative to the male over the course of
%   courtship. 
%
%   The visual field is divvied up between the total number of LC10a
%   neurons, with 15 degrees of binocular overlap. LC10a neurons only
%   receive input current for times when visual motion occured in their
%   spatial receptive fields. At the end of the simulation, total spiking
%   across the LC10a population in the left versus right hemisphere is
%   compared in each 30ms time-bin. Turning is proportional to this
%   difference. 
%
%   The script can be set to run with the assumption of no motion-direction
%   selectivity ('none'), selectivity for progressive motion ('progressive'),
%   or selectivity for regressive motion ('regressive');
%
%   Sample behavioral, neural, and stimulus data can be found in
%   'Dmel_freeBehavior_1fem1male_pair1_sampleData'.mat'. Please note that
%   the Circular Statistics Toolbox by Philipp Berens must be installed. 
%   (https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-
%   statistics-toolbox-directional-statistics)
%
%--------------------------------------------------------------------------
%   Paper: Sexual arousal gates visual processing during 
%          Drosophila Courtship (Nature, 2021)
%   Author: Tom Hindmarsh Sten, 2021
%   e-mail: thindmarsh@rockefeller.edu
%   Version: 1.0
%--------------------------------------------------------------------------
%% Set parameters and load behavioral data

%load data
load('Dmel_freeBehavior_1fem1male_pair1_sampleData')

%set direction selectivity: 'progressive','regressive', or 'none'
DIR_SELECTIVITY = 'progressive';

%set number of LC10a neurons in model ? even numbers (nLC/2 per hemisphere).
nLC = 40;

%set the degree of binocular overlap
binocRegionSize = pi/12;

%set the hemifield sizes
hemifieldSizes = [-pi/4 pi/2+binocRegionSize;pi/2-binocRegionSize 5*pi/4];

%divvy up the right hemifield equally between the LC10a neurons
rightRFBounds = linspace(hemifieldSizes(1,1),hemifieldSizes(1,2),nLC/2+1)';
rightRFBounds(1:end-1,2) = rightRFBounds(2:end,1);
rightRFBounds(end,:) = [];

%divvy up the left hemifield equally between the LC10a neurons
leftRFBounds = linspace(hemifieldSizes(2,1),hemifieldSizes(2,2),nLC/2+1)';
leftRFBounds(1:end-1,2) = leftRFBounds(2:end,1);
leftRFBounds(end,:) = [];

%concatenate left and right hemisphere neurons into one vector
spatialRFs = cat(1,rightRFBounds,leftRFBounds);

%Set model neuron parameters
tau = 0.01;                                  % Time-constant (s)
vRest = -65e-3;                              % Resting potential (V)
vReset = -65e-3;                             % Reset potential (V)
vThresh = -50e-3;                            % Spike threshold (V)
R_m = 10e6;                                  % Total input resistance (Ohm)
r_m = 10e6;                                  % Leaky K channel resistance (Ohm)
deltaGsra = 14e-9;                           % change in adaptation conductance (g)
tau_sra = 0.2;                               % time constant of SR adaptation (s)
E_K = -70e-3;                                % Potassium reversal potential (V)

spikeMag = 40e-3;                            % Magnitude of evoked spikes (V)

%Set integration properties
t = max(bTime);                              % Model run time (s)
dt = 0.003;                                  % Step-size for Euler's (s)

%Initialize vectors
N = floor(t/dt);                             % Number of time-steps
V_exc = zeros(N,nLC);                        % Voltage vector
V_exc(1,:) = vRest;                          % Set neurons to resting potential
Gsra_exc = zeros(N,nLC);                     % Spike-rate adaptation currents

%% Run model
  
for i = 1:N-1
    
    %Track and print progress
    if rem(dt*i,3) == 0
        clc
        fprintf('%ds have passed out of %ds total',dt*i,t);
    end
    
    %find current time and index all time-points up until it.
    currTime = dt*i;
    hist = find(bTime <= currTime);
    corrTime = bTime(hist) - currTime;
    
    %check which past time points stimulus was moving in the spatial
    %receptive field of each neuron
    for k = 1:nLC

        %get distance of female from upper and lower bound of RFs
        d1 = abs(circ_dist(femaleAngle(hist,:),spatialRFs(k,1)));
        d2 = abs(circ_dist(femaleAngle(hist,:),spatialRFs(k,2)));
        
        %get width of RF
        td = abs(circ_dist(spatialRFs(k,1),spatialRFs(k,2)));
        
        %check whether female was within spatial RF for each frame
        presence = (d1 < td & d2 < td)';

        %compute motion direction of target
        motionDirection = sign(diff(femaleAngle(hist)));
        
        %constrict to motion direction depending on setting
        if strcmpi(DIR_SELECTIVITY,'progressive')
            
            %if neuron is on right side, include when stim moves right
            if k <= nLC/2
                inField(k,:) = presence.*[0 motionDirection'<0];
            %%if neuron is on left side, include when stim moves left
            else
                inField(k,:) = presence.*[0 motionDirection'>0];
            end
        
        %flip contingencies if regressive selective
        elseif strcmpi(DIR_SELECTIVITY,'regressive')
            
            %if neuron is on right side, include when stim moves left
            if k <= nLC/2
                inField(k,:) = presence.*[0 motionDirection'>0];
            %%if neuron is on left side, include when stim moves right
            else
                inField(k,:) = presence.*[0 motionDirection'<0];
            end
            
        %if no selectivity, just use spatial position
        elseif  strcmpi(DIR_SELECTIVITY,'none')
            inField(k,:) = presence;
        end
    end
    
    %get input current to each LC10a neuron 
    filtInput = temporalRF_LC10(inField,corrTime,1);
    
    % compute voltage and record spiking in each model neuron
    for j = 1:nLC
        
        %if last time point was a spike, reset voltage to baseline
        if V_exc(i,j) == spikeMag                          
            V_exc(i+1,j) = vReset;
        
        else

           exc = R_m*filtInput(j); %input
           adpt = r_m*Gsra_exc(i,j)*(V_exc(i,j) - E_K); %leak 
           fV = (vRest - V_exc(i,j)) + exc - adpt; %forward voltage
           
           V_exc(i+1,j) = V_exc(i,j) + (dt/tau)*fV; %new voltage
            
            %if threshold is exceeded, ~~ SPIKE! ~~
            if V_exc(i+1,j) > vThresh
                V_exc(i+1,j) = spikeMag;
            end
        end
        
        %modify adaptation current: decay if neuron didn't spike, add if it
        %did. 
        if V_exc(i,j) == spikeMag
            Gsra_exc(i+1,j) = Gsra_exc(i,j) + deltaGsra;
        else
            Gsra_exc(i+1,j) = Gsra_exc(i,j) + (dt/tau_sra)*(-Gsra_exc(i,j));
        end
        
    end
    
    clear inField

end
%% Transform spiking to turning ? compare to behavior

binSize = 30; %ms
modelTimeStamp = [0:dt:t-dt];
deltaHeading = zeros(length(modelTimeStamp)-binSize,1);
for i = 1:length(modelTimeStamp)-binSize
    %compute spikes in right and left hemisphere; compare
    rightSpikes = sum(sum(V_exc(i:i+binSize,1:nLC/2) == spikeMag));
    leftSpikes = sum(sum(V_exc(i:i+binSize,nLC/2+1:nLC) == spikeMag));
    
    deltaHeading(i) = leftSpikes - rightSpikes;
end

plot(bTime(2:end),hDiff,'g'); yyaxis right; 
plot(modelTimeStamp(1:end-binSize),deltaHeading,'b'); 
legend('Animal Turning (rad/s)','Model Turning (spikes)')


