function filtIn = temporalRF_LC10(inRField,times,tcValue)
%--------------------------------------------------------------------------
%   Temporal component of LC10a receptive field for spiking model (v1.0).
%--------------------------------------------------------------------------
%
%   This function computes the input current to LC10a neurons in a
%   spiking network model of courtship behavior, based on the
%   receptive fields estimated by Ribeiro et al. (2018). 
%
%   Inputs: 
%       inRField: NxT binary matrix indicating whether the visual 
%                 stimulus was in the receptive field of each model
%                 neuron in each frame (N = number of model
%                 neurons, T = number of past frames)
%                 
%       times: vector of time stamps for each past frame
%
%       tcValue: ?F/F0 of P1 neurons at the closest imaging frame. 
%                This value is set to unity if the time series is not
%                incorporated into the model
%
%
%   Outputs: 
%       filtIn: total input current to each LC10a neuron in amperes. 
%
%--------------------------------------------------------------------------
%   Paper: An arousal-gated visual circuit controls pursuit during
%          Drosophila courtship (Nature, 2021)
%   Author: Tom Hindmarsh Sten, 2021
%   e-mail: thindmarsh@rockefeller.edu
%   Version: 1.0
%--------------------------------------------------------------------------
%
%% Compute input

%find number of LC10a neurons in model
nNeurons = size(inRField,1); 

%set receptive-field parameters;
kappa = 0.84962;
sigma = 5.5273;
alpha = -0.1859;
beta = 15.5884;
Sf = 25e-10;

%compute potential total input current to each LC10a neuron
sumCurve = exp(sigma.*(kappa + times))...
           ./((exp(beta.*(times-alpha))+1)...
           .*(exp(sigma.*(kappa + times))+1));

%sum input current when stimulus was in spatial RF of each neuron
%and convert to amps. 
filtIn = sum(inRField.*repmat(sumCurve',nNeurons,1),2).*Sf.*tcValue;

end