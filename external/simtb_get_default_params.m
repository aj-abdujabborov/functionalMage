function  P = simtb_get_default_params(sourceType)
%% This is a SimTB function.

if sourceType == 1
    %% 'Convolution with canonical HRF';
    % HRF parameters vary slightly between components/subjects
    % P: 1x7 vector of parameters for the response kernel (difference of two gamma functions)
    P(1) = 4+3*abs(randn(1));     % delay of response (relative to onset)
    P(2) = 12+3*abs(randn(1));    % delay of undershoot (relative to onset)
    P(3) = 1+0.2*randn(1);        % dispersion of response
    P(4) = 1+0.2*randn(1);        % dispersion of undershoot
    P(5) = 2+3*abs(randn(1));     % ratio of response to undershoot
    P(6) = 0;                     % onset (seconds)
    P(7) = 32;                    % length of kernel (seconds)

elseif sourceType == 2
    %% 'Hemodynamic nonlinear model (Windkessel Balloon Model)';
    % model parameters vary slightly between subjects/components.
    % P: 1x7 vector of parameters; see Figure 7 of Friston et al., 2000.
    P(1) = 1/(1.54+0.15*randn(1));           % signal decay
    P(2) = 1/(2.46+0.15*randn(1));           % autoregulation
    %P(3) = 2+0.15*randn(1);                  % transit time
    P(3) = 0.98 + 0.15*randn(1);                  % transit time
    P(4) = 0.36;                             % stiffness
    P(5) = 0.34;                             % OEF
    P(6) = 0.04;                             % TE
    P(7) = 0.54;                             % neural efficacy
elseif sourceType == 3
    %% 'Convolution with fast spike';
    % 'Spike' parameters vary slightly between components/subjects
    % P: 1x7 vector of parameters for the response kernel (difference of two gamma functions)
    P(1) = 2+.05*randn(1);        % delay of response (relative to onset)
    P(2) = 6+0.05*randn(1);       % delay of undershoot (relative to onset)
    P(3) = 0.8+0.02*randn(1);     % dispersion of response
    P(4) = 1+0.02*randn(1);       % dispersion of undershoot
    P(5) = 4;                     % ratio of response to undershoot
    P(6) = 0;                     % onset (seconds)
    P(7) = 20;                    % length of kernel (seconds)
else
    P = [];
end