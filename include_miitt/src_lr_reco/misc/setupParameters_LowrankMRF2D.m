function params = setupParameters_LowrankMRF2D( )

params = struct();
params.updateFreq =         1; % how often to display intermediate images
params.lambdaLLR =          0.02; % locally low-rank image patch denoising
% if N==192
params.block_dim =          [6,6]; % the patch size should divide evenly into the matrix size (e.g., 192 is divisible by 6)
params.block_step =         6;
% elseif N==256
%     params.block_dim =          [8,8]; % the patch size should divide evenly into the matrix size (e.g., 192 is divisible by 6)
%     params.block_step =         8;
% end
params.lambdaWav =          0; % wavelet denoising
params.lambdaSpatialTV=     5e-4; % total variation denoising (along x and y dimensions)
params.lambdaTemporalTV=    0; % total variation denoising (MRF contrast dimension)
params.lambdaTemporalFT=    0; % total variation denoising (time dimension)
params.numIter =            100; % total number of iterations
params.beta =               1;
params.stopThresh =         0; % stopping criterion
params.t0 =                 1;
params.alpha =              1.2; % adjusts the OGM1 step size