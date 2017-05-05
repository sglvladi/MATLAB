function outIndex=resampling(w,N)

% %INPUT:  
% % w:         the weights with sum(w)=1
% % N:         number of resamples you want to resample
% %OUTPUT:
% % outIndex:  indices for resampled particles
% 
if nargin == 1
   N = length(w);
end
% here you can use any unbiased resampling scheme that output
% equally-weighting paticles. Ref:

% T. Li, M. Bolic, P. Djuric, Resampling methods for particle filtering, 
% IEEE Signal Processing Magazine,vol.32, no.3, pp.70-86

% See also:
% T. Li, G. Villarrubia, S. Sun, J. M. Corchado, J. Bajo. 
% Resampling methods for particle filtering: identical distribution, 
% a new method and comparable study, Frontiers of Information Technology 
% & Electronic Engineering (invited), vol.16, no.11, pp.969-984

% Website: https://sites.google.com/site/tianchengli85/

%% this following scheme uses the systematic resampling (faster) 
M = length(w);
w = w ./ sum(w);
Q = cumsum(w);
outIndex = zeros(1, N);
T = linspace(0,1-1/N,N) + rand(1)/N;

i=1;
j=1;
while (i<=N && j<=M),
    while Q(j)<T(i)
        j=j+1;
    end
    outIndex(i)=j;
    i=i+1;
end