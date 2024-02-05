function pphi = reconstruct_zeroforms(PHI,hEG)

% only for dual grid
% reconstruct 0-forms

global N

phi = reshape(PHI,N,N);

pphi = hEG'*phi*hEG;