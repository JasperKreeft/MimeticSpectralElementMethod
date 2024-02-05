function pphi = reconstruct_zeroforms(PHI,hGL)

% only for single grid
% reconstruct 0-forms

global globalnr_0

phi = PHI(globalnr_0);

pphi = hGL'*phi*hGL;