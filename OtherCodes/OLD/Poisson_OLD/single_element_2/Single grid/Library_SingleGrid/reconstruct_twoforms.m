function pphi = reconstruct_twoforms(PHI,e,J)

% only for single grid
% reconstruct 2-forms

global globalnr_2

phi = PHI(globalnr_2);

pphi = (e'*phi*e).*J;