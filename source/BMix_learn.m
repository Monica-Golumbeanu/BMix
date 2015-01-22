% Wrapper for the learn_model function, applies the model on the input data
%
% Input: Signals_folder - path of the directory where the mismatch
%                         frequency profiles are located
%        Output_folder - path of the folder where the results are saved
%        CovTHR - minimum coverage to consider
%        PosteriorTHR - posterior lower bound used for classification
%
% Written by Monica Golumbeanu
% monica.golumbeanu@bsse.ethz.ch
%
function BMix_learn(Signals_folder, Output_folder, CovTHR, PosteriorTHR)

% learn parameters for the forward strand
tc_struct = dir([Signals_folder 'TC_f.txt']);
ta_struct = dir([Signals_folder 'AC_f.txt']);
tg_struct = dir([Signals_folder 'GC_f.txt']);

tc_files = strcat(Signals_folder, {tc_struct.name});
ta_files = strcat(Signals_folder, {ta_struct.name});
tg_files = strcat(Signals_folder, {tg_struct.name});

learn_model(tc_files, [ta_files tg_files], Output_folder, CovTHR, PosteriorTHR);

% learn parameters for the reverse strand
ag_struct = dir([Signals_folder 'AG_r.txt']);
at_struct = dir([Signals_folder 'TG_r.txt']);
ac_struct = dir([Signals_folder 'CG_r.txt']);

ag_files = strcat(Signals_folder, {ag_struct.name});
at_files = strcat(Signals_folder, {at_struct.name});
ac_files = strcat(Signals_folder, {ac_struct.name});

learn_model(ag_files, [at_files ac_files], Output_folder, CovTHR, PosteriorTHR);

end
