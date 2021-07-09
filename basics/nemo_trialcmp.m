function [goodtrials,badtrials] = nemo_trialcmp(dataorig,dataclean)
% compares two versions of a data structure before and after trial rejection
% add returns a list of "good trials" and "bad trials".
% Useful for re-applying the existing trial rejection to the data after revising a small pre-processing detail

badtrials = find(ismember(dataorig.sampleinfo,dataclean.cfg.artfctdef.summary.artifact,'rows'));
goodtrials = find(~ismember(dataorig.sampleinfo,dataclean.cfg.artfctdef.summary.artifact,'rows'));
