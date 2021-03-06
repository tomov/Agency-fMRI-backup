function bspmview_wrapper(map)

EXPT = optCon_expt();

V = spm_vol(fullfile('../Momchil', 'spmT_0001.nii'));

% hacks to make it save the t-map as a t-map
V.fname = fullfile('../Momchil', ['temp_map.nii']); % change immediately!
V.dt = [16 0];
V.private.dat.dtype = 'FLOAT32-LE';
V.private.dat.fname = V.fname;

% save map
V.fname
spm_write_vol(V, map);

% view map
struc = fullfile(EXPT.modeldir,'mean.nii');
if exist(struc,'file')
    bspmview(V.fname, struc);
else
    bspmview(V.fname);
end
