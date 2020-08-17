

[mask, Vmask, Ymask] = load_mask('masks/I_POS_kda10_testStat_mean_cMass_testStat_thresh.nii');
Vmask.fname = 'masks/bartra_pos.nii';
spm_write_vol(Vmask, mask);
spm_write_vol(Vmask, mask);
view_mask('masks/bartra_pos.nii');