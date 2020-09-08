[Pu,V] = ccnl_load_mask('masks/Pu.nii');
Ca = ccnl_load_mask('masks/Ca.nii');
Na = ccnl_load_mask('masks/NAC.nii');
Str = Pu | Ca | Na;

Str_L = Str;
Str_R = Str;

Pu_L = Pu;
Pu_R = Pu;

Na_L = Na;
Na_R = Na;

L_xs = 1:round(size(Str,1)/2); % one hemisphere
R_xs = round(size(Str,1)/2)+1:size(Str,1); % one hemisphere
Str_R(L_xs,:,:) = 0;
Str_L(R_xs,:,:) = 0;
Pu_R(L_xs,:,:) = 0;
Pu_L(R_xs,:,:) = 0;
Na_R(L_xs,:,:) = 0;
Na_L(R_xs,:,:) = 0;

V.fname = 'masks/Str_R.nii';
spm_write_vol(V, Str_R);
V.fname = 'masks/Str_L.nii';
spm_write_vol(V, Str_L);
V.fname = 'masks/Pu_R.nii';
spm_write_vol(V, Pu_R);
V.fname = 'masks/Pu_L.nii';
spm_write_vol(V, Pu_L);
V.fname = 'masks/NAC_R.nii';
spm_write_vol(V, Na_R);
V.fname = 'masks/NAC_L.nii';
spm_write_vol(V, Na_L);


V.fname = '../Momchil/Pu_post.nii'; % !!!!!!!!!!!
Pu_y = squeeze(sum(sum(Pu,1),3));
miny = min(find(Pu_y));
maxy = max(find(Pu_y));
midy = round((miny * 0.8 + maxy * 0.2) );
Pu(:,midy:end,:) = 0;

spm_write_vol(V, Pu);
