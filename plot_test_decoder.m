 %[dec, deconv, x, pmod] = ccnl_decode_regressor(optCon_expt(), 6, 'RPE', 'masks/left_VS.nii', 1, 1)
 [dec, deconv, x, pmod] = ccnl_decode_regressor(optCon_expt(), 11, 'psi', 'masks/left_VS.nii', 1, 1)     
 
      dec = mean(dec{1},2);
      deconv = mean(deconv{1}, 2);
      x = x{1};
      pmod = pmod{1};
      figure; plot([dec x]); legend({'decoded x HRF', 'original x HRF'}); xlabel('TR');
      figure; scatter(dec, x); xlabel('decoded x HRF'); ylabel('original x HRF');
      [r, p] = corr(dec, x)
      figure; plot([deconv pmod]); legend({'deconvolved', 'original'}); xlabel('TR');
      figure; scatter(deconv, pmod); xlabel('deconvolved'); ylabel('original');
      [r, p] = corr(deconv, pmod)