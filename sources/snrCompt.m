function snr=snrCompt(I,uu)
%snr=-20*log10(norm(I - real(uu), 'fro')/norm(real(uu), 'fro'));
%snr=-20*log10(norm(I - (uu), 'fro')/norm(uu, 'fro'));
snr=-20*log10(norm(I(:) - (uu(:)))/norm(I(:)));

end