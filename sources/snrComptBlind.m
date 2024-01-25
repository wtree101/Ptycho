function [snr,scale]=snrComptBlind(I,uu)
%% uu as reference

scale=sum(uu(:).*conj(I(:)))/norm(I,'fro')^2;
I=scale* I;
snr=-20*log10(norm(I - (uu), 'fro')/norm(I, 'fro'));

end