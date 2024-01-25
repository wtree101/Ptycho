function [snr,scale]=snrComptC(I,uu)
%% uu as reference

scale=exp(-1i*angle(trace(uu'*I)));
I=scale* I;
snr=-20*log10(norm(I - (uu), 'fro')/norm(I, 'fro'));

end