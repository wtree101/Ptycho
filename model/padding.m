function ker = padding(ker0,Nx)
[kerh,kerw] = size(ker0);

select = (ker0~=0);
ker0(select) = 1;
h_st = fix(Nx/2 - kerh/2);
w_st = fix(Nx/2 - kerw/2);

ker = zeros(Nx,Nx);
ker(h_st+1:h_st+kerh,w_st+1:w_st+kerw) = ker0;

end