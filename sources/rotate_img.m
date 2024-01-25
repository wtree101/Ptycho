function J=rotate_img(I,theta,grey)
J = imrotate(I,theta,'crop');
J(J==0) = grey;
%imshow(J);
% kernel=kernel/sum(kernel,'all');
end