function J=translate_img(I,x,y,grey)
J = imtranslate(I,[x,y],'FillValues',grey);
%J = imtranslate(I,[x,y]);
%imshow(J);
% kernel=kernel/sum(kernel,'all');
end

