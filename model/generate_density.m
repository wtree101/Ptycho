function density = generate_density(masks)
[px py modes] = size(masks);
q = reshape(masks,[px*py modes]);

density = q * q';

end