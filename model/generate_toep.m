function density = generate_toep(px,J) %J - J(0,x)

N = px*px;
density = zeros(N,N);
[x,y] = size(J);

dx = floor(x/2) + mod(x,2);
dy = floor(y/2) + mod(y,2);
    
for i=1:N
    for j=1:N
        
        [x1,y1] = ind2sub([px px],i);
        [x2,y2] = ind2sub([px px],j);
        %x1 = x1 - N/2; x2 = x2 - N/2; y1 = y1 - N/2; y2 = y2 - N/2;
        %distx = x2 - x1 + px/2; disty = y2 - y1 + px/2;
        distx = x2 - x1 + dx; disty = y2 - y1 + dy;
        if (1<=distx && distx<= x) && (1<=disty && disty<= y)
            indx = sub2ind([x y],distx,disty);
            density(i,j) = J(indx);
            %density(j,i) = J(indx);
        end
    end
end




end