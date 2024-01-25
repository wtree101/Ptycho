function co = coherence(masks)
[px,py,modes] = size(masks);
q = reshape(masks,[px*py modes]);
co = 0;
for i=1:modes
    for j=1:modes
        if (i~=j)
            co = max( (abs(q(:,j)'*q(:,i)))/(norm(q(:,i),2)*norm(q(:,j),2)),co );
        end
    end
end