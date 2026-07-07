function D = pairwiseDist(X,Y)
%PAIRWISEDIST: computes pairwise euclidean distance between sets X and Y
%   X and Y must have sizes nPoints*nDimensions

nX = size(X,1);
nY = size(Y,1);
D = zeros(nX,nY);
for i = 1:nX
    for j = 1:nY
        D(i,j) = norm(X(i,:) - Y(j,:),2);
    end
end

end

