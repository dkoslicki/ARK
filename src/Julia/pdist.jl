#Pairwise euclidean distance function
function pdist(X::Matrix, Y::Matrix)
	m = size(X,1);
	n = size(Y,1);
	Yt = Y';
	XX = sum(X.*X,2);
	YY = sum(Yt.*Yt,1);
	D = XX[:,ones(1,n)] + YY[ones(1,m), :] - 2*X*Yt;
	return D
end