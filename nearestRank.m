function [newA] = nearestRank(A,r)
    [n,m,l] = size(A);
    newA = zeros(n,m,l);
    disp(n);
    A_red = double(A(:,:,1));
    A_green = double(A(:,:,2));
    A_blue = double(A(:,:,3));
    [U,S,V] = svd(A_red,'econ');
    disp(size(U))
    disp(size(V))
    for i=1:r
       newA(:,:,1) = newA(:,:,1) + S(i,i)*U(:,i)*V(:,i)';
    end
    [U,S,V] = svd(A_green,'econ');
    for i=1:r
       newA(:,:,2) = newA(:,:,2) + S(i,i)*U(:,i)*V(:,i)';
    end
    [U,S,V] = svd(A_blue,'econ');
    for i=1:r
       newA(:,:,3) = newA(:,:,3) + S(i,i)*U(:,i)*V(:,i)';
    end
    newA = uint8(newA);
end