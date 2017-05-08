function [newA] = nearestRankFast(A,r)
    [n,m,l] = size(A);
    newA = zeros(n,m,l);
    disp(n);
    A_red = double(A(:,:,1));
    A_green = double(A(:,:,2));
    A_blue = double(A(:,:,3));
    [V,~] = powerMethod(A_red, r, 20);
    [U,T] = powerMethod(A_red', r, 20);
    for i=1:r
       [~, lt] = eig(T);
       newA(:,:,1) = newA(:,:,1) + sqrt(lt(i,i))*U(:,i)*V(:,i)';
    end
    [V,~] = powerMethod(A_green, r, 20);
    [U,T] = powerMethod(A_green', r, 20);
    for i=1:r
       [~, lt] = eig(T);
       newA(:,:,2) = newA(:,:,2) + sqrt(lt(i,i))*U(:,i)*V(:,i)';
    end
    [V,~] = powerMethod(A_blue, r, 20);
    [U,T] = powerMethod(A_blue', r, 20);
    for i=1:r
       [~, lt] = eig(T);
       newA(:,:,3) = newA(:,:,3) + sqrt(lt(i,i))*U(:,i)*V(:,i)';
    end
    newA = uint8(newA);
end