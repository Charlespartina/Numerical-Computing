function testGooglematrix(M, url)
    list = [0.85,0.5,0.25];
    colors = {'rx','b*','go'};
    
    figure();
    hold on;
    for i=1:3
        [A,G,X,lambda] = googlematrix(M,list(i));
        % Get the eigenvector corresponding to the eigenvalue one
        [Y,I] = max(lambda);
        x = X(:,I);
        % Normalize the vector
        x = x/norm(x, 1);
        % Check the sign
        if x(1) < 0
           x = -x;
        end
        % Sort the components
        [xsort, index] = sort(x, 'descend');
        fprintf('alpha = %d\n',list(i));
        
        % Get the highest 10 ranked web pages
        for j=1:10
            fprintf('index = %d ',index(j));
            disp(url(index(j)));
        end
        
        % Plot eigenvalues
        plot(real(lambda),imag(lambda),colors{i});
    end
    axis equal;
    legend('alpha=0.85','alpha=0.5','alpha=0.25');
    xlabel('Real part of eigenvalues');
    ylabel('Imagin part of eigenvalues');
    title('Eigenvalues of G for different damping parameter alpha');
    hold off;
    
    % Compare eigenvalues using power method
    for i=1:3
        % Use eig
        [A,G,X,lambda] = googlematrix(M,list(i));
        % Use power method
        [v] = powerMethod(A,list(i));
        % Find the eigenvector corresponding to eigenvalue 1
        [Y,I] = max(lambda);
        x = X(:,I);
        % Normalize the vector
        x = x/norm(x, 1);
        if x(1) < 0
           x = -x;
        end
        % Output the second largest eigenvalue
        sort(lambda,'descend');
        fprintf('The second largest eigenvalue for alpha=%d is %d\n',list(i),lambda(2));
    end
end