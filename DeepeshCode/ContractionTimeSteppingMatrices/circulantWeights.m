function weights = circulantWeights(weightsVector)

    % Convert to a row vector
    weightsVector = (weightsVector(:))';

    % Find dimensions of circulant matrix
    Dim = length(weightsVector);
    
    % pint + 1
    pintP1 = length(find(weightsVector));
    
    % Full circulant matrix
    c = zeros(Dim);
    c(mod(Dim*Dim*(Dim+1)-1:-Dim-1:0,Dim*Dim)+1)=rectpulse(weightsVector,Dim);
    c = c';
    
    % Extract relevant rows
    weights = c(1:pintP1:end,:);

end