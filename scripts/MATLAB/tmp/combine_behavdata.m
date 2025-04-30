function R = combine_behavdata(S1, S2)
% Assuming S1 and S2 are already defined
fields = fieldnames(S1);

% Initialize a new structure S
R = struct();

% Define the desired dimensions
desiredDimensions = [8, 68];

% Loop through each field and concatenate the corresponding fields of S1 and S2 if dimensions match
for i = 1:numel(fields)
    field = fields{i};
    
    % Check if the dimensions of the field in both S1 and S2 match the desired dimensions
%     if isequal(size(S1.(field)), desiredDimensions) && isequal(size(S2.(field)), desiredDimensions)
        % Concatenate along the first dimension (vertically)
    R.(field) = [S1.(field); S2.(field)];
%     end
end

% Display the result (for verification)
