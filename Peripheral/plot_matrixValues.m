function plot_matrixValues(values)

if isa(values,'double')
    textStrings = num2str(values(:), '%0.3f');    % Create strings from the matrix values
else
    textStrings = values(:);
end

textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
textStrings = strrep(textStrings(:),'NaN','-'); % replace NaN with
[x, y] = meshgrid(1:(size(values,2)), 1:(size(values,1)));  % Create x and y coordinates for the strings
    
if isa(values,'double')
    hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
        'HorizontalAlignment', 'center', 'FontWeight', 'Bold');
else
    hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
        'HorizontalAlignment', 'center');
end
textColors = repmat(0.00,length(values(:)), 3);
if isa(values,'double')
    diffCol = find(values(:) < 0.05);
    vals = ones(length(diffCol),3);
    textColors(diffCol,:) = vals;
end
% Choose white or black for the text color of the strings so
set(hStrings, {'Color'}, num2cell(textColors, 2));  % Change the text colors

end