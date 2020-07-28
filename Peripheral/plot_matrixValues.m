function plot_matrixValues(values)

textStrings = num2str(values(:), '%0.3f');       % Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  % Remove any space padding
[x, y] = meshgrid(1:(size(values,2)), 1:(size(values,1)));  % Create x and y coordinates for the strings
hStrings = text(x(:), y(:), textStrings(:), ...  % Plot the strings
                'HorizontalAlignment', 'center');
%textColors = zeros(size(image,2)*size(image,1),3); 
textColors = repmat(values(:) < 0.05, 1, 3);  % Choose white or black for the text color of the strings so
set(hStrings, {'Color'}, num2cell(textColors, 2));  % Change the text colors

     
end