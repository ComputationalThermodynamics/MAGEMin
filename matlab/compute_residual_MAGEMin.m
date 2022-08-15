clear;
close all;

data = fileread("./output/_residual_norm.txt");
data = cellstr(data);
data = cellfun(@(newline) strsplit(newline, '\n'), data, 'UniformOutput', false);
data = [data{:}];
data = data';

% str2num(string(data{1}))



figure(1)

for i=1:1:size(data,1)
    C = str2num(string(data{i}));
    
    up = min(size(C,2),128);
    plot(log10(C(1,2:up)),'-');
ylim([-6 1])
xlim([0 256])
    hold on
end
