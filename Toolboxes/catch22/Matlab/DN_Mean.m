function out = DN_Mean(y)

    % no combination of single functions
    coder.inline('never');
    
    out = mean(y);%DN_Mean(y); % took out the 5 since only takes one argument