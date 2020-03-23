function [splitVeh] = SplitFunc(nov,tr)

noc = length(tr); % Number of Choices (noc), the number of directions the cars can go in.

splitVeh = zeros(noc+1,1); % Pre determine size of output matrix

for i = 1:noc
    
    ranges = cumsum(binopdf(0:nov,nov,tr(i))); % Calculate the cumulative ranges for the binomial distribution
    
    randVar = rand; % Generate random variable for a random number of vehicles 
    
    splitVeh(i) = nov + 1 - sum(ranges >= randVar); % Determine the number of vehicles turning
    
    nov = nov - splitVeh(i); % Update remaining number of vehicles
    
end

splitVeh(end) = nov; % Add remaining number of vehicles for main direction

end

