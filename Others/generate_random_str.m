function rand_str = generate_random_str

% this function generates a string made of digits, usually of length 10,
% sometimes a bit less than 10.
% the first 5 digits are based on the current time,
% the second 5 digits are based on randomization
% the function's purpose is that files containing such a random string will
% be very very unlikely to have the same name as another existing file.

rand_str = strcat(num2str(floor(mod(now*10^10,10^5))) ,num2str(floor(rand*10^5)));

end