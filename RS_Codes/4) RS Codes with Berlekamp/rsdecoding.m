clear;
close all;
%% Parameters
p = 8;     %GF(2^p)-field
pp = primpoly(p,'min','nodisplay');
n = 2^p -1;  %Size of code word
t = 16;      %Error Correcting Capacity
k = n - 2*t; %Size of Message

%% Finding all possible primitive elements
ff = 0:1:2^p -1;
ff = gf(ff,p,pp);
ff_star = ff(1,2:end);

factor_n = find_factors(n);
ord_element = zeros(1,length(ff_star));

for i = 1:length(ff_star)
    for j = 1:2^p-1
        if ff_star(1,i)^j == 1
            ord_element(1,i) = j;
            break;
        end
    end
end

% Putting All primitive element into a list
prim_ele = gf(zeros(1,0),p,pp);
for i = 1:length(ord_element)
    if ord_element(1,i)==n
        prim_ele(1,end+1) = ff_star(1,i);
    end
end

%% Power Table;
% alpha - Chosen primitive Element
alpha = prim_ele(1,1);
power_ele= gf(zeros(1,2^p),p,pp);

for i = 1:2^p
    ele = alpha^(i-1);
    power_ele(1,i) = ele;
end

power = 0:1:n;
power_ele = power_ele.x;

power_to_ele_dict = containers.Map(power,power_ele);
ele_to_power_dict = containers.Map(power_ele,power);

%% Generating Data
message = gf(zeros(1,0),p,pp);

for i = 1:k
    indx = randi([1 length(ff)],1);
    r_m  = ff(1,indx);
    message(1,end+1) = r_m;
end

%% Generator Polynomial

g = [1];
for i = 1:2*t
    g = conv(g,[alpha^i 1]);
end

% Generator Polynomial is expressed in terms of array g

%% Encoding (systematic)

% Implemented systematic encoding
% Top k Symbols can be directly taken for the message

% message_app = [zeros(1,n-k) message];
% 
% 
% % Beware : deconvolution is taking higher power to be at zero at lowest
% % power at end for Both Inputs and Outputs
% 
% [q,r] = deconv(flip(message_app),flip(g));
% fprintf("\nQuotient and Remainder\n")
% r_rev = flip(r);
% 
% code_word = message_app + r_rev;

message_app = [message zeros(1,n-k)];

% Beware : deconvolution is taking higher power to be at zero at lowest
% power at end for Both Inputs and Outputs

[q,r] = deconv(message_app,flip(g));
fprintf("\nQuotient and Remainder\n")
% disp(r.x);
% r_rev = flip(r);
% disp(message_app.x);

code_word = message_app + r;
% disp(code_word.x);

% now highest degree is present in 0th index and lowest degree coeffient in
% present in index end.

%% Adding Errors

% Randomly choosing size_of_error error locations
num_list = 1:1:n;
size_of_error = t; % Select How much error to add
i_2t = datasample(num_list, size_of_error, 'Replace', false);
i_2t = sort(i_2t);

% Randomly selecting error value for each error location
err_val_vec = randi([0 n],1,length(i_2t));
err_val_vec = gf(err_val_vec,p,pp);

rx_word = code_word;
for i=1:length(i_2t)
    indx = i_2t(1,i);
    rx_word(1,indx) = rx_word(1,indx) + err_val_vec(1,i);
end
%% Decoding
[codeword_rx,flag] = rsDecode(rx_word,alpha,n,t,p,pp,power_to_ele_dict,ele_to_power_dict);