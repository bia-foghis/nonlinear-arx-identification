clear; close all;
load("iddata-07.mat");

u_id = id.InputData;
y_id = id.OutputData;

u_val = val.InputData;
y_val = val.OutputData;

%plot(id);
%title("Identification Data");
%figure
%plot(val);
%title("Validation Data");

%%
max_na_nb = 3;
max_m = 5;

MSE_hat_final = zeros(max_m);
MSE_sim_final = zeros(max_m);

for i = 1:max_m
    for j = 1:max_na_nb
        [MSE_hat,MSE_sim] = ARX(j,j,i);
        MSE_sim_final(i,j) = MSE_sim;
        MSE_hat_final(i,j) = MSE_hat;
    end
end

%%
function [MSE_hat,MSE_sim] = ARX(na, nb, m)

load("iddata-07.mat");

u_id = id.InputData;
y_id = id.OutputData;

u_val = val.InputData;
y_val = val.OutputData;

nk = 1; % delay

% generate all combinations of values between 0 and m
combinations = [];

n = na + nb; % number of elements in the vector
help_comb = zeros(1, n);

while true
    combinations = [combinations; help_comb];
    
    k = n;
    while k > 0 && help_comb(k) == m
        k = k - 1;
    end
    
    if k == 0
        break;
    end
    
    help_comb(k) = help_comb(k) + 1;
    help_comb(k + 1:end) = 0;
end

combinations = combinations(sum(combinations, 2) <= m, :);

%disp(combinations);

PHI_id = zeros(length(y_id),length(combinations));

[nr_of_lines,nr_of_columns] = size(combinations);

for i = 1:length(y_id)
    PHI_line = zeros(nr_of_columns,1);

    for j = 1:nr_of_lines
        val = 0;

        for x = 1:nr_of_columns
            % for y
            if(x <= nr_of_columns/2)
                if(val == 0 && combinations(j,x) ~= 0)
                    if(i - (mod(x-1,na)+1) >= 1)
                        val = val + y_id(i - (mod(x-1,na)+1))^combinations(j,x);
                    end
                else
                    if(i - (mod(x-1,na)+1) >= 1)
                        val = val * y_id(i - (mod(x-1,na)+1))^combinations(j,x);               
                    end
                end
            
            % for u
            elseif(x > nr_of_columns/2)
                if(val == 0 && combinations(j,x) ~= 0)
                    if(i - (mod(x-1,nb)+1) >= 1)
                        val = val + u_id(i - (mod(x-1,nb)+1) - nk + 1)^combinations(j,x);   
                    end
                else
                    if(i - (mod(x-1,nb)+1) >= 1)
                        val = val * u_id(i - (mod(x-1,nb)+1) - nk + 1)^combinations(j,x);    
                    end
                end
            end
        end
        PHI_line(j) = val;
    end
    PHI_id(i,:) = PHI_line;
end

PHI_id(:,1) = 1;
PHI_id(1,1) = 0;

theta = PHI_id\y_id;

y_hat_id = PHI_id*theta;

PHI_val = zeros(length(y_val),length(combinations));

[nr_of_lines,nr_of_columns] = size(combinations);

for i = 1:length(y_val)
    PHI_line = zeros(nr_of_columns,1);

    for j = 1:nr_of_lines
        val = 0;

        for x = 1:nr_of_columns
            % for y
            if(x <= nr_of_columns/2)
                if(val == 0 && combinations(j,x) ~= 0)
                    if(i - (mod(x-1,na)+1) >= 1)
                        val = val + y_val(i - (mod(x-1,na)+1))^combinations(j,x);
                    end
                else
                    if(i - (mod(x-1,na)+1) >= 1)
                        val = val * y_val(i - (mod(x-1,na)+1))^combinations(j,x);               
                    end
                end
            
            % for u
            elseif(x > nr_of_columns/2)
                if(val == 0 && combinations(j,x) ~= 0)
                    if(i - (mod(x-1,nb)+1) >= 1)
                        val = val + u_val(i - (mod(x-1,nb)+1) - nk + 1)^combinations(j,x);   
                    end
                else
                    if(i - (mod(x-1,nb)+1) >= 1)
                        val = val * u_val(i - (mod(x-1,nb)+1) - nk + 1)^combinations(j,x);    
                    end
                end
            end
        end
        PHI_line(j) = val;
    end
    PHI_val(i,:) = PHI_line;
end

PHI_val(:,1) = 1;
PHI_val(1,1) = 0;

y_hat_val = PHI_val*theta;

y_sim_id = zeros(length(y_id),1);
theta_idx = 0;

for i = 1:length(y_id)
    y_sim_id_value = 0;

    if(theta_idx == length(theta))
        theta_idx = 0;
    end

    for j = 1:length(combinations)
        val = 0;
        
        for x = 1:nr_of_columns
            % for y
            if(x <= nr_of_columns/2)
                if(val == 0 && combinations(j,x) ~= 0)
                    if(i - (mod(x-1,na)+1) >= 1)
                        val = val + y_sim_id(i - (mod(x-1,na)+1))^combinations(j,x);                  
                    end
                else
                    if(i - (mod(x-1,na)+1) >= 1)
                        val = val * y_sim_id(i - (mod(x-1,na)+1))^combinations(j,x);                
                    end
                end
            
            % for u
            elseif(x > nr_of_columns/2)
                if(val == 0 && combinations(j,x) ~= 0)
                    if(i - (mod(x-1,nb)+1) >= 1)
                        val = val + u_id(i - (mod(x-1,nb)+1) - nk + 1)^combinations(j,x);  
                    end
                else
                    if(i - (mod(x-1,nb)+1) >= 1)
                        val = val * u_id(i - (mod(x-1,nb)+1) - nk + 1)^combinations(j,x);                  
                    end
                end
            end
        end
        theta_idx = theta_idx + 1;
        y_sim_id_value = y_sim_id_value + (val * theta(theta_idx));        
    end
    y_sim_id(i) = y_sim_id_value;
end

y_sim_val = zeros(length(y_val),1);
theta_idx = 0;

for i = 1:length(y_val)
    y_sim_val_value = 0;

    if(theta_idx == length(theta))
        theta_idx = 0;
    end

    for j = 1:length(combinations)
        val = 0;
        
        for x = 1:nr_of_columns
            % for y
            if(x <= nr_of_columns/2)
                if(val == 0 && combinations(j,x) ~= 0)
                    if(i - (mod(x-1,na)+1) >= 1)
                        val = val + y_sim_val(i - (mod(x-1,na)+1))^combinations(j,x);                  
                    end
                else
                    if(i - (mod(x-1,na)+1) >= 1)
                        val = val * y_sim_val(i - (mod(x-1,na)+1))^combinations(j,x);                
                    end
                end
            
            % for u
            elseif(x > nr_of_columns/2)
                if(val == 0 && combinations(j,x) ~= 0)
                    if(i - (mod(x-1,nb)+1) >= 1)
                        val = val + u_val(i - (mod(x-1,nb)+1) - nk + 1)^combinations(j,x);  
                    end
                else
                    if(i - (mod(x-1,nb)+1) >= 1)
                        val = val * u_val(i - (mod(x-1,nb)+1) - nk + 1)^combinations(j,x);                  
                    end
                end
            end
        end
        theta_idx = theta_idx + 1;
        y_sim_val_value = y_sim_val_value + (val * theta(theta_idx));        
    end
    y_sim_val(i) = y_sim_val_value;
end

% MSE
MSE_hat = 0;
MSE_sim = 0;

for i = 1:length(y_hat_val)
    MSE_hat = MSE_hat + (y_val(i) - y_hat_val(i))^2;
    MSE_sim = MSE_sim + (y_val(i) - y_sim_val(i))^2;
end

MSE_hat = 1/length(y_hat_val) * MSE_hat;
MSE_sim = 1/length(y_sim_val) * MSE_sim;
