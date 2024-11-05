function [kr1, kr2, R1, R2, z] = opt_data(m,opt,N)

% These are the different wave-vector ratios used for the excitation beam 
Ex = [0.6, 0.65, 0.7, 0.75 1]; 

% These are the different wave-vector ratios used for the depletion beam 
Dp = [0.7, 0.75, 0.8, 0.85, 1];

    if opt == 1 % excitation case
        % The optimised kr1, kr2, A1 and A2 for excitation beam are available in following csv file 
        A = readmatrix('excitation.csv');
        for i = 1 : 5
            if Ex(i) == m
                kr1 = A(i,2);
                kr2 = kr1 * m;
                rng('default');
                R1 = rand(N);
                R1(R1>A(i,5))=1;
                R1(R1<=A(i,5))=0;
                R2 = abs(1-R1);
                z = A(i,6) * 10^-6; % propagation distance is in micrometers
            end
        end
    elseif opt == 2 % depletion case
        % The optimised kr1, kr2, A1 and A2 for depletion beam are available in following csv file 
        A = readmatrix('depletion.csv');
        for i = 1 : 5
            if Dp(i) == m
                kr1 = A(i,2);
                kr2 = kr1 * m;
                rng('default');
                R1 = rand(N);
                R1(R1>A(i,5))=1;
                R1(R1<=A(i,5))=0;
                R2 = abs(1-R1);
                z = A(i,6) * 10^-6; % propagation distance is in micrometers
            end
        end
    end

end
