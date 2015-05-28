%% Generate_Experiment.m
%
% This function will generate a computational experiment that can screen
% the parameters of a model for their importance. The analysis of the
% results is performed by a second function.
%
% args:
%   k = # parameters
%   p = # divisions for each parameter
%   r = # random orientations
%
% returns:
%

function experiments = Generate_Experiment(k,p,r)

if rem(k,2)
    error('Only valid for EVEN number of parameters')
end

% initialise design matrices
B = tril(ones(k+1,k),-1);
D_ = diag(randi(2,1,k)*2-3); % random diagonal
J = ones(k+1,k);

% generate random orientation
for i=1:r
    x_ = (randi(p-1,1,k)-1)./(p-1); % pick random parameter vector
    P_ = eye(k);
    P_ = P_(randperm(k),:); % random permutation matrix

    % random orientation of B
    experiments{i} = (J(:,1)*x_ + (1/(2*(p-1))).*((2*B-J)*D_+J))*P_;
end

end