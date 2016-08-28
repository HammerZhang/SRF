function [ M ] = srf( Me,omega,dec_fac, err_threshold, num_loop,rank_val )
% SRF: provide a solution for matrix completion problem with a novel
% algorithm named smoothed rank function. The result is a approximation of
% the input maxtrix Me with minium rank value.
% Usage: 
% >> [M] = srf(me, dec_fac, err_threshold,num_loop,omega,rank_val)
% Inputs:
% Me: - initial matrix with full rank as the starting matrix;
% omega: - position set of revealed entries, the value of elements in omega
% are ones if the corresponding postition of Me has entities, the rest are
% zeros
% dec_fac: - decreasing factor of delta
% err_threshold: - err threshold to stop the algorithm
% num_loop: - number of inner loop
% rank_val: - rank value of initial matrix
% 
% Outputs:
% M: results of matrix completion
%
% @ Author: ...
% @ Time: 2016.6.1 23:46
%
%% ======================================================================

    if nargin < 2
        error('Initial matrix and entries position set must be given....\n');
    end

    if nargin < 3
        dec_fac = 0.9;
    end

    if nargin < 4
        err_threshold = 1e-11;
    end

    if nargin < 5
        num_loop = 8;
    end

    if nargin < 6
        % indicates that rank_val is not given
        rank_val = -1;
    end


    % Initialization
    [~,Sm,~] = svd(Me);                         % svd of initial matrix Me
    delta = Sm(1,1);                            % practically initial value of delta is 3 or 4 times of largest sigular value of Me;
    d = inf;
    X = Me;   [n1,n2] = size(X);                 % initial X and size
    %mu = 250;                                     % step size
    
    % start loop
    while d > err_threshold
        Y = X;
        % Inner loop
        for ii = 1:num_loop
            % compute svd of X, if rank_val is given, use svds, else use
            % svd
            if rank_val == -1
                [U,S,V] = svd(X);
            else
                [U,S,V] = svds(X);
            end
            
            DiagF = zeros(size(S));
            mu = sum(sum(S,1),2)/min(size(S))*0.75;
            for ni = 1:size(S,1)
                for mi = 1:size(S,2)
                    if ni == mi
                        DiagF(ni,mi) = S(ni,mi) - mu*S(ni,mi)/delta^2 * exp(-S(ni,mi)^2/(2*delta^2));
                    end
                end
            end
            
            Xt = U*DiagF*V';
            % compute projection
            X = computeProj(Me,Xt,omega);
        end
        
        % update parameters
        delta = delta * dec_fac;
        d = norm(Y-X,'fro');
        fprintf('dev is %s, delta is %s\n',d,delta);
        %d = d/sqrt(n1*n2);
    end
    
    % final result
    M = X;


end


function [Xp] = computeProj(M,X,W)
% COMPUTEPROJ compute a projection from matrix X to matrix Xp, the project
% is defined in that paper
%
% @Author;...
% @Time: 2016.6.2.0:38
%
%% =========================================================================

    if size(M) ~= size(X)
        error('matrix M and X must have the same length...\n');
    end

    entities = find(W == 1);
    noentities = find(W == 0);

    Xp = zeros(size(M));
    Xp(entities) = M(entities);
    Xp(noentities) = X(noentities);

end




