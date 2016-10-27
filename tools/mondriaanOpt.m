function [I, s] = mondriaanOpt(A, Imbalance, Volume)
% mondriaan       Uses the MondriaanOpt algorithm to assign the nonzeros
%                 of a given sparse matrix A to 2 processors. The computed
%                 partitioning has the lowest communication volume among
%                 all partitionings that have imbalance at most epsilon.
%         
%         Required arguments:
%             A         = Sparse matrix
%             Imbalance = Maximum allowed imbalance
%             Volume    = The initial upper bound for the communication volume
%
%         Return values:
%             I = Copy of the matrix A where the numerical value of each
%                 nonzero is replaced by the index of the processor to
%                 which this nonzero was assigned to by the MondriaanOpt
%                 algorithm:
%                  1 = processor 0
%                  2 = processor 1
%                  3 = free nonzero
%             s = The statistics vector with information about the
%                 partitioning: duration, imbalance, communication volume
%         
%         Usage:
%             [I, s] = mondriaanOpt(A, Imbalance, Volume);
%
    
    % Check input
    if (nargin < 3)
        disp 'Error: the three arguments to mondriaanOpt function are mandatory!';
        return;
    end
    if (Imbalance <= 0)
        disp('Invalid value for Imbalance');
        return;
    end
    if (issparse(A) ~= 1)
        disp('Input matrix A is not sparse');
        return;
    end
    
    % Run MondriaanOpt
    tic;
    comVol = MatlabMondriaanOpt(A, Imbalance, Volume);
    elapsedTime = toc;
    
    if(comVol > Volume)
        throw(MException('mondriaanOpt:NoSol', ['No solution with a volume at most ' num2str(Volume) ' exists!']));
    end
    
    % Read computed partitioning from disk
    I = mmread('MatlabMex.mtx-I2f');
    
    % Compute the load on the processors
    loadCUT = sum(sum(I==3));
    loadP1  = sum(sum(I==1));
    loadP2  = sum(sum(I==2));
    
    % Distribute the cut elements evenly
    % First, make both parts even
    if(loadP1 < loadP2)
        add = min(loadCUT, loadP2-loadP1);
        loadP1 = loadP1 + add;
        loadCUT = loadCUT - add;
    elseif(loadP2 < loadP1)
        add = min(loadCUT, loadP1-loadP2);
        loadP2 = loadP2 + add;
        loadCUT = loadCUT - add;
    end
    % Now cut the rest in two
    loadP1 = loadP1 + floor(loadCUT/2);
    loadP2 = loadP2 + ceil(loadCUT/2);
    
    % Determine the best max load we can achieve
    maxLoad = max(loadP1,loadP2);
    
    % Compute epsilon
    nonz = nnz(I);
    if mod(nonz,2) == 0
        epsilon = 2*maxLoad/nonz - 1;
    else
        epsilon = 2*maxLoad/(nonz+1) - 1;
    end
    
    s = [elapsedTime; epsilon; comVol];

