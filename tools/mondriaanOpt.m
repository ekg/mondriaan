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
%                 which this nonzero was assigned to by the Mondriaan
%                 algorithm.
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
    MatlabMondriaanOpt(A, Imbalance, Volume);
    elapsedTime = toc;
    
    % Read computed partitioning from disk
    I = mmProcIndexRead('MatlabMex_P2.mtx');
    
    % Compute the load on the processors
    loadCUT = sum(sum(I==3));
    loadP1  = sum(sum(I==1));
    loadP2  = sum(sum(I==2));
    
    % Distribute the cut elements evenly
    loadP1 = loadP1 + floor(loadCUT/2);
    loadP2 = loadP2 + floor(loadCUT/2);
    if mod(loadCUT, 2) == 1
       if(loadP1 < loadP2)
           loadP1 = loadP1 + 1;
       elseif(loadP2 < loadP1)
           loadP2 = loadP1 + 1;
       end
    end
    
    % Determine the best max load we can achieve
    maxLoad = max(loadP1,loadP2);
    
    % Compute epsilon
    nonz = nnz(I);
    if mod(nonz,2) == 0
        epsilon = 2*maxLoad/nonz - 1;
    else
        epsilon = 2*maxLoad/(nonz+1) - 1;
    end
    
    % Compute communication volume
    comVol = sum((sum(I == 1, 1)>0)+(sum(I == 2, 1)>0) > 1) + ...
             sum((sum(I == 1, 2)>0)+(sum(I == 2, 2)>0) > 1);
    
    s = [elapsedTime; epsilon; full(comVol)];



% mmProcIndexRead Function based on mmread(), reading a matrix-market file,
%                 and turning it into a processor index matrix. Compared to
%                 mmread(), this fixes two things:
%                 1) mmread() does not recognise the 'integer' field,
%                 2) add 1 to all values of the matrix, because otherwise
%                    all entries of processor 0 disappear by sparsity of
%                    the matrix.
%         
%         Required arguments:
%             filename  = The name of the file to read
%
%         Return values:
%             A = The matrix A where the numerical value of each
%                 nonzero is the index of the processor to which this
%                 nonzero was assigned to by the MondriaanOpt algorithm:
%                 1 = Processor 0
%                 2 = Processor 1
%                 3 = A cut element, may be assigned to either 0 or 1
%                     without affecting the communication volume
%         
%         Usage:
%             [A] = mmProcIndexRead(filename);
%
function [A] = mmProcIndexRead(filename)
    mmfile = fopen(filename,'r');
    if ( mmfile == -1 )
     disp(filename);
     error('File not found');
    end;

    header = fgets(mmfile);
    if (header == -1 )
      error('Empty file.')
    end

    % Read through comments, ignoring them

    commentline = fgets(mmfile);
    while length(commentline) > 0 && commentline(1) == '%',
      commentline = fgets(mmfile);
    end

    % Read size information
    [sizeinfo,count] = sscanf(commentline,'%d%d%d');
    while ( count == 0 )
     commentline =  fgets(mmfile);
     if (commentline == -1 )
       error('End-of-file reached before size information was found.')
     end
     [sizeinfo,count] = sscanf(commentline,'%d%d%d');
     if ( count > 0 & count ~= 3 )
       error('Invalid size specification line.')
     end
    end
    rows = sizeinfo(1);
    cols = sizeinfo(2);
    entries = sizeinfo(3);

    [T,~] = fscanf(mmfile,'%f',3);
    T = [T; fscanf(mmfile,'%f')];
    if ( size(T) ~= 3*entries )
       message = ...
       str2mat('Data file does not contain expected amount of data.',...
               'Check that number of data lines matches nonzero count.');
       disp(message);
       error('Invalid data.');
    end
    T = reshape(T,3,entries)';
    A = sparse(T(:,1), T(:,2), T(:,3)+1, rows , cols);


