function [best_design] = ahp_PDD()
%% AHP analytical hierarchy process
%
%  To run, just load script into editor and hit the run key!
%
% This simple by example function (with default values) gives the basic 
% elements of the Analytical Hierarchial Process (AHP) for decision making 
% to include matrix formulations, pairwise analysis, calculating 
% eigenvectors, and determining the final 'best' decision based on criteria. 

%% Example problem:
% Situation: Selecting the best design for PDD based on criteria 
% created by our group through the creation of a tree diagram. 
% decision, I will use AHP with the following:
%         alternatives:  Design 1, 2, 3, 4, 5 ...
%         criteria:      tbd
%%---------------------------------------------

%% Problem formulation:
clear all; close all; clc;
%% Step 1: Criteria Matrix and Criteria Eigenvector
% 
% This part is subjective. We must rank the criteria (subjective 
% qualitative) using Scale:
%   1-equal, 3-moderate, 5-strong, 7-very strong, 9-extreme
i = [];
PCM = [];
PCM(:,:) = input(['Please input the pairwise comparison matrix.']);
% See example:
% pairwise comparision of each criteria to each criteria denoted as the
% matrix PCM
    % (S)   (R)  (FE)
    % PCM= [ 1/1   1/2   3/1; ...   (S-style)
    % 2/1   1/1   4/1; ...    (R-reliability)
    % 1/3   1/4   1/1 ]       (FE-fuel economy)
disp('Criteria Pairwise Comparison Matrix PCM');
disp(PCM);

   %ePCM=eig(PCM)
   
   ePCM=calc_eig(PCM);
   disp(ePCM);

   % Take the top 11 highest ranked criteria and put them in order of score

   % ePCM = unique(ePCM);
   % ePCM = ePCM((end-10):end)';
   % ePCM = fliplr(ePCM);
   % ePCM = ePCM'
%%
   

%% Step 2: Alternatives Matrix and Alternatives Eigenvectors 
% Ranking the laternative designs 
% Compair each design to each design denoted as the alternatives 
% (Design 1, 2, 3, 4, 5, etc...)
% In the terms of style, reliability, and fuel economy using
% pairwise comparison for qualitative and normalization for quantitative
% 
%   an example of a style comparison between cars:
%
%                    civic        focus         corolla         BMW318
% civic               1/1         1/4            4/1            1/6
% focus               4/1         1/1            4/1            1/4
% corolla             1/4         1/4            1/1            1/5
% BMW318              6/1         4/1            5/1            1/1

%% Qualitative Comparison matrices

% See example below:
% ACM_Re = [ 1/1 2/1 5/1 1/1; ...
%           1/2 1/1 3/1 2/1; ...
%           1/5 1/3 1/1 1/4; ...
%           1/1 1/2 4/1 1/1]

N_cr_Qual = input('Please input the number of  qualitative criteria');

for i = 1:1:N_cr_Qual
    ACM{i} = input(['Please input comparison matrix for qualitative criteria ' num2str( i)]);
    disp('Alternative Qualitative Pairwise');
    disp(ACM{i});
end


for i = 1:1:N_cr_Qual
    eACM{i} = calc_eig(ACM{i}); % calculate eigenvector on qualitative matrix
end


%% Quantitative Comparison Matrices

% Example format of quantitative comparison matrix:
% given MPG data for each vehicle, create a fuel economy matrix
% MPG matrix    
% cv = 34;   % civic
% sa = 27;   % focus
% es = 24;   % corolla
% cl = 28;   % BMW318
 
% ACM_Fe = [ cv; ...
%            sa; ...
%            es; ...
%            cl]

N_cr_quant = input('Please input the number of  quantitative criteria');

if N_cr_quant > 0;

    disp(i);

    if length(i)<1

        for i = 1:1:N_cr_quant
            disp(i)
            ACM{i} = input(['Please input comparison matrix for quantitative criteria ' num2str( i(end)-N_cr_Qual)]);
            disp('Alternatives Quantitative Pairwise');
            disp(ACM{i});
        end
    else
        for i = i(end)+1:1:i(end) + N_cr_quant
            disp(i)
            ACM{i} = input(['Please input comparison matrix for quantitative criteria ' num2str( i(end)-N_cr_Qual)]);
            disp('Alternatives Quantitative Pairwise');
            disp(ACM{i});
        end
    end

    for i = i(end) + 1 - N_cr_quant:1:i(end)
        eACM{i} = calc_norm(ACM{i}); % normalize quantitative type data
    end


else

end
     
%% Step 3:  Calculate Final Answer and Determine winner
% construct a matrix of eigenvectors calculated above for each criteria

eM = cell2mat(eACM);



% multiply eigenvector matrix by eigenvector of criteria to obtain
% scores for each car based on criteria and car-to-car comparisons
disp('Scores for alternatives 1, 2, 3, 4, 5, etc ...')
ePCM = transpose(ePCM);
Design_Scores = eM .* ePCM

scores_total = [];

% Best design based on criteria
for i = 1:height(Design_Scores)
    scores_total(i) = sum(Design_Scores(i,:));
end
disp(scores_total);
[pos,~] = find(max(scores_total));
disp(['Winning design is design ' num2str(pos)]);
best_design = max(scores_total);  


%% sub-function on calculating eigenvectors
    function [eigvect ] = calc_eig(M) 
        %Convert pairwise matrix (PCM) into ranking of criteria (RCM) using
        %eignevectors (reference: the analytical hierarchy process, 1990,
        % Thomas L. Saaty
        
        % Note: A simple/fast way to solve for the eigenvectors are:
        % 1. Raise pairwise matrix to powers that are successively squared
        % 2. Sum the rows and normalize summed rows.
        % 3. Stop when the difference between the sums in two consecutive
        %    iterations is smaller than tolerance.
        c=1;
        [m n]=size(M);
        nrM(m,:)=10000; tolmet=0; tolerance=.035;
        while tolmet<1 
            c=c+1;                                        % counter
            M=M^2;                                        % pairwise matrix^2
            sr1M = sum(M,2);                              % sum rows
            sr2M = sum(sr1M);                             % sum of sum rows
            nrM(:,c) = sr1M./sr2M;                        % normalize
            tol(c)=sum(abs(nrM(:,c) - nrM(:,c-1)));       % calc. tolerance
             if tol < tolerance                    % tolerance met?
                tolmet=1;                          % tolerance met, stop iterations
             elseif sum(sum(M))>=10e30 
                 tolmet=1;                         % tolerance unlikely, stop iterations
             end
        end
        disp('Eigenvector of matrix');
        eigvect = nrM(:,end); % eigenvector of PCM
    end

 
end



   
  