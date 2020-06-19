% This is a function to scramble network, keeping various properties the
% same.
% INPUTS:
% NetworkProp is a structure containing the following fields:
% - N: number of types genes (including 1 for ribosomes)
% - numTFs: number of TFs
% - numInteractions: number of interactions
% ScramblingMethod: scalar representing type of scrambling method:
% if 1: we keep the same number of TFs and number of regulations
% 2: we only keep the number of regulations
% 3: keep in and out degrees the same
% 4: keep # of TFs and # of autoregulatory vs
% non-self interactions.
% 5: keep #TFs,#self vs Non-self, #TF-TF vs
% TF-nonTF interactions

function IntTypeVec = ScrambleNetwork(NetworkProp, ScramblingMethod)

    % Extract parameters
    N = NetworkProp.N;
    numTFs = NetworkProp.numTFs;
    numActInts = NetworkProp.numActInts;
    numInhibInts = NetworkProp.numInhibInts;
    numInteractions = numActInts + numInhibInts;
    % numInteractions = NetworkProp.numInteractions;

    % Main scrambling
    if ScramblingMethod == 1
        % Create another random interaction matrix with the same N, numTFs and 
        % number of interactions.
        TFinds_rand2 = randsample((1:numTFs)',numInteractions,true);
        regulatedinds_rand2 = randsample((1:N)',numInteractions,true);
        IntPairs_rand2 = [TFinds_rand2,regulatedinds_rand2];
        [~,ia] = unique(IntPairs_rand2,'rows');
        while length(ia) < numInteractions
            repeatInd = setdiff(1:numInteractions,ia);
            IntPairs_rand2(repeatInd,:) = [randsample((1:numTFs)',length(repeatInd),true),...
                randsample((1:N)',length(repeatInd),true)];
            [~,ia] = unique(IntPairs_rand2,'rows');
        end
        if length(unique(IntPairs_rand2(:,1)))<numTFs
            warning('drawn number of TFs less than desired!');
        end
        IntInds_rand2 = sub2ind([N,N],IntPairs_rand2(:,1),IntPairs_rand2(:,2));
        IntTypeVec = zeros(N*N,1);
        IntTypeVec(IntInds_rand2(1:numInhibInts)) = -1;
        IntTypeVec(IntInds_rand2(1+numInhibInts:end)) = 1;
        
    elseif ScramblingMethod == 2
        % Create another random interaction matrix with the same N and number of
        % interactions.
        inds = randperm(N*(N-1),numInhibInts+numActInts);
        IntTypeVec = zeros(N*(N-1),1);        
        IntTypeVec(inds(1:numInhibInts)) = -1;
        IntTypeVec(inds(1+numInhibInts:end)) = 1;   
        IntTypeMat = reshape(IntTypeVec,[N-1,N]);
        IntTypeMat(end+1,:) = 0;
        IntTypeVec = IntTypeMat(:);
        
%         IntTypeVec = zeros(N*N,1);
%         inds = randperm(N*N,numInhibInts+numActInts);
%         IntTypeVec(inds(1:numInhibInts)) = -1;
%         IntTypeVec(inds(1+numInhibInts:end)) = 1;   
        
    elseif ScramblingMethod == 3
        % Maintain more specific features of the original network,
        % including:
        % - number of TFs
        % - number of TF-TF interactions vs TF-nonTF interactions
        % - number of self-regulations and total number of interactions
        % - out-degrees of each TF (number of targets each TF has)
        % - in-degree distribution (number of regulators each gene has)
        numSelfActInts = NetworkProp.numSelfActInts;
        numSelfInhibInts = NetworkProp.numSelfInhibInts;
        numSelfInts = numSelfActInts + numSelfInhibInts;
        numNonSelfActInts = numActInts - numSelfActInts;
        numNonSelfInhibInts = numInhibInts - numSelfInhibInts;
        assert(numSelfInts <= numTFs,...
            'number of self interactions must be less than the number of TFs');
        outdeg_TFs = NetworkProp.outdeg_TFs;
        indeg_TFs = NetworkProp.indeg_TFs;
        indeg_nonTFs = NetworkProp.indeg_nonTFs;
        
        TFinds = (1:numTFs)';
%         nonTFinds = (numTFs+1:N-1)';
        selfregInds = randperm(numSelfInts)';
        selfacts = selfregInds(1:numSelfActInts);
        selfinhibs = selfregInds(1+numSelfActInts:end);
        IntPairs_autoreg = [[selfacts;selfinhibs],[selfacts;selfinhibs]];
        IntType_autoreg = [ones(numSelfActInts,1);-ones(numSelfInhibInts,1)];
        
        % shuffle out- and in- degrees
        outdeg_TFs = sort(outdeg_TFs,'descend');
        %outdeg_TFs = outdeg_TFs(randperm(length(outdeg_TFs)));
        indeg_TFs = indeg_TFs(randperm(length(indeg_TFs)));
        indeg_nonTFs = indeg_nonTFs(randperm(length(indeg_nonTFs)));
        indeg_all = [indeg_TFs;indeg_nonTFs];
        assert(sum(outdeg_TFs) == sum(indeg_all),...
            'total number of out and in degrees must be the same!');
        
        % Create interaction pairs by iterating through each TF
        IntPairs_TFinds = repelem(TFinds,outdeg_TFs);
        IntPairs_targetinds = zeros(sum(outdeg_TFs),1);
        currindx = 0;
        remainingdeg = indeg_all;
        for TFind = 1:numTFs
            remainingTargets = find(remainingdeg>0);
            assert(length(remainingTargets)>=outdeg_TFs(TFind),...
                'insufficient number of remaining targets!');
            targets = randsample(remainingTargets,outdeg_TFs(TFind));
            % weights = remainingdeg(remainingTargets)./sum(remainingdeg);
            % targets = randsample(remainingTargets,outdeg_TFs(TFind),false,weights);
            remainingdeg(targets) = remainingdeg(targets)-1;
            IntPairs_targetinds(currindx+1:currindx+outdeg_TFs(TFind)) = targets;
            currindx = currindx+outdeg_TFs(TFind);
        end
        IntPairs = [IntPairs_TFinds,IntPairs_targetinds];
        
        % create interaction pairs
%         IntPairs_targetinds = [repelem(TFinds,indeg_TFs);...
%             repelem(nonTFinds,indeg_nonTFs)];
%         IntPairs_targetinds = IntPairs_targetinds(randperm(length(IntPairs_targetinds)));
%         IntPairs = [IntPairs_TFinds,IntPairs_targetinds];
%         
%         uniqueIntPairs = unique(IntPairs,'rows');
%         numself = sum(IntPairs(:,1)==IntPairs(:,2));
%         while numself > 0 || size(uniqueIntPairs,1)<size(IntPairs,1)
%             IntPairs_targetinds = IntPairs_targetinds(randperm(length(IntPairs_targetinds)));
%             IntPairs = [IntPairs_TFinds,IntPairs_targetinds];
% 
%             uniqueIntPairs = unique(IntPairs,'rows');
%             numself = sum(IntPairs(:,1)==IntPairs(:,2));
%         end

        IntType_nonself = [ones(numNonSelfActInts,1);-ones(numNonSelfInhibInts,1)];
        IntType_nonself = IntType_nonself(randperm(length(IntType_nonself)));
        
        % combine interactions
        IntPairs_all = [IntPairs_autoreg;IntPairs];
        IntType_all = [IntType_autoreg;IntType_nonself];
        IntTypeMat = sparse(IntPairs_all(:,1),IntPairs_all(:,2),IntType_all,N,N);
        IntTypeVec = IntTypeMat(:); 
        
    elseif ScramblingMethod == 4
        % like ScramblingMethod = 1 except that we also take into account
        % the fraction of autoregulation vs non-self interactions (and
        % their respective number of activating and inhibiting
        % interactions).
        numSelfActInts = NetworkProp.numSelfActInts;
        numSelfInhibInts = NetworkProp.numSelfInhibInts;
        numSelfInts = numSelfActInts + numSelfInhibInts;
        numNonSelfActInts = numActInts - numSelfActInts;
        numNonSelfInhibInts = numInhibInts - numSelfInhibInts;
        numNonSelfInts = numNonSelfActInts + numNonSelfInhibInts;
        assert(numSelfInts <= numTFs,...
            'number of self interactions must be less than the number of TFs');

        TFinds = (1:numTFs)';
        selfregInds = randperm(numSelfInts)';
        selfacts = selfregInds(1:numSelfActInts);
        selfinhibs = selfregInds(1+numSelfActInts:end);
        
        % Draw non-self interactions:
        IntPair_TFinds = randsample(TFinds,numNonSelfInts,true);
        IntPair_regulatedinds = randsample((1:N)',numNonSelfInts,true);
        IntPairs = [IntPair_TFinds,IntPair_regulatedinds];
        IntPair_nonself = IntPairs(IntPairs(:,1)~=IntPairs(:,2),:);
        numnonselfdrawn = size(IntPair_nonself,1);
        [~,ia] = unique(IntPair_nonself,'rows');
        while length(ia) < numNonSelfInts
            if numnonselfdrawn < numNonSelfInts
                numtoadd = numNonSelfInts-numnonselfdrawn;
                IntPair_toadd = [randsample((1:numTFs)',numtoadd,true),...
                    randsample((1:N)',numtoadd,true)];
            else
                IntPair_toadd = [];
            end
            repeatInd = setdiff(1:numnonselfdrawn,ia);
            IntPair_nonself(repeatInd,:) = [randsample((1:numTFs)',length(repeatInd),true),...
                randsample((1:N)',length(repeatInd),true)];
            IntPairs = [IntPair_nonself;IntPair_toadd];
            
            IntPair_nonself = IntPairs(IntPairs(:,1)~=IntPairs(:,2),:);
            numnonselfdrawn = size(IntPair_nonself,1);
            [~,ia] = unique(IntPair_nonself,'rows');
        end
        IntInds = sub2ind([N,N],IntPairs(:,1),IntPairs(:,2));
        IntTypeVec = zeros(N*N,1);
        IntTypeVec(IntInds(1:numNonSelfInhibInts)) = -1;
        IntTypeVec(IntInds(1+numNonSelfInhibInts:end)) = 1;
        
        % self-regulation:
        IntTypeVec(sub2ind([N,N],selfacts,selfacts)) = 1;
        IntTypeVec(sub2ind([N,N],selfinhibs,selfinhibs)) = -1;
        
        % check that number of TFs drawn is indeed what we wanted
        numTFs_actual = length(unique([IntPairs(:,1);selfregInds]));
        if numTFs_actual<numTFs
            warning('drawn number of TFs less than desired!');
        end
        
    
    elseif ScramblingMethod == 5
        % like ScramblingMethod = 4 except that we also take into account
        % the fraction of non-self interactions that are TF-otherTF vs 
        % TF-nonTF interactions.
        
        numSelfActInts = NetworkProp.numSelfActInts;
        numSelfInhibInts = NetworkProp.numSelfInhibInts;
        numSelfInts = numSelfActInts + numSelfInhibInts;
        assert(numSelfInts <= numTFs,...
            'number of self interactions must be less than the number of TFs');
        numNonSelfInts = numInteractions - numSelfInts;
        numNonSelf_TFTF_ActInts = NetworkProp.numNonSelf_TFTF_ActInts;
        numNonSelf_TFTF_InhibInts = NetworkProp.numNonSelf_TFTF_InhibInts;
        numNonSelfInts_TFTF = numNonSelf_TFTF_ActInts + numNonSelf_TFTF_InhibInts;
%         numNonSelf_TFnonTF_ActInts = numActInts - numSelfActInts - numNonSelf_TFTF_ActInts;
        numNonSelf_TFnonTF_InhibInts = numInhibInts - numSelfInhibInts - numNonSelf_TFTF_InhibInts;
        numNonSelfInts_TFnonTF = numNonSelfInts - numNonSelfInts_TFTF;
        
        TFinds = (1:numTFs)';
        selfregInds = randperm(numTFs,numSelfInts)';
        selfacts = selfregInds(1:numSelfActInts);
        selfinhibs = selfregInds(1+numSelfActInts:end);
        
        % Draw non-self TF-TF interactions:
        IntPair_TFinds = randsample(TFinds,numNonSelfInts,true);
        % TF-otherTF interactions:
        IntPair_regulatedinds_TF = randsample(TFinds,numNonSelfInts_TFTF,true);
        IntPairs_TFTF = [IntPair_TFinds(1:numNonSelfInts_TFTF),IntPair_regulatedinds_TF];
        IntPair_nonself = IntPairs_TFTF(IntPairs_TFTF(:,1)~=IntPairs_TFTF(:,2),:);
        numnonselfdrawn = size(IntPair_nonself,1);
        [~,ia] = unique(IntPair_nonself,'rows');
        while length(ia) < numNonSelfInts_TFTF
            if numnonselfdrawn < numNonSelfInts_TFTF
                numtoadd = numNonSelfInts_TFTF-numnonselfdrawn;
                IntPair_toadd = [randsample((1:numTFs)',numtoadd,true),...
                    randsample((1:numTFs)',numtoadd,true)];
            else
                IntPair_toadd = [];
            end
            repeatInd = setdiff(1:numnonselfdrawn,ia);
            IntPair_nonself(repeatInd,:) = [randsample((1:numTFs)',length(repeatInd),true),...
                randsample((1:numTFs)',length(repeatInd),true)];
            IntPairs_TFTF = [IntPair_nonself;IntPair_toadd];
            
            IntPair_nonself = IntPairs_TFTF(IntPairs_TFTF(:,1)~=IntPairs_TFTF(:,2),:);
            numnonselfdrawn = size(IntPair_nonself,1);
            [~,ia] = unique(IntPair_nonself,'rows');
        end
        IntInds = sub2ind([N,N],IntPairs_TFTF(:,1),IntPairs_TFTF(:,2));
        IntTypeVec = zeros(N*N,1);
        IntTypeVec(IntInds(1:numNonSelf_TFTF_InhibInts)) = -1;
        IntTypeVec(IntInds(1+numNonSelf_TFTF_InhibInts:end)) = 1;
        % TF-nonTF interactions:
        IntPair_regulatedinds_nonTF = randsample((1+numTFs:N)',numNonSelfInts_TFnonTF,true);
        IntPairs_TFnonTF = [IntPair_TFinds(1+numNonSelfInts_TFTF:end),IntPair_regulatedinds_nonTF];
        [~,ia] = unique(IntPairs_TFnonTF,'rows');
        while length(ia) < numNonSelfInts_TFnonTF
            repeatInd = setdiff(1:numNonSelfInts_TFnonTF,ia);
            IntPairs_TFnonTF(repeatInd,:) = [randsample((1:numTFs)',length(repeatInd),true),...
                randsample((1+numTFs:N)',length(repeatInd),true)];
            [~,ia] = unique(IntPairs_TFnonTF,'rows');
        end
        IntInds = sub2ind([N,N],IntPairs_TFnonTF(:,1),IntPairs_TFnonTF(:,2));
        IntTypeVec(IntInds(1:numNonSelf_TFnonTF_InhibInts)) = -1;
        IntTypeVec(IntInds(1+numNonSelf_TFnonTF_InhibInts:end)) = 1;
        
        % self-regulation:
        IntTypeVec(sub2ind([N,N],selfacts,selfacts)) = 1;
        IntTypeVec(sub2ind([N,N],selfinhibs,selfinhibs)) = -1;
        
        % check that number of TFs drawn is indeed what we wanted
        IntPairs_all = [IntPairs_TFTF;IntPairs_TFnonTF];
        numTFs_actual = length(unique([IntPairs_all(:,1);selfregInds]));
        if numTFs_actual<numTFs
            warning('drawn number of TFs less than desired!');
        end
        
        
    elseif ScramblingMethod == 6 % (relaxation of method 3)
        % Maintain the following features:
        % - number of TFs
        % - number (and type) of TF-TF interactions vs TF-nonTF interactions
        % - number (and type) of self-regulations and total interactions
        % - k of the largest out-degrees 
        % [note that the in-degrees are not kept constant]
        
        numSelfActInts = NetworkProp.numSelfActInts;
        numSelfInhibInts = NetworkProp.numSelfInhibInts;
        numSelfInts = numSelfActInts + numSelfInhibInts;
        assert(numSelfInts <= numTFs,...
            'number of self interactions must be less than the number of TFs');
        numNonSelfActInts = numActInts - numSelfActInts;
        numNonSelfInhibInts = numInhibInts - numSelfInhibInts;
        outdeg_TFs = NetworkProp.outdeg_TFs;
        
        
    end

end