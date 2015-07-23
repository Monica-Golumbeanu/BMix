% 
% Infers the model parameters based on maximum likelihood
% Input: DataPathArray - the path of the file containing TC mismatch
%                        frequency profile
%        NegDataPathArray - path of the file containing the mismatch
%                           frequency for the non TC mismatches (AC and GC)
%        OutDir - path of the directory where the output will be saved
%        CovTHR - minimum coverage to consider
%        PosteriorTHR - posterior threshold used for classification
%
% Written by Pejman Mohammadi and Monica Golumbeanu
% pejman.mohammadi@bsse.ethz.ch
% monica.golumbeanu@bsse.ethz.ch
%
%
function learn_model(DataPathArray, NegDataPathArray, OutDir, CovTHR, PosteriorTHR)

try
    if nargin < 1
        disp('Incorrect number of arguments! \n Usage: learn_model(data,neg_data,out_dir,coverage_threshold');
    else
        %% load positive data (T->C)
        cX=[];cY=[];cW=[];
        for fi = 1: length(DataPathArray)
            DataPath = DataPathArray{fi};
            Data = sRead_Table(DataPath);
            Data = sStruct_RowDel(Data, Data.Y<CovTHR);
            [cX, cY, cW] = Compress_Data(Data.X,Data.Y,cX, cY, cW);         % compress duplicates
        end
        %% load negative data (A->C, G->C)
        cXn=[];cYn=[];cWn=[];
        for fi = 1: length(NegDataPathArray)
            DataPathNeg = NegDataPathArray{fi};
            DataN = sRead_Table(DataPathNeg);
            DataN = sStruct_RowDel(DataN, DataN.Y<CovTHR);
            [cXn, cYn, cWn] = Compress_Data(DataN.X,DataN.Y,cXn, cYn, cWn); % compress duplicates
        end
        
        %% Create output directory
        if ~exist(OutDir,'dir')
            mkdir(OutDir);
        end

        if ~exist([OutDir '/Figures/'],'dir')
            mkdir([OutDir '/Figures/']);
        end
        
        %% Estimate the parameters using gradient descent
        tic
        try
            % New Matlab versions
            OPToptions = optimoptions(@fmincon, 'MaxFunEvals', 3000, 'MaxIter', 1000, 'GradObj', 'on', 'TolX', 1E-5 , 'TolFun', 1E-6, 'Display', 'off');
        catch
            % Older Matlab versions
            OPToptions = optimset('MaxFunEvals', 3000, 'MaxIter', 1000, 'GradObj', 'on', 'TolX', 1E-5 , 'TolFun', 1E-6, 'Display', 'off');
        end
        Problem2.objective = @(Z2)Get_total_Like(cX,cY,cW,cXn,cYn,cWn, Z2);
        Problem2.x0 = [.1; .1;  .1; .1];
        Problem2.solver = 'fmincon';
        Problem2.options= OPToptions;
        Problem2.lb = [ 1E-3 1E-3 1E-6 1E-3];
        Problem2.ub = [0.25 1 1 1];  % .499
        [OPTtheta,~,exitflag] = fmincon(Problem2);
        toc
        
        %% Print the estimated parameters
        file_param  = fopen([OutDir,'parameters.txt'],'a');
        MLE_epsilon = OPTtheta(1);
        MLE_gamma   = OPTtheta(2);
        MLE_nu      = OPTtheta(3);
        MLE_psi     = OPTtheta(4);
        fprintf('The learned parameters are:\n epsilon = %f\n gamma = %f\n nu = %f\n psi = %f\n', MLE_epsilon, MLE_gamma, MLE_nu,MLE_psi);
        fprintf(file_param,'The learned parameters are:\n epsilon = %f\n gamma = %f\n nu = %f\n psi = %f\n', MLE_epsilon, MLE_gamma, MLE_nu,MLE_psi);
        fclose(file_param);

        %% Classify the genomic loci using the posterior and the inferred parameters
        [cMLE_ClassI, cMLE_pC] = sClassify(cX,cY, MLE_epsilon, MLE_gamma, MLE_nu, MLE_psi);
        
        %% Save the results to file
        cXYall = [cX, cY];
        for fi = 1: length(DataPathArray)
            DataPath = DataPathArray{fi};
            [~,FileName,~] = fileparts(DataPath);
            Data = sRead_Table(DataPath);
            Data = sStruct_RowDel(Data, Data.Y<CovTHR);
            [Cxy, ~, iC] = unique([Data.X Data.Y],'rows');
            [~, iCxy, iCxyall] = intersect(Cxy, cXYall, 'rows');
            Cxy_MLE_ClassI(iCxy) = cMLE_ClassI(iCxyall);
            Cxy_cMLE_pC(iCxy)    = cMLE_pC(iCxyall);
            Data.MLE_ClassI = Cxy_MLE_ClassI(iC);
            Data.MLE_pC     = Cxy_cMLE_pC(iC);
            sWrite_Table([OutDir '/' FileName '.results'], Data)
        end
        
        %% Plot results
        Fig2 = figure;
        F = (cMLE_ClassI==1);
        semilogx(cY(F),cX(F)./cY(F)*100, '.', 'Displayname', ' SNP', 'color' , [.9 .3 .1], 'markersize', 10);
        hold on
        F = (cMLE_ClassI==2);
        semilogx(cY(F),cX(F)./cY(F)*100, '.r', 'Displayname', '~RNP & ~SNP', 'markersize', 10,'color', [.7 .05 .1]);
        F = (cMLE_ClassI==3);
        semilogx(cY(F),cX(F)./cY(F)*100, '.k', 'Displayname', ' RNP & ~SNP', 'markersize', 10, 'color',[.05 .05 .2]);
        F = (cMLE_pC<PosteriorTHR);
        semilogx(cY(F),cX(F)./cY(F)*100, '.', 'Displayname', ' RNP & ~SNP', 'markersize', 10, 'color',[0.7 0.7 0.7]);
        legend ('variant','background','cross-link loci','not classified');
        set(gcf, 'position', [1 1 350 350])
        set(gca, 'YTick', 0:25:100);
        grid on
        ylim([-1 101])
        sSavePlot(Fig2 ,[OutDir '/Figures/' 'Classification_' FileName]);
    end
catch Er
    disp(Er.message)
    disp(Er.getReport)
end
end


% calculate the total data likelihood
function [nLL, GnLL] = Get_total_Like(X,Y,W, Xn,Yn,Wn, Z)
[nLL1, GnLL1] =  sGet_like(X,Y,W, Z(1), Z(2), Z(3), Z(4));
if ~isempty(Xn)
    [nLL2, GnLL2] =  sGet_likeN(Xn,Yn,Wn, Z(1), Z(3));
else
    nLL2 = 0;
    GnLL2 = zeros(size(GnLL1));
end
nLL  =  nLL1 + nLL2;
GnLL = GnLL1 + GnLL2;
end

function [cX, cY, cW] = Compress_Data(X,Y, cX0, cY0, cW0)
% this function gets a pair of X and Y then comresses it so that
% [cX0, cY0] is assumed to have unique rows and cW0 stores their frequency.

[Cxy, ~, ic] = unique([X,Y], 'rows');
CxyW = hist(ic, 1:size(Cxy,1))';

if nargin==2 || isempty(cX0)
    cX=Cxy(:,1);
    cY=Cxy(:,2);
    cW=CxyW;
else
    Cxy0 = [cX0, cY0]; clear cX0 cY0;
    [~,ia,ib] = intersect(Cxy0,Cxy,'rows');
    if(~isempty(ia) && ~isempty(ib))
        cW0(ia) = cW0(ia)+CxyW(ib);
    end
    Fb = true(size(CxyW)); Fb(ib) = false;
    cX = [Cxy0(:,1); Cxy(Fb,1)];
    cY = [Cxy0(:,2); Cxy(Fb,2)];
    cW = [cW0; CxyW(Fb)];
end
end
