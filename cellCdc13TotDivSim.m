function [cellTotPar,matCellGrowth,matCellCdc13]=cellCdc13TotDivSim(growthRate,numCells,simTime,numNoise,synRate,degRate,meanDivLen,meanDivCdc13,synNoise)
%use with adder function

threshTime=0;

guessMean=meanDivLen/2;
matCellGrowth=normrnd(guessMean,abs(0.1*guessMean),[numCells,1]);
matCellCdc13=matCellGrowth*0.1;
matDivCdc13=zeros(numCells,1);

for i=1:simTime;
    i
    indMax=sum(any(~isnan(matCellGrowth),1));
    matCellGrowth=matCellGrowth(:,1:indMax);
    matCellCdc13=matCellCdc13(:,1:indMax);
    
    matGrow=~any(matCellGrowth==0,2);
    %find the last valid point of a cell
    currCellEnd=sum(~isnan(matCellGrowth),2);
    %check which of these last points are valid (is cell still growing)
    valEnds=currCellEnd(matGrow);
    
    %valid growing cell names;
    valGrowCells=find(matGrow);
    %find growth vals that are valid
    indGrowEnd=sub2ind(size(matCellGrowth),valGrowCells,valEnds);
    indNewGrowEnd=sub2ind([size(matCellGrowth,1),size(matCellGrowth,2)+1],valGrowCells,valEnds+1);
    
    %Current growing or dividing ends
    currSizes=matCellGrowth(indGrowEnd);
    currCdc13=matCellCdc13(indGrowEnd);
    
    %get dimensions right
    addZero=zeros(size(matCellGrowth,1)-size(matDivCdc13,1),1);
    matDivCdc13=[matDivCdc13;addZero];
    matDivCdc13(currCellEnd==1)=normrnd(meanDivCdc13,numNoise,sum(currCellEnd==1),1);
    synNoiseMat=normrnd(synRate,synNoise,size(currSizes,1),1);
        if i>500
            matDivCdc13(currCellEnd==1)=normrnd(100000,numNoise,sum(currCellEnd==1),1);
    end

    cDivCdc13=matDivCdc13(matGrow);
    %Include if using a timer
    currpDivEv=cDivCdc13.*(valEnds>threshTime);
    %Generate growth amounts
    currGrow=currSizes+currSizes.*growthRate; %0.006 works well
    %Generate amount of Cdc13 added
    currIncCdc13=currCdc13+currSizes.*synNoiseMat-currCdc13*degRate;
    
    %generate random numbers to choose which cells grow/don't
    
    currDiv=(currCdc13./currSizes)>currpDivEv & currpDivEv~=0;
    %note a division, by naming it 0
    currGrow(currDiv)=0;
    currIncCdc13(currDiv)=0;
    %make new horz cat to deal with current cells:
   
    matCellGrowth=[matCellGrowth,nan(size(matCellGrowth,1),1)];
    matCellGrowth(indNewGrowEnd)=currGrow;
    
     matCellCdc13=[matCellCdc13,nan(size(matCellCdc13,1),1)];
    matCellCdc13(indNewGrowEnd)=currIncCdc13;
    
    %count number of division, mult by 2 to initiate new cells
    numDiv=sum(currDiv);
    
    matNewPar=[];
    if numDiv>0
        valPars=valGrowCells(currDiv);    
        matNewPar=nan(2*numDiv,2);
        %double up cells
        matNewPar(1:2:end,1)=valPars;
        matNewPar(2:2:end,1)=valPars;
        newChildNames=[size(matCellGrowth,1)+1:size(matCellGrowth,1)+numDiv*2];

        %new Children
        matNewPar(1:end,2)=newChildNames';

        %new cell holder
        newCurrCells=nan(numDiv*2,size(matCellGrowth,2));
        newCurrCellsCdc13=newCurrCells;
        %calcBirthSize
        
        newBirthSizes=nan(numDiv*2,1);
        baseSize=(currSizes(currDiv)/2)*1.05;
        assym=normrnd(0,0.1,numDiv,1);
        
        newBirthSizes(1:2:end)=baseSize-assym;
        newBirthSizes(2:2:end)=baseSize+assym;
        newBirthSizes=abs(newBirthSizes);
        newCurrCells(:,1)=newBirthSizes;
        newCurrCellsCdc13(:,1)=newBirthSizes*synRate*10;
        
        matCellGrowth=[matCellGrowth;newCurrCells];
        matCellCdc13=[matCellCdc13;newCurrCellsCdc13];
    end
    matTotPar=[valGrowCells(~currDiv),valGrowCells(~currDiv)];
    matTotPar=[matTotPar;matNewPar];
    
    cellTotPar{i}=matTotPar;
    
end


