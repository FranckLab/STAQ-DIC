% ==============================================
% function funSmoothDispCrack
% ==============================================
function FLocal = funSmoothStrain_adapt(FLocal,coordinatesFEM,elementsFEM,winstepsize,LevelNo,FilterSizeInput,FilterStd)

% close all; Plotdisp_show(ULocal,elementsFEM,coordinatesFEM);

% prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
DoYouWantToSmoothOnceMore = 0; % DoYouWantToSmoothOnceMore = input(prompt);
if DoYouWantToSmoothOnceMore == 0  
    if isempty(FilterStd) == 1
        prompt = 'Choose filter standard deviation(0-default): ';
        FilterStd = input(prompt);
        if FilterStd == 0
            FilterStd = 0.5; 
        end
    else
        if FilterStd == 0
            FilterStd = 0.5;
        end
    end
    if isempty(FilterSizeInput) == 1
        prompt = 'Choose Gaussian filter size(0-default): ';
        FilterSizeInput = input(prompt);
        if FilterSizeInput == 0
            FilterSizeInput = 2*ceil(2*FilterStd)+1; 
        end
    else
        if FilterSizeInput == 0
            FilterSizeInput = 2*ceil(2*FilterStd)+1;
        end
    end
end

SmoothTimes = 1;
while (DoYouWantToSmoothOnceMore==0)
    Coordxnodes = [min(coordinatesFEM(:,1)):winstepsize/(2^(LevelNo-1)):max(coordinatesFEM(:,1))]'; 
    Coordynodes = [min(coordinatesFEM(:,2)):winstepsize/(2^(LevelNo-1)):max(coordinatesFEM(:,2))]';
    Iblur_11 = gridfit(coordinatesFEM(:,1), coordinatesFEM(:,2), FLocal(1:4:end),Coordxnodes,Coordynodes,'regularizer','springs'); 
    Iblur_11=Iblur_11';
    Iblur_22 = gridfit(coordinatesFEM(:,1), coordinatesFEM(:,2), FLocal(4:4:end),Coordxnodes,Coordynodes,'regularizer','springs');  
    Iblur_22=Iblur_22';
    Iblur_21 = gridfit(coordinatesFEM(:,1), coordinatesFEM(:,2), FLocal(2:4:end),Coordxnodes,Coordynodes,'regularizer','springs'); 
    Iblur_21=Iblur_21';
    Iblur_12 = gridfit(coordinatesFEM(:,1), coordinatesFEM(:,2), FLocal(3:4:end),Coordxnodes,Coordynodes,'regularizer','springs');
    Iblur_12=Iblur_12';
    % -------------------------------------------------------
    imageFilter=fspecial('gaussian',FilterSizeInput,FilterStd);
    Iblur_1 = nanconv(Iblur_11,imageFilter,'edge','nanout');
    Iblur_4 = nanconv(Iblur_22,imageFilter,'edge','nanout');
    Iblur_2 = nanconv(Iblur_21,imageFilter,'edge','nanout');
    Iblur_3 = nanconv(Iblur_12,imageFilter,'edge','nanout');
    
    for tempi = 1:size(coordinatesFEM,1)
        [row1,~] = find(Coordxnodes==coordinatesFEM(tempi,1));
        [row2,~] = find(Coordynodes==coordinatesFEM(tempi,2));
        FLocal(4*tempi-3) = Iblur_1(row1,row2);
        FLocal(4*tempi)   = Iblur_4(row1,row2);
        FLocal(4*tempi-2) = Iblur_2(row1,row2);
        FLocal(4*tempi-1) = Iblur_3(row1,row2);
    end
     
    % prompt = 'Do you want to smooth displacement once more? (0-yes; 1-no)';
    % DoYouWantToSmoothOnceMore = input(prompt); 
    SmoothTimes = SmoothTimes+1;
    if SmoothTimes > 2
        DoYouWantToSmoothOnceMore = 1;
    end
    
end