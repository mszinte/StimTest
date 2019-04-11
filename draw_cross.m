function draw_cross(const,cross)
% -------------------------------------------------------------------------
% draw_cross(const,cross)
% -------------------------------------------------------------------------
% Goal of the function :
% Draw cross stimulus masks (squares), test item (cross with vertical line 
% offset), and distractor items (symmetrical crosses).
% -------------------------------------------------------------------------
% Input(s) :
% const : general settings
% letter : cross stimulus specific settings
% -------------------------------------------------------------------------
% Output(s):
% -
% -------------------------------------------------------------------------
% Function created by Nina M. Hanning (hanning.nina@gmail.com) 
% and Martin Szinte (martin.szinte@gmail.com)
% Last update : 09 / 04 / 2019
% Project :     StimtTest
% Version :     1.0
% -------------------------------------------------------------------------


%% Prepare stimuli
colCross = const.colBlack; % black cross stimulus


%% Loop and show stimuli

for nbf = const.start_fr : const.max_fr
    
    % draw test and distractor items
    if nbf >= const.start_test_fr && nbf <= const.end_test_fr
        
        
        % distractor crosses
        if ~isempty(const.posD)
        crossMat = const.posD;
        crosslocMatH = nan(2,(size(crossMat,1)*2));
        crosslocMatV = nan(2,(size(crossMat,1)*2));
        
        for i = 1:size(crossMat,1)
            crosslocMatH(1,(i*2-1)) = crossMat(i,1)-const.stimRad; 	crosslocMatV(1,(i*2-1)) = crossMat(i,1);
            crosslocMatH(1,(i*2))   = crossMat(i,1)+const.stimRad; 	crosslocMatV(1,(i*2))   = crossMat(i,1);
            crosslocMatH(2,(i*2-1)) = crossMat(i,2);                crosslocMatV(2,(i*2-1)) = crossMat(i,2)-const.stimRad;
            crosslocMatH(2,(i*2))   = crossMat(i,2);                crosslocMatV(2,(i*2))   = crossMat(i,2)+const.stimRad;
        end
        crosslocMat = [crosslocMatH,crosslocMatV];
        
        Screen('DrawLines',const.screen,crosslocMat,cross.lineW,colCross,[],0);
        end
        
        % test cross 
        crossMat = const.posT;
        crosslocMatH = nan(2,(size(crossMat,1)*2));
        crosslocMatV = nan(2,(size(crossMat,1)*2));
        
        for i = 1:size(crossMat,1)
            crosslocMatH(1,(i*2-1)) = crossMat(i,1)-const.stimRad; 	crosslocMatV(1,(i*2-1)) = crossMat(i,1)+cross.offset;
            crosslocMatH(1,(i*2))   = crossMat(i,1)+const.stimRad; 	crosslocMatV(1,(i*2))   = crossMat(i,1)+cross.offset;
            crosslocMatH(2,(i*2-1)) = crossMat(i,2);                crosslocMatV(2,(i*2-1)) = crossMat(i,2)-const.stimRad;
            crosslocMatH(2,(i*2))   = crossMat(i,2);                crosslocMatV(2,(i*2))   = crossMat(i,2)+const.stimRad;
        end
        crosslocMat = [crosslocMatH,crosslocMatV];
        
        Screen('DrawLines',const.screen,crosslocMat,cross.lineW,colCross,[],0);

        
    % draw mask items
    else
        Screen('FrameRect',const.screen,colCross,const.recMat,cross.lineW);
        
    end
    
    Screen('Flip',const.screen);
end

end