%% Stimulus launcher
% -------------------------------------------------------------------------
% Last edit : 11/04/19
% -------------------------------------------------------------------------
% Description: 
% This script launches one out of six stimuli to measure visuospatial
% attention. 
% -------------------------------------------------------------------------
% corresponding publication: 
% Hanning N. M., Deubel H., & Szinte M. 
% Sensitivity measures of visuospatial attention.
% -------------------------------------------------------------------------
close all;

%% Stimulus type
% Select the stimulus type
stim_cond = input(sprintf('\n\t(1) Digital letters\n\t(2) Static gabor patches\n\t(3) Crosses\n\t(4) Pink noise\n\t(5) RDK\n\t(6) Dynamic gabor streams:\n\n\t>>> Choose stimulus: '));
% 1: Digital letters      	% 4: Pink noise
% 2: Static gabor patches 	% 5: Random dot kinematograms
% 3: Crosses                % 6: Dynamic gabor streams

%% General settings
% Ajust screen size and specify item positions and trial timing

% Spatial configurations                
[scrX,scrY] = Screen('WindowSize',0);
spaceDT = 300;                              % space between stimuli                      
const.posT = [scrX/2,scrY/2-spaceDT/2];     % [x,y] coordinates of test item (in pixels)          
const.posD = [scrX/2-spaceDT/2,scrY/2;...   % [x1,y1; x2,y2; ... xN,yN] coordinates of N distractor item(s)
              scrX/2+spaceDT/2,scrY/2;...
              scrX/2,scrY/2+spaceDT/2];     

% Temporal configurations 
trial_dur_t = 2.000;                        % trial duration (in seconds)
test_onset_t = 0.400;                       % test onset time relative to trial start (in seconds)

const.frame_dur = 1/60;                     % frame duration in seconds (e.g. 1/60 for screen refresh rate of 60 Hz)

%% Stimulus-specific parameters
% Specify item size, test identity and task difficulty

switch stim_cond
    case 1 
        % Digital letters
        letter.heigth   = 35;               % item heigth (width is is half)
        letter.lineW    = 4;                % item line width
        letter.testID   = 1;                % test item "E" (1) or "mirror-E" (2)
        test_dur_t      = 5/60;             % test presentation duration (seconds)
        const.stimRad   = letter.heigth/2;  % item size
        
    case 2
        % Static gabor patches
        gabor.rad       = 30;               % item radius 
        gabor.pixpc     = 10;               % gabor frequency (in pixels per cycle; e.g. 10 for 1 cicle per 10 pixel)
        gabor.ang       = -10;              % test item tilt angle (relative to horizontal; e.g. -10 for 10? left tilt)
        test_dur_t      = 0.025;            % test presentation duration (seconds)
        const.stimRad   = gabor.rad;        % item size
        
    case 3 
        % Crosses
        cross.rad       = 20;               % item radius
        cross.lineW     = 4;                % item line width
        cross.offset    = -5;               % test item vertical bar offset (relative to midline; e.g. -15 for 5 pixel offset to the left)
        test_dur_t      = 0.1;              % test presentation duration (seconds)
        const.stimRad   = cross.rad;        % item size 
        
    case 4 
        % Pink noise
        pinkn.rad       = 40;               % item radius
        pinkn.tilt      = 1;                % test item left (-1) or right (1) tilt
        pinkn.filt      = 20;               % test item filter strength (the smaller, the clearer)
        test_dur_t      = 5/60;             % test presentation duration (seconds)
        const.stimRad   = pinkn.rad;        % item size
                
    case 5 
        % Random dot kinematograms
        rdk.rad         = 60;               % item radius
        rdk.dirSignal   = 0;                % test item motion direction (0:right, 90:up, 180:left, 270:down)
        rdk.kappa       = 10;               % test item motion coherence
        test_dur_t      = 0.1;              % test presentation duration (seconds)
        const.stimRad   = rdk.rad;          % item size
        
    case 6 
        % Dynamic gabor streams
        gabstr.rad      = 30;               % gabor patch radius
        gabstr.pixpc    = 10;               % gabor frequency (in pixels per cycle; e.g. 10 for 1 cicle per 10 pixel)
        gabstr.ang      = -10;              % tilt angle relative to horizontal (e.g. -10 for 10? left tilt)
        test_dur_t      = 0.0250;           % test presentation duration (seconds)
        const.stimRad = gabstr.rad;         % item size
end
  
%% Time and locations conversion

% Trial timing (in frames)
const.start_fr       = 1;                                                           % trial start
const.start_test_fr  = round(test_onset_t/const.frame_dur);                         % test presentation start
const.test_dur_fr    = round(test_dur_t/const.frame_dur);                           % test presentation duration
const.end_test_fr    = round(test_onset_t/const.frame_dur)+const.test_dur_fr-1;     % test presentation end
const.max_fr         = round(trial_dur_t/const.frame_dur);                         	% trial end

% Stimulus locations
allPos = [const.posT;const.posD];
const.recMat = nan(4,size(allPos,1));
for tPos = 1:size(allPos,1)
    const.recMat(1,tPos) = allPos(tPos,1) - const.stimRad;
    const.recMat(2,tPos) = allPos(tPos,2) - const.stimRad;
    const.recMat(3,tPos) = allPos(tPos,1) + const.stimRad;
    const.recMat(4,tPos) = allPos(tPos,2) + const.stimRad;
end


%% Open screen and draw stimuli
const.colWhite = [255, 255, 255];
const.colBlack = [  0,   0,   0];
const.colBG = mean([const.colWhite,const.colBlack]);

Screen('Preference','SyncTestSettings', 0.01, 50, 0.25);
Screen('Preference','SuppressAllWarnings', 1);
Screen('Preference','Verbosity', 0);
Screen('Preference','SkipSyncTests',1);
[const.screen,~]=Screen('OpenWindow',0,const.colBG,[0,0,scrX,scrY]);
Screen('BlendFunction',const.screen,GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

DrawFormattedText(const.screen,'Press a button to start.','center','center',[0,0,0]);
Screen('Flip',const.screen);
KbWait;

% Draw specified stimuli
switch stim_cond
    case 1; draw_letter(const,letter);   
    case 2; draw_gabor(const,gabor);
    case 3; draw_cross(const,cross);
    case 4; draw_pinkn(const,pinkn);
    case 5; draw_rdk(const,rdk);
    case 6; draw_gabstr(const,gabstr);
end

Screen('FillRect',const.screen,const.colBG);
Screen('Flip',const.screen);

% Close the screen
WaitSecs(2);sca;
close all ;