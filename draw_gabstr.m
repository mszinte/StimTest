function draw_gabstr(const,gabstr)
% -------------------------------------------------------------------------
% draw_gabstr(const,gabstr)
% -------------------------------------------------------------------------
% Goal of the function :
% Draw gabor stream stimulus masks (alternating white pixel noise and 
% vertical gabor patches), test item (tilted gabor patch), and distractor 
% items (vertical gabor patches).
% -------------------------------------------------------------------------
% Input(s) :
% const : general settings
% gaborstr : gabor stream stimulus specific settings
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

% General settings
gabstr.contrast          = 1;
gabstr.angleDistractor   = 0;
gabstr.sigma_period      = 0.9;
gabstr.sigma             = gabstr.pixpc*gabstr.sigma_period;
gabstr.bgluminance       = const.colBG(1)/255;
gabstr.noisePixSize      = gabstr.rad/5; %%%% GOOD PROXY? %%% was 5 before

% Appertures for gabors / noise
apertureGabor = GenerateGaussian(gabstr.rad*2, gabstr.rad*2, gabstr.sigma, gabstr.sigma, 0, 0, 0);
apertureGabor = apertureGabor.*255;
apertureNoise = GenerateGaussian(round(gabstr.rad*2/gabstr.noisePixSize), round(gabstr.rad*2/gabstr.noisePixSize), gabstr.sigma/gabstr.noisePixSize, gabstr.sigma/gabstr.noisePixSize, 0, 0, 0);
apertureNoise = apertureNoise.*255;

% Grating zero contrast / Blank
grating0 = GenerateGrating(gabstr.rad*2, gabstr.rad*2, gabstr.angleDistractor, gabstr.pixpc, 0, 0);
grating0 = grating0 + gabstr.bgluminance;grating0(grating0>1)=1;grating0(grating0<0)=0;grating0 = grating0.*255;
blank3D(:,:,1) = grating0;
blank3D(:,:,2) = grating0;
blank3D(:,:,3) = grating0;
blank3D(:,:,4) = apertureGabor;
texBlank    = Screen('MakeTexture',const.screen,blank3D);


%% Specify for each frame wheter to draw noise masks (0), gabor masks (1), test&distractors (2) or blank (3)
timeNGB_vec = zeros(const.start_fr,const.max_fr);

% test-time: test and distractors
timeNGB_vec(const.start_test_fr:const.end_test_fr) = 2;

% pre-test time: alternate noise and gabors
gab_vec =  const.start_test_fr - const.test_dur_fr*2 : -const.test_dur_fr*2 : const.start_fr;
for idx_gab = 1:numel(gab_vec)
    timeNGB_vec(gab_vec(idx_gab):gab_vec(idx_gab)+const.test_dur_fr-1) = 1;
end

% post-test time: alternate noise and blanks
blank_vec =  const.end_test_fr + const.test_dur_fr+1 : const.test_dur_fr*2 : const.max_fr;
for idx_gab = 1:numel(blank_vec)
    timeNGB_vec(blank_vec(idx_gab):blank_vec(idx_gab)+const.test_dur_fr-1) = 3;
end


%% Loop and show stimuli

% time loop
for nbf = const.start_fr : const.max_fr
    
    % Change noise every step of the noise stream
    if nbf == 1 || (timeNGB_vec(nbf)==0 && timeNGB_vec(nbf-1)~=0)
        noiseImg=round((rand(round(gabstr.rad*2/gabstr.noisePixSize),round(gabstr.rad*2/gabstr.noisePixSize)))*const.colWhite(1));
        mask3D(:,:,1) = noiseImg; 
        mask3D(:,:,2) = noiseImg; 
        mask3D(:,:,3) = noiseImg;
        mask3D(:,:,4) = apertureNoise;
    end
    
    % Change phase every step of the stim stream
    if nbf == 1 || (timeNGB_vec(nbf)==1 && timeNGB_vec(nbf-1)~=1)
        phase = rand;grating = GenerateGrating(gabstr.rad*2,gabstr.rad*2,gabstr.angleDistractor,gabstr.pixpc, phase, gabstr.contrast);
        grating = grating + gabstr.bgluminance;grating(grating>1)=1;grating(grating<0)=0;grating = grating.*255;
        gab3D(:,:,1) = grating;
        gab3D(:,:,2) = grating;
        gab3D(:,:,3) = grating;
        gab3D(:,:,4) = apertureGabor;
    end
    texMask = Screen('MakeTexture',const.screen,mask3D);
    texGabor = Screen('MakeTexture',const.screen,gab3D);
    
    rotAngle_vec = nan(1,size(const.recMat,2));
    texPtr_vec = nan(1,size(const.recMat,2));
    
    % position loop
    for tPos = 1:size(const.recMat,2)

        if timeNGB_vec(nbf) == 0                    % draw noise mask
            texPtr_vec(tPos) = texMask;
            rotAngle_vec(tPos) = 0;
        elseif timeNGB_vec(nbf) == 1 
            texPtr_vec(tPos) = texGabor;            % draw gabor mask
            rotAngle_vec(tPos) = 0;
        elseif timeNGB_vec(nbf) == 2    
            if tPos == 1                          	% draw test gabor
                texPtr_vec(tPos) = texGabor;
                rotAngle_vec(tPos) = gabstr.ang;
            else
                texPtr_vec(tPos) = texGabor;        % draw distractor gabor
                rotAngle_vec(tPos) = 0;
            end
         elseif timeNGB_vec(nbf) == 3               % draw blank
            texPtr_vec(tPos) = texBlank;
            rotAngle_vec(tPos) = 0;              
        end
        
    end
    
    % Draw all textures, then close them
    Screen('DrawTextures',const.screen,texPtr_vec,[],const.recMat,rotAngle_vec);
    Screen('Close',[texGabor,texMask]);
    
    Screen('Flip',const.screen);
    
end

% Close open textures
Screen('Close',texBlank);

end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function gaussian = GenerateGaussian(m, n, sigma1,sigma2, theta,xoff,yoff)
% generates a 2D gaussian with a given patch size, sd, orientation and offset
% gaussian = GenerateGaussian(256,256, 32,32, pi/2,0,0); ishow(gaussian)
% J Greenwood 2009

m = round(m);
n = round(n);

[X,Y]=meshgrid(-m/2:m/2-1,-n/2:n/2-1);
X  = X-xoff;
Y  = Y-yoff;

% speedy variables
c1=cos(pi-theta);
s1=sin(pi-theta);
sig1_squared=2*sigma1*sigma1;
sig2_squared=2*sigma2*sigma2;

% co-ordinate matrices
Xt=zeros(m,n);
Yt=zeros(m,n);

% rotate co-ordinates
Xt = X.*c1 + Y.*s1;
Yt = Y.*c1 - X.*s1;

gaussian = exp(-(Xt.*Xt)/sig1_squared-(Yt.*Yt)/sig2_squared);
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function grating = GenerateGrating(m,n,angle,lambda,ph,con)
% grating = GenerateGrating(256, 256, pi/2, 16, pi/2, 1);
% m/n = x/y patch size in pixels; angle = orientation of grating in degrees; 
% lambda = spatial period in pixels; p=phase (0 to 1); con=contrast;
% J Greenwood 2009 & D Jonikaitis 2011

m = round(m);
n = round(n);

[X,Y] = meshgrid(-m/2:m/2-1,-n/2:n/2-1); % Grid size
angle=-angle+90; % Make grating vertical (instead of horizontal)
theta=angle*pi/180; % Orientation in radians
phase=ph*2*pi; % Phase converted to radians

% rotate co-ordinates
Xt2 = X.*(cos(pi/2-theta)) + Y.*(sin(pi/2-theta));
grating = (0.5*con)*cos(Xt2.*((2.0*pi)/lambda)+phase); % use 0.5*contrast to set max and min values around zero
end