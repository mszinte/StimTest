function draw_gabor(const,gabor)
% -------------------------------------------------------------------------
% draw_gabor(const,gabor)
% -------------------------------------------------------------------------
% Goal of the function :
% Draw gabor stimulus masks (white pixel noise), test item (tilted gabor
% patch), and distractor items (vertical gabor patches).
% -------------------------------------------------------------------------
% Input(s) :
% const : general settings
% gabor : gabor stimulus specific settings
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
gabor.contrast          = 1;
gabor.angleDistractor   = 0;
gabor.sigma_period      = 0.9;
gabor.sigma             = gabor.pixpc*gabor.sigma_period;
gabor.bgluminance       = const.colBG(1)/255;
gabor.noisePixSize      = gabor.rad/5; %%%% GOOD PROXY? %%% was 5 before

% Appertures for gabors / noise
apertureGabor = GenerateGaussian(gabor.rad*2, gabor.rad*2, gabor.sigma, gabor.sigma, 0, 0, 0);
apertureGabor = apertureGabor.*255;
apertureNoise = GenerateGaussian(round(gabor.rad*2/gabor.noisePixSize), round(gabor.rad*2/gabor.noisePixSize), gabor.sigma/gabor.noisePixSize, gabor.sigma/gabor.noisePixSize, 0, 0, 0);
apertureNoise = apertureNoise.*255;

% Distractor/test gabor
phase = rand;
grating = GenerateGrating(gabor.rad*2, gabor.rad*2, gabor.angleDistractor, gabor.pixpc, phase, gabor.contrast);
grating = grating + gabor.bgluminance;grating(grating>1)=1;grating(grating<0)=0;grating = grating.*255;
gab3D(:,:,1) = grating;
gab3D(:,:,2) = grating;
gab3D(:,:,3) = grating;
gab3D(:,:,4) = apertureGabor;
texGabor = Screen('MakeTexture',const.screen,gab3D);

% Noise mask (one static noise patch for all locations)
noiseImg = round((rand(round(gabor.rad*2/gabor.noisePixSize), round(gabor.rad*2/gabor.noisePixSize)))*const.colWhite(1));
mask3D(:,:,1) = noiseImg;
mask3D(:,:,2) = noiseImg;
mask3D(:,:,3) = noiseImg;
mask3D(:,:,4) = apertureNoise;
texMask = Screen('MakeTexture',const.screen,mask3D);



%% Loop and show stimuli

% time loop
for nbf = const.start_fr : const.max_fr
    
    texPtr_vec = nan(1,size(const.recMat,2));
    rotAngle_vec = nan(1,size(const.recMat,2));
    
    % position loop
    for tPos = 1:size(const.recMat,2)
        
        % draw distractors / test gabors
        if nbf >= const.start_test_fr && nbf <= const.end_test_fr 
            texPtr_vec(tPos) = texGabor;
            if tPos == 1
                rotAngle_vec(tPos) = gabor.ang;
            else
                rotAngle_vec(tPos) = 0;
            end
        % draw mask (static noise)  
        else
            texPtr_vec(tPos)    = texMask;
            rotAngle_vec(tPos)  = 0;
        end
    end
    
    % Draw all textures
    Screen('DrawTextures',const.screen,texPtr_vec,[],const.recMat,rotAngle_vec);
    Screen('Flip',const.screen);
end

% Close textures
Screen('Close',[texMask,texGabor]);

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