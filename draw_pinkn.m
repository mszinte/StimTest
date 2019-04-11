function draw_pinkn(const,pinkn)
% -------------------------------------------------------------------------
% draw_pinkn(const,pinkn)
% -------------------------------------------------------------------------
% Goal of the function :
% Draw pink noise stimulus masks (unoriented pink noise), test item 
% (orientation filtered pink noise), and distractor items (unoriented pink 
% noise).
% -------------------------------------------------------------------------
% Input(s) :
% const : general settings
% pinkn : letter stimulus specific settings
% -------------------------------------------------------------------------
% Output(s):
% -
% -------------------------------------------------------------------------
% Function created by Nina M. Hanning (hanning.nina@gmail.com) 
% and Martin Szinte (martin.szinte@gmail.com)
% Last update : 10 / 04 / 2019
% Project :     StimtTest
% Version :     1.0
% -------------------------------------------------------------------------

% Stimulus settings
const.noiseRR_t     = 1/60;                                     % noise image refresh rate [frames]
const.noiseRR_fr    = round(const.noiseRR_t/const.frame_dur);   % noise image refresh rate [sec]
const.test_rot      = 40;                                       % tilt angle
const.sizeIm        = const.stimRad*2;
const.sigma         = round(const.sizeIm/3);                    % sigma of the image window

if rem(const.test_dur_fr,const.noiseRR_fr); 
    error('Test duration nees to be multiple of noise refresh rate.');
end

% Make raised cosine mask
imMask = makeRaisedCosineMask(const.sizeIm,const.sizeIm,const.sigma);

% Orientation filter
rSigma2 = 0.15; aSigma = pinkn.filt;

filt_size = 2.^ceil(log2([const.sizeIm,const.sizeIm]));
ImSize_x = filt_size(2); ImSize_y = filt_size(1);

Filter = zeros(ImSize_y,ImSize_x);
for id_X = 1:ImSize_x
    for id_Y = round(ImSize_y/2):ImSize_y
        
        x = id_X - ImSize_x/2; y = id_Y - ImSize_y/2;
        alpha = atan2d(y,x);
        
        r2 = x^2 + y^2; r2 = r2/(ImSize_x*ImSize_y);
        fVal = exp(-(alpha-90)^2/aSigma^2) * exp(-r2/rSigma2);
        
        Filter(id_Y,id_X) = fVal;
        Filter(ImSize_y-id_Y+1,ImSize_x-id_X+1) = fVal;
    end
end
Filter=rot90(Filter);

% Specify when to refresh the noise images
timeRef_vec = zeros(const.start_fr,const.max_fr);
timeRef_vec([const.start_test_fr : -const.noiseRR_fr : const.start_fr, const.start_test_fr+const.noiseRR_fr : +const.noiseRR_fr : const.max_fr]) = 1;


%% Loop and show stimuli

texPtr_vec = [];

% time loop
for nbf = const.start_fr : const.max_fr
    
    % when time to refresh, make new textures
    if nbf == 1 || timeRef_vec(nbf)
        Screen('Close',texPtr_vec);
        
        % position loop
        texPtr_vec = nan(1,size(const.recMat,2));
        rotAngle_vec = zeros(1,size(const.recMat,2));
    
        for tPos = 1:size(const.recMat,2)
            
            % pink noise image
            noiseIm = make_pinkn(const.stimRad);
            
            % if test time filter test noise image
            if nbf >= const.start_test_fr && nbf <= const.end_test_fr && tPos == 1

                meanSub = mean(noiseIm(:));
                noiseIm_fft = fftshift(fft2(noiseIm-meanSub,filt_size(1),filt_size(2)));    % fft and shift

                noiseIm_fft_filt = Filter .* noiseIm_fft;                                   % apply filter
                
                noiseIm_filt = real(ifft2(ifftshift(noiseIm_fft_filt)));                    % shift back
                
                noiseIm_filt = noiseIm_filt(1:size(noiseIm,1),1:size(noiseIm,2));
                noiseIm_filt = noiseIm_filt+meanSub;
                noiseIm_filt = noiseIm_filt-min(noiseIm_filt(:));
                noiseIm_filt = noiseIm_filt./max(noiseIm_filt(:));
                noiseIm      = noiseIm_filt;
                
                rotAngle_vec(tPos) = pinkn.tilt*const.test_rot;                             % test rotation angle 
            end
            
            % put noise in raised cosine mask
            noiseIm_patch(:,:,1) = noiseIm;
            noiseIm_patch(:,:,2) = imMask;

            % make textures
            texPtr_vec(tPos) = Screen('MakeTexture',const.screen,noiseIm_patch.*255);
        end
    end

    % Draw all textures
    Screen('DrawTextures',const.screen,texPtr_vec,[],const.recMat,rotAngle_vec);
    Screen('Flip',const.screen);
end

end


% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function raisedCosMask = makeRaisedCosineMask(imX,imY,nCosSteps,apHeight,apWidth)
% function to make a raised cosine mask for a 2D visual stimulus. 
% [imX,imY]: image size; nCosSteps: step number; [apHeight,apWidth]: aperture size;
% by D Aagten-Murphy (2015), edited by NM Hanning (2017)

if nargin == 3;     apHeight = imX; apWidth = apHeight; % aperture has same diameter as Image
elseif nargin == 4; apWidth = apHeight;                 % aperture is smaller than Image and round
end
HWratio = apHeight/apWidth;

imX = ceil(imX/2)*2;imY = ceil(imY/2)*2; % ensure even size
    
if min([imX imY])/2 < nCosSteps; error('Error cosine mask: sigma too big'); end

[X,Y] = meshgrid(-imX/2+1:imX/2,-imY/2+1:imY/2);
radii = sqrt((X/HWratio).^2 + Y.^2);

% Linear transformation to scale the radii values so that the value 
% corresponding to the inner edge of the ramp is equal to (zero x pi) and 
% the value for the outer edge is equal to (1 x pi). 
% The cosine of these values will be 1.0 and zero, respectively.
 
% set inner edge to zero
radii = radii - radii(round(end/2),round(end/2+apWidth/2-nCosSteps)); 

% Do linear transform to set outer edge to pi
outerVal = radii(round(end/2),round(end/2+(apWidth/2-1)));
radii = radii * pi/outerVal ;

radii((radii<=0)) = 0;  	% set values more central than the soft aperture to 0 (ie, cos(0) = 1)
radii((radii>=pi)) = pi;  	% set values more beyond soft aperture to pi (ie, cos(pi) = 0)

raisedCosMask = .5 + .5 * cos(radii);	% cos of all transformed radial values 
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function [noiseIm] = make_pinkn(imR)
% Create 1/f noise image of size [2*imR by 2*imR]
% by D Aagten-Murphy (2015), edited by NM Hanning (2019)

a = rand(imR*2) ;
F = fftshift(fft2(a) );                 % do FFT and shift origin to centre
[x,y] = meshgrid(-(imR):((imR)-1));     % find radial frequencies
radialFreq = sqrt(x.^2 + y.^2);

radialFreq(end/2+1,end/2+1) = 1;        % make sure centre freq is set to 1.
F = F .* (1./radialFreq) ;              % multiply F with inverse of rad Freq.
noiseIm = real( ifft2(fftshift(F)) );   % shift F back and do Inverse FFT

noiseIm = noiseIm - mean(noiseIm(:)); 	% normalise amplitude range
maxVal = max(abs(noiseIm(:)));
noiseIm = noiseIm/(2*maxVal) + .5;

noiseIm = ((noiseIm-.5).*0.9)+0.5;      % slightly reduce contrast to avoid it clipping by chance
end