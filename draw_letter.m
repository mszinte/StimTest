function draw_letter(const,letter)
% -------------------------------------------------------------------------
% draw_letter(const,letter)
% -------------------------------------------------------------------------
% Goal of the function :
% Draw letter stimulus masks (digital "8"s), test item ("E" or "mirror-E"),
% and distractor items (digital "2"s and "5"s).
% -------------------------------------------------------------------------
% Input(s) :
% const : general settings
% letter : letter stimulus specific settings
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

color = const.colBlack; % black letter stimulus

% compute letter heigth width ratio
letSVal = 2*letter.heigth;  
lw_rat = letSVal/letter.lineW;
letS = 1000;
letW = round(letS/lw_rat);
letB = round(letS/5);

% RGB values (layer 1:3) & Transparency (layer 4)
RGB_ones = ones(letS,letS);
cbar = cat(3,RGB_ones.*color(1),RGB_ones.*color(2),RGB_ones.*color(3));
LetIm(:,:,1:3) = (cbar);
LetIm(:,:,4) = zeros(size(RGB_ones,1),size(RGB_ones,2));
LetIm([1:letW,floor(letS/2-letW/2)+1:floor(letS/2+letW/2),(letS-letW)+1:letS],(letB+1:letS-letB),4) = 255;

% texture "mirror-E" 
LetIm3 = LetIm;
LetIm3([letW+1:floor(letS/2-letW/2)],[letS-letB-letW+1:letS-letB],4) = 255;
LetIm3([floor(letS/2+letW/2)+1:letS-letW],[letS-letB-letW+1:letS-letB],4) = 255;
texLet_3 = Screen('MakeTexture',const.screen,LetIm3);

% texrure "E"
LetImE = LetIm;
LetImE([letW+1:floor(letS/2-letW/2)],[letB+1:letB+letW],4) = 255;
LetImE([floor(letS/2+letW/2)+1:letS-letW],[letB+1:letB+letW],4) = 255;
texLet_E = Screen('MakeTexture',const.screen,LetImE);

% texture "8"
LetIm8 = LetIm;
LetIm8([letW+1:floor(letS/2-letW/2)],[letB+1:letB+letW],4) = 255;
LetIm8([letW+1:floor(letS/2-letW/2)],[letS-letB-letW+1:letS-letB],4) = 255;
LetIm8([floor(letS/2+letW/2)+1:letS-letW],[letB+1:letB+letW],4) = 255;
LetIm8([floor(letS/2+letW/2)+1:letS-letW],[letS-letB-letW+1:letS-letB],4) = 255;
texLet_8 = Screen('MakeTexture',const.screen,LetIm8);

% texture "2"
LetIm2 = LetIm;
LetIm2([letW+1:floor(letS/2-letW/2)],[letS-letB-letW+1:letS-letB],4) = 255;
LetIm2([floor(letS/2+letW/2)+1:letS-letW],[letB+1:letB+letW],4) = 255;
texLet_2 = Screen('MakeTexture',const.screen,LetIm2);

% texture "5"
LetIm5 = LetIm;
LetIm5([letW+1:floor(letS/2-letW/2)],[letB+1:letB+letW],4) = 255;
LetIm5([floor(letS/2+letW/2)+1:letS-letW],[letS-letB-letW+1:letS-letB],4) = 255; 
texLet_5 = Screen('MakeTexture',const.screen,LetIm5);


% Select radom distractor identities
d_ids = [texLet_2,texLet_5]; d_vec = nan(1,size(const.posD,1));
for idx_d = 1:size(const.posD,1)
    rand_idx = randi(2);
    d_vec(idx_d) = d_ids(rand_idx);
end

if letter.testID == 1;      texLet_vecDT = [texLet_E,d_vec]; % vector containing test and distractor textures
elseif letter.testID == 2;  texLet_vecDT = [texLet_3,d_vec];
end
texLet_vecM = repmat(texLet_8,size(texLet_vecDT)); % vector containing mask textures


%% Loop and show stimuli

for nbf = const.start_fr : const.max_fr
    
    % draw test and distractor items
    if nbf >= const.start_test_fr && nbf <= const.end_test_fr
        Screen('DrawTextures',const.screen,texLet_vecDT,[],const.recMat,[]); 
        
    % draw mask items    
    else
        Screen('DrawTextures',const.screen,texLet_vecM,[],const.recMat,[]);
    end
    Screen('Flip',const.screen);
end

% Close textures
 Screen('Close',[texLet_vecM,texLet_vecDT]);

end