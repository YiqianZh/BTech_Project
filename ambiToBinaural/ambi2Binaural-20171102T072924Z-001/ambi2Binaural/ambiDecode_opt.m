clc;clear;

% user is facing along Y-axis

%% polar co-ordinates
% el_p = -90:90     wrt X-Y plane
% az_p = -180:180   projection on X-Y plane wrt Y-axis.

%% binaural co-ordinates
% el_b = -90:270    projection on Y-Z plane wrt Y axis
% az_b = -90:90     wrt Y-Z plane 

% The azimuth and elevation are in binaural polar co-ords hereon i.e.,
% x_b = sin(az_b);
% y_b = cos(az_b)*cos(el_b);
% z_b = cos(az_b)*sin(el_b);

az_range = [-80 -65 -55 -45:5:45 55 65 80];
el_range = -45 + 5.625*(0:49);
numSources = length(az_range) * length(el_range);

% Speaker Configuration
src_az = repmat(az_range, 1, length(el_range));
src_el = repmat(el_range, length(az_range), 1);
src_el = src_el(:)';

% If p is the vector denoting input to the sources; 1st order ambisonics B format is given as:
% B = [W X Y Z]' = Cp
% where,
% C = [... [0.707; x_b; y_b; z_b]_p ...]
%
% Then p is given as:
% p = pinv(C)*B

ambOrder = 1;   % ambisonics encoding order
numChannels = (ambOrder+1)^2;
C = zeros(numChannels, numSources);

% channel ordering
Wch = 1;
Xch = 2;
Ych = 3;
Zch = 4;

C(Wch,:) = ones(1,numSources)/sqrt(2);
C(Xch,:) = sind(src_az);
C(Ych,:) = cosd(src_az).*cosd(src_el);
C(Zch,:) = cosd(src_az).*sind(src_el);

[B, Fs, nBits] = wavread('Anderson-Nearfield.wav');
B = B';

% Fs = 48000;
% s = (wavread('input_mono.wav'))';

% Ambisonics encoding

%% Steady source
% s = s(1,1:10*Fs);
% B = zeros(numChannels, length(s));
% az = -45;
% el = 0;
% B(Wch,:) = s/sqrt(2);
% B(Xch,:) = s*sind(az);
% B(Ych,:) = s*cosd(az)*cosd(el);
% B(Zch,:) = s*cosd(az)*sind(el);

%% Movement along azimuth
% s = s(1,1:length(az_range)*Fs);
% B = zeros(numChannels, length(s));
% el = 0;
% for i = 1:length(az_range)
%     az = az_range(i);
%     B(Wch,(i-1)*Fs+1:i*Fs) = s(1,(i-1)*Fs+1:i*Fs)/sqrt(2);
%     B(Xch,(i-1)*Fs+1:i*Fs) = s(1,(i-1)*Fs+1:i*Fs)*sind(az);
%     B(Ych,(i-1)*Fs+1:i*Fs) = s(1,(i-1)*Fs+1:i*Fs)*cosd(az)*cosd(el);
%     B(Zch,(i-1)*Fs+1:i*Fs) = s(1,(i-1)*Fs+1:i*Fs)*cosd(az)*sind(el);
% end

%% Movement along elevation
% s = s(1,1:length(el_range)*Fs/2);
% B = zeros(numChannels, length(s));
% az = -45;
% for i = 1:length(el_range)/2
%     el = el_range(2*i);
%     B(Wch,(i-1)*Fs+1:i*Fs) = s(1,(i-1)*Fs+1:i*Fs)/sqrt(2);
%     B(Xch,(i-1)*Fs+1:i*Fs) = s(1,(i-1)*Fs+1:i*Fs)*sind(az);
%     B(Ych,(i-1)*Fs+1:i*Fs) = s(1,(i-1)*Fs+1:i*Fs)*cosd(az)*cosd(el);
%     B(Zch,(i-1)*Fs+1:i*Fs) = s(1,(i-1)*Fs+1:i*Fs)*cosd(az)*sind(el);
% end

% Prepare the HRIR data
load('..\CIPIC_hrtf_database\standard_hrir_database\subject_147\hrir_final.mat');    % 147
szHrir = size(hrir_l);
if (szHrir(1) ~= length(az_range) || szHrir(2) ~= length(el_range))
    error('The source grid and hrir dimensions are not consistent.');
end

hrir_l = reshape(hrir_l, numSources, szHrir(3));
hrir_r = reshape(hrir_r, numSources, szHrir(3));

lowCut = 20;
highCut = 100;
hrir_l = hrir_l(:,lowCut:highCut);
hrir_r = hrir_r(:,lowCut:highCut);
% plot(hrir_l');


% Decoder for the given speaker configuration
D = pinv(C);

D_l = hrir_l'*D;
D_r = hrir_r'*D;

% rotate the sound field
% rot = zeros(4,4);
% rot(Wch,Wch) = 1;
% rot(Xch,Zch) = -1;
% rot(Ych,Ych) = 1;
% rot(Zch,Xch) = 1;
% B = rot*B;

y_l = zeros(1, size(B,2));
y_r = zeros(1, size(B,2));
for ch = 1:numChannels
    y_l = y_l + filter(D_l(:,ch),1,B(ch, :));
    y_r = y_r + filter(D_r(:,ch),1,B(ch, :));
end

y = [y_l; y_r];
sound(y',Fs);
wavwrite(y', Fs, 32, 'binauralDecoded.wav');



