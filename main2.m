
% Check if the given wav file exists:
wavFileName = 'F:\Kuliah\Bismillah\matlab\data\Tes1.wav';

fp = fopen(wavFileName, 'rb');
if (fp<0)
	fprintf('The file %s has not been found!\n', wavFileName);
	return;
end 
fclose(fp);

% Check if .wav extension exists:
if  (strcmpi(wavFileName(end-3:end),'.wav'))
    % read the wav file name:
    [x,fs] = wavread(wavFileName);
else
    fprintf('Unknown file type!\n');
    return;
end


% Convert mono to stereo
if (size(x, 2)==2)
	x = mean(x')';
end

% Window length and step (in seconds):
win = 0.050;
step = 0.050;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  THRESHOLD ESTIMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%

Weight = 5; % used in the threshold estimation method

% Compute short-time energy and spectral centroid of the signal:
Eor = ShortTimeEnergy(x, win*fs, step*fs);
Cor = SpectralCentroid(x, win*fs, step*fs, fs);

% Apply median filtering in the feature sequences (twice), using 5 windows:
% (i.e., 250 mseconds)
E = medfilt1(Eor, 5); E = medfilt1(E, 5);
C = medfilt1(Cor, 5); C = medfilt1(C, 5);

% Get the average values of the smoothed feature sequences:
E_mean = mean(E);
Z_mean = mean(C);

% Find energy threshold:
[HistE, X_E] = hist(E, round(length(E) / 10));  % histogram computation
[MaximaE, countMaximaE] = findMaxima(HistE, 3); % find the local maxima of the histogram
if (size(MaximaE,2)>=2) % if at least two local maxima have been found in the histogram:
    T_E = (Weight*X_E(MaximaE(1,1))+X_E(MaximaE(1,2))) / (Weight+1); % ... then compute the threshold as the weighted average between the two first histogram's local maxima.
else
    T_E = E_mean / 2;
end

% Find spectral centroid threshold:
[HistC, X_C] = hist(C, round(length(C) / 10));
[MaximaC, countMaximaC] = findMaxima(HistC, 3);
if (size(MaximaC,2)>=2)
    T_C = (Weight*X_C(MaximaC(1,1))+X_C(MaximaC(1,2))) / (Weight+1);
else
    T_C = Z_mean / 2;
end

% Thresholding:
Flags1 = (E>=T_E);
Flags2 = (C>=T_C);
flags = Flags1 & Flags2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  SPEECH SEGMENTS DETECTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count = 1;
WIN = 5;
Limits = [];
while (count < length(flags)) % while there are windows to be processed:
	% initilize:
	curX = [];	
	countTemp = 1;
	% while flags=1:
	while ((flags(count)==1) && (count < length(flags)))
		if (countTemp==1) % if this is the first of the current speech segment:
			Limit1 = round((count-WIN)*step*fs)+1; % set start limit:
			if (Limit1<1)	Limit1 = 1; end        
		end	
		count = count + 1; 		% increase overall counter
		countTemp = countTemp + 1;	% increase counter of the CURRENT speech segment
	end

	if (countTemp>1) % if at least one segment has been found in the current loop:
		Limit2 = round((count+WIN)*step*fs);			% set end counter
		if (Limit2>length(x))
            Limit2 = length(x);
        end
        
        Limits(end+1, 1) = Limit1;
        Limits(end,   2) = Limit2;
    end
	count = count + 1; % increase overall counter
end

%%%%%%%%%%%%%%%%%%%%%%%
% POST - PROCESS      %
%%%%%%%%%%%%%%%%%%%%%%%

% A. MERGE OVERLAPPING SEGMENTS:
RUN = 1;
while (RUN==1)
    RUN = 0;
    for (i=1:size(Limits,1)-1) % for each segment
        if (Limits(i,2)>=Limits(i+1,1))
            RUN = 1;
            Limits(i,2) = Limits(i+1,2);
            Limits(i+1,:) = [];
            break;
        end
    end
end

% B. Get final segments:
segments = {};
for (i=1:size(Limits,1))
    segments{end+1} = x(Limits(i,1):Limits(i,2)); 
end

ThresholdVad = 0.02;

segments = cell2mat(segments);
cuttedSignal = [];

for i=1:size(segments)
    if segments(i) > ThresholdVad
        cutStart = i;
        break;
    end
end

for i=size(segments):-1:1
    if segments(i) > ThresholdVad
        cutEnd = i;
        break;
    end
end

for i=1:(cutEnd-cutStart)
    cuttedSignal(i) = segments(i+cutStart);
end

cuttedSignal = transpose(cuttedSignal);

Tw = 25;           % analysis frame duration (ms)
Ts = 10;           % analysis frame shift (ms)
alpha = 0.97;      % preemphasis coefficient
R = [ 300 3700 ];  % frequency range to consider
M = 20;            % number of filterbank channels 
C = 13;            % number of cepstral coefficients
L = 22;            % cepstral sine lifter parameter

% hamming window (see Eq. (5.2) on p.73 of [1])
hamming = @(N)(0.54-0.46*cos(2*pi*[0:N-1].'/(N-1)));

% Read speech samples, sampling rate and precision from file
[ speech, fs, nbits ] = wavread(wavFileName);

% Feature extraction (feature vectors as columns)
[ MFCCs, FBEs, frames ] = ...
                           mfcc( segments, fs, Tw, Ts, alpha, hamming, R, M, C, L );
                       
% Plot cepstrum over time
figure('Position', [30 100 800 200], 'PaperPositionMode', 'auto', ... 
                  'color', 'w', 'PaperOrientation', 'landscape', 'Visible', 'on' ); 
       
imagesc( [1:size(MFCCs,2)], [0:C-1], MFCCs ); 
axis( 'xy' );
xlabel( 'Frame index' ); 
ylabel( 'Cepstrum index' );
title( 'Mel frequency cepstrum' );                       