function [ DFF, Time, group_delay, Baseline] = DC_removal( DATA, varargin)
% syntax: [DATA, BASELINE] = DC_removal(DATA,'Fsampling', Fsampling, 'Fpass', Fpass, 'Fcut', Fcut, 'Ap', RippleInthePassBand, 'Channels', Channels 'Debug',Debug);
% example Vm = DC_removal( single(Vm),'Fsampling', Sampling, 'Fpass', 0.5, 'Fcut', 1.0,'Ap', 0.1, 'Channels', 128, 'Debug', 0);
% written by Ilija Uzelac, email: uzelaci@gmail.com
% DATA is represented with 3D matrix 2D+time 
%
% With this function every pixel trace is high-pass filtered
%
% Input data (fluorescence = change + baseline (DF + F) ) is first low-pass filtered calculating the baseline (F) which is subsequently subtracted from the input raw trace and then divided by baseline.
% Baseline is obtained by low-pass Kaiser filter which passes frequencies from 0Hz to Fpacing/2 and blocks all frequencies from Fpacing and above with 40db transitioning slope between Fpacing/2 and Fpacing 
% Filter is specified with the allowable ripple in the pass pand 0.1 dB
% ripple in the pass-band (Baseline) RippleInthePassBand 
% 
% GPU computation only !!
% DATA HAS TO BE IN THE SINGLE/DOUBLE PRECISION FORMAT
% F_sampling, sampling frequency in Hz, has to be specified. 
% F_cut, the hihest frequency at which DC removal falls to 40db . Has to be less than the pacing frequency in Hz, has to be specified.Greatly
% affects filter delay
% F_pass the lowest frequency to be removed, has t0 be 1/2 the F_pacing to
% persere alternans
% Channels, number of channels to be processed in parallel on GPU at the
% same time
% if Ap is not specified, 0.1 is used
% if Channels is not specified, 1000 is used
% if Debug is not specified, 0 is used
% 
% Note: cutting frequency higher than 1/2 Fpacing for Baseline determination affects the signal if alterans are present, making APs equal. Alternans will not
% be visible! Basically with alternans present frequency is effectivelly Fpacing/2 so
% filter needs to be contructed as (Fpacing/2 ) /2 = Fpacing/4
% 
% since filtered DC part is inherently delayed through the filter, input (DF+F) signal and baseline signal
% (F) needs to be alligned in time for a proper subtraction. 
% This is how analog filter works in matlab, which produces output before
% the output is valid advancing it for one group delay, which has to be
% trimmed
%         -------------------------------------- F + DF (input)
%  ------------------------------------- F (output)
% but one group delay of the output signal at the begining, and input data
% at the end
%https://www.mathworks.com/help/signal/ug/compensate-for-delay-and-distortion-introduced-by-filters.html
% Computeed DF/F = ((F + DF) - F) / F is additionally trimmed at the begining by one filter delay
% to remove filter artifacts/ settlement time. 
% In total, input data is trimmed by 2 filter group delays, at the end and
% the begining
% Baseline is saved as uint16 to save on space
% filtfilt function would intriduce zero-lag across all frequencies
% including no delay at DC but it is not implemented for GPU computation


% checking input function arguments
    if nargin<2,  display('Not enough input arguments'), return, end
    Fsampling=[];
    Fpacing=[];
    Ap=[];
    ROI=[];
    Debug=[];
    Channels = [];
    for i=1:length(varargin)
        if strcmp(varargin{i},'Fsampling'),Fsampling = varargin{i+1}; end
        if strcmp(varargin{i},'Fpass'),Fpass = varargin{i+1}; end
        if strcmp(varargin{i},'Fcut'), Fcut  = varargin{i+1}; end
        if strcmp(varargin{i},'Ap'   ),  Ap      = varargin{i+1}; end
        if strcmp(varargin{i},'Channels'  ),  Channels     = varargin{i+1}; end
        if strcmp(varargin{i},'Debug'),  Debug   = varargin{i+1}; end
    end
    
    if isempty(Fsampling),display('Fsampling frequency needs to be specified'); return, end
    if isempty(Fpass),display('Fpacing frequency needs to be specified'); return, end
    if isempty(Fcut),display('Fpacing frequency needs to be specified'); return, end
    if isempty(Ap), Ap = .1; end 
    if isempty(Channels), Channels = 64; end 
    if isempty(Debug), Debug = 0; end 
    
    if ndims(DATA) == 3
        RR = size(DATA,1);    
        CC = size(DATA,2);
        DATA = reshape(DATA, [RR*CC, size(DATA,3)]);
    end
   tic 
% constructing LP filter ( pass = [0, Fpacing/2], stop = [Fpacing, < ], transitionwindow = [Fpacing/2, Fpacing])
    %Fp  - passband frequency
    %Fst - Stopband Frequency
    %Ap  - Passband Ripple (dB)
    %Ast - Stopband Attenuation (dB)
    d = fdesign.lowpass('Fp,Fst,Ap,Ast', Fpass, Fcut, Ap, 40, Fsampling);
    % for kaiserwindow filter Ast parameter does not change the filter
    % length
    % 20,30 or 40db the filter length is the same
    LP = design(d,'kaiserwin');
    %LP = design(d,'equiripple');
    if Debug==1
        figure
        fvtool(LP);
        xlim([0 5]);
        text(Fpass/2,0+2,'Fpass','fontsize',14)
        text(Fcut,-20+2,'Fcut','fontsize',14)
    end
    group_delay = round(length(LP.Numerator)/2);
    display( [ 'Group filter delay = ', num2str(group_delay), ' samples'] )
    
    % checking minimal input DATA length > 2* filter group delays
        if 2*group_delay > size(DATA,2)
            display(strcat('length of your input data is too short. The minimal length needed is : ', num2str(2*group_delay),'frames. Filter length is greatly affected with specified pacing frequency specified'))
            return
        end
 
  DFF = zeros(size(DATA,1),size(DATA,2)-2*group_delay,'single');
 % check on output parameters
 if nargout>3
     Baseline = zeros(size(DATA,1),size(DATA,2)-2*group_delay,'uint16');
     % Baseline will be time aligned with DF
 else
     Baseline = [];
 end
 

% Filtering using CPU
    LP = LP.Numerator;
   
    h = waitbar(0); drawnow; waitbar(0,h, 'DC removal: ');
    NN = ceil(size(DATA,1)/Channels);
    size(DATA,1)/Channels;
    
    CPU_time = 0;
    CPU_transfer_time = 0;
    for i=1:NN
    % CPU processing
        Start = 1+(i-1)*Channels;
        End   = i*Channels;
        if i==NN
           End = size(DATA,1);
        end
                                                                            LapTime = toc;
        DC_withAC = single( DATA( Start:End ,:) );
        DC = DC_withAC; % pre initialization to be single
                                                                            CPU_transfer_time = CPU_transfer_time + toc - LapTime;
             
                                                                            LapTime = toc;
        DC = filter( LP,1, DC_withAC,[],2); % 2 = operate along second dimension
                                                                            CPU_time = CPU_time+toc - LapTime;
     
      % FILTERING, subtraction with time allignment.
                                                                            LapTime = toc;
        DC_withAC(:, end-group_delay+1:end ) = [];
        DC(:, 1:group_delay ) = [];
        DC_withAC = ( DC_withAC - DC) ./ DC;
                                                                           
        DFF( Start:End , : ) = single( DC_withAC(:, group_delay+1:end)); % trimming one filter delay in addition at begining
        if ~isempty(Baseline)
           Baseline( Start:End , : ) = uint16(( DC(:, group_delay+1:end) ));
        end
                                                                            CPU_transfer_time = CPU_transfer_time+toc - LapTime;
        
      % update time bar
        TotalTime = toc;
        temp = round(  (NN - i)/i* TotalTime ); % remaining number of rows times mean time per row
        waitbar(i/NN,h, strcat( ['DC removal: ' ''],num2str(temp),' s'));         
    end
    
    
    DFF = reshape(DFF,  [RR, CC, size(DFF,2)]);
    if ~isempty(Baseline)
        Baseline = reshape( Baseline, [RR, CC, size(Baseline,2)]);
    end
    
    close(h)  
    display( [ 'total CPU time = ', num2str(CPU_time)]  )
    display( [ 'total CPU transfer time = ', num2str(CPU_transfer_time)]  )
    display( ['total time =', num2str(TotalTime) ] )
    
    Time.TotalTime = TotalTime;
    Time.CPU_transfer_time = CPU_transfer_time;
    Time.CPU_time = CPU_time;