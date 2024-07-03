lambda = 10.^(-0.5:-0.5:-12.5);
order=0; 
reg_param_method = 'g';
compute_params=1;
model='SR'
fs=500;
SNR=1;
Nsamples=2500
f_low=0.5;
f_high=30;
switch model
    case 'SR'
        samples = 2:Nsamples;   %Sinusal sampling. Instant 1 is erased.
    case 'SAF'
        samples = length(EGMs(1,:))-Nsamples-1500+2:length(EGMs(1,:))-1500;
    case 'CAF'
        samples = length(EGMs(1,:))-Nsamples-2000+2:length(EGMs(1,:))-2000;
        end
x = EGMs(:,samples);
[x, y, Cn, A] = preprocess (x, MTransfer, model, SNR, fs, f_low, f_high);
x=x(:,500:999);
y=y(:,500:999);
[AA, L, LL] = precompute_matrices (A, order, atrial_model);
[x_hat, param_opt] = tikhonov (A, L, AA, LL, y, lambda, SNR, order, reg_param_method, compute_params);
reg_corner=sqrt(param_opt);
for i=1:2048
    for j=1:500
Err(i,j)=(norm(x(i,j)-x_hat(i,j))^2)/500;

    end
end
Err_root=sqrt(Err);