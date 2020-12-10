% one column with one excitatory population and one PV+ population, with multiple
% neurons: prediction regarding variance explained by divisive model and subtractive model
% as a function of magnitude of heterogeneity of PV+ activation
%

clear;

% time parameters
t0=0;
tmax=0.4;
tstep=0.001;
t=t0:tstep:tmax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% network parameters
tau_ex=0.01;         % synaptic time constant (10ms)
tau_pv=0.01;         % synaptic time constant (10ms)

supra_n = 2.;        % supra non-linearity

% functions
heavisd = @(x) x>0;                 % heaviside function definition, including h(0)=0
thresh_lin= @(x) heavisd(x).*x;     % threshold linear with no saturation function definition, including h(0)=0
gaus_fun = @(x,mu,sigm) exp(-((x-mu).^2)/(2*sigm^2))/sqrt(2*pi*sigm^2);
iofunc_E = @(x) thresh_lin(x).^supra_n;
iofunc_I = @(x) thresh_lin(x).^supra_n;

% number of excitatory and inhibitory neurons
n_exc = 400;
n_inh = 100;

seed = 73;
rng(seed);

% external input to excitatory and inhibitory cells
h_exc = 2*rand(n_exc,1);
h_inh = 2*rand(n_inh,1);

% connectivity
sparsity = 0;
wexex_1 = rand(n_exc);
wpvex_1 = rand(n_inh,n_exc);
wexpv_1 = rand(n_exc,n_inh);
wpvpv_1 = rand(n_inh);
wexex_2 = .2/n_exc;
wpvex_2 = .6/n_exc;
wexpv_2 = -.3/n_inh;
wpvpv_2 = -1.7/n_inh;
wexex_3 = rand(n_exc);
wpvex_3 = rand(n_inh,n_exc);
wexpv_3 = rand(n_exc,n_inh);
wpvpv_3 = rand(n_inh);

wexex = wexex_2*wexex_1.*(wexex_3>sparsity);
wpvex = wpvex_2*wpvex_1.*(wpvex_3>sparsity);
wexpv = wexpv_2*wexpv_1.*(wexpv_3>sparsity);
wpvpv = wpvpv_2*wpvpv_1.*(wpvpv_3>sparsity);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E1 = 10;    % auditory stimulus maximum magnitude

% number of stimulus frequencies
m=15;

% characteristic frequency
charact_freq = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SSFO stimulation of PV+: degree of heterogeneity (mn-a)/sqrt(3)
ssfo_min = [0 .1:.2:1.9];
ssfo_mn = [0 1.9*ones(1,length(ssfo_min)-1)];
ssfo_mat = ones(n_inh,1)*ssfo_min + 2*rand(n_inh,1)*(ssfo_mn-ssfo_min);
ssfo_std = std(ssfo_mat);

num_exp = length(ssfo_min)-1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% magnitude of neural or synaptic noise
nois_fact = 0.03;

num_rep = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate network
X = zeros(n_exc,m,num_rep);
Y = zeros(n_inh,m,num_rep);
X_PV_SSFO = zeros(n_exc,m,num_exp,num_rep);
Y_PV_SSFO = zeros(n_inh,m,num_exp,num_rep);

% column i has the highest thalamic input
for rep=1:num_rep
    nois_exc_mat = nois_fact*randn(n_exc,m,length(t),2)/sqrt(tstep);
    nois_inh_mat = nois_fact*randn(n_inh,m,length(t),2)/sqrt(tstep);
    
    for kk=1:1:num_exp+1
        for i=1:m
            x=zeros(n_exc,length(t));
            y=zeros(n_inh,length(t));

            burst=zeros(1,length(t));
            burst(round(0.1/tstep):round(0.2/tstep)) = gaus_fun(i,charact_freq,3)*exp(-(0:tstep:0.1)/.05);

            % noise for control or SSFO-on condition
            if kk==1
                kkk = 1;
            else
                kkk = 2;
            end

            for j=2:length(t)
                nois_exc = nois_exc_mat(:,i,j,kkk);
                nois_inh = nois_inh_mat(:,i,j,kkk);
                x_input = wexex*x(:,j-1) + wexpv*y(:,j-1) + E1*burst(j-1) + h_exc + nois_exc;
                y_input = wpvex*x(:,j-1) + wpvpv*y(:,j-1) + E1*burst(j-1) + ssfo_mat(:,kk) + h_inh + nois_inh;
                x(:,j)=x(:,j-1)+tstep*(-x(:,j-1)+iofunc_E( x_input ))/tau_ex;
                y(:,j)=y(:,j-1)+tstep*(-y(:,j-1)+iofunc_I( y_input ))/tau_pv;
            end

            if kk==1
                X(:,i,rep) = mean(x(:,round(0.1/tstep):round(0.19/tstep)),2);
                Y(:,i,rep) = mean(y(:,round(0.1/tstep):round(0.19/tstep)),2);

            else
                X_PV_SSFO(:,i,kk-1,rep) = mean(x(:,round(0.1/tstep):round(0.19/tstep)),2);
                Y_PV_SSFO(:,i,kk-1,rep) = mean(y(:,round(0.1/tstep):round(0.19/tstep)),2);

            end

        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute slope and intercept of SSFO-on tuning vs control tuning, and
% error of purely divisive model (intercept constrained at zero), or
% error of purely subtractive model (slope constrained at one)

% excitatory neurons
error_area_div = zeros(n_exc,num_exp,num_rep);
error_area_sub = zeros(n_exc,num_exp,num_rep);
subselect_suppressed = zeros(n_exc,num_exp,num_rep);
for rep=1:num_rep
    for i = 1:n_exc
        for j = 1:num_exp
            tuning_off = X(i,:,rep);
            tuning_on = X_PV_SSFO(i,:,j,rep);
            max_tuning = max(tuning_off);
            
            % suppressed or enhanced
            subselect_suppressed(i,j,rep) = mean(X_PV_SSFO(i,:,j,rep)-X(i,:,rep))<0;

            % normalize
            tuning_off = tuning_off/max_tuning;
            tuning_on = tuning_on/max_tuning;

            % exclude low activities. Note that if PV+ activation is high,
            % activities can become very low and thus many samples might be
            % eliminated, which will bias these regressions
            thresh_act = 0.05;
            invalid = tuning_on < thresh_act | tuning_off < thresh_act;
            tuning_off(invalid) = [];
            tuning_on(invalid) = [];
            
            % regression
            [slope_on1,intercept_on1,corr1,sm1,sb1] = lsqfitma(tuning_off',tuning_on');
            error_area_div(i,j,rep) = abs(tuning_off-(1/slope_on1)*tuning_on)*abs(tuning_on'-slope_on1*tuning_off');
            error_area_sub(i,j,rep) = abs(tuning_off-tuning_on+intercept_on1)*abs(tuning_on'-tuning_off'-intercept_on1);
            if isnan(slope_on1)
                error_area_div(i,j,rep) = 0;
                subselect_suppressed(i,j,rep) = 0;
            end
            
            if isnan(intercept_on1)
                error_area_sub(i,j,rep) = 0;
                subselect_suppressed(i,j,rep) = 0;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mean and std of errors for divisive model, for suppressed cells (and not NaNs)
mn_error_area_div1 = sum(error_area_div.*subselect_suppressed,1)./sum(subselect_suppressed,1);
mn_error_area_div = mean(mn_error_area_div1,3);
std_error_area_div = std(mn_error_area_div1,0,3);

% mean and std of errors for subtractive model, for suppressed cells (and not NaNs)
mn_error_area_sub1 = sum(error_area_sub.*subselect_suppressed,1)./sum(subselect_suppressed,1);
mn_error_area_sub = mean(mn_error_area_sub1,3);
std_error_area_sub = std(mn_error_area_sub1,0,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting
clf(figure(1))
lw_plt = 1;
fontsize_plt = 13;
plot(ssfo_std(2:end),mn_error_area_div,'-k','linewidth',lw_plt); hold on;
plot(ssfo_std(2:end),mn_error_area_sub,'--k','linewidth',lw_plt); hold on;
h11=patch(horzcat(ssfo_std(2:end),fliplr(ssfo_std(2:end))),...
    horzcat(mn_error_area_div-std_error_area_div, fliplr(mn_error_area_div+std_error_area_div)), -ones(1,2*length(ssfo_std(2:end))),'k'); hold on;
h12=patch(horzcat(ssfo_std(2:end),fliplr(ssfo_std(2:end))),...
    horzcat(mn_error_area_sub-std_error_area_sub, fliplr(mn_error_area_sub+std_error_area_sub)), -ones(1,2*length(ssfo_std(2:end))),'k'); hold on;
set(h11,'FaceAlpha',0.3,'edgecolor','none');
set(h12,'FaceAlpha',0.3,'edgecolor','none');
set(gca,'TickDir','out','linewidth',lw_plt,'fontsize',fontsize_plt);
xlabel('PV-activation std')
ylabel('model error')
box off
axis tight
legend('Div exc','Sub exc')
legend('boxoff')