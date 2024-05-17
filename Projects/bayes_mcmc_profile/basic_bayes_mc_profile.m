function result = basic_bayes_mc_profile()
    clear
    path_to_top_level_vbr='../../';
    addpath(path_to_top_level_vbr)
    vbr_init

    model_settings.MCMC_max_steps = 100000;
    model_settings.progress_plots = 1;
    model_settings.MCMC_active_z_memory = 10000;
    model_settings.output_every_n = 500;
    model_settings.nz = 20;
    model_settings.zLAB = 100;
    model_settings.zmax = 200;
    model_settings.zcrust = 10;
    model_settings.Tpot = 1425; % deg C
    model_settings.dTdz_ad = 0.3; % deg/km adiabat
    model_settings.phi0 = 0.01;
    model_settings.melt_layer_thickness = 30;
    model_settings.dg_lith = 0.001; % m
    model_settings.dg_asth = 0.01; % m
    model_settings.rho_crust = 2800;
    model_settings.rho_mantle = 3300;
    model_settings.sig_MPa = .1;
    model_settings.P_surf_Pa = 1e5;
    model_settings.freq = [0.01];%, 0.1, 10];
    model_settings.noise_magnitude = 0.001;
    model_settings.VBR.in.elastic.methods_list={'anharmonic';'anh_poro'};
    model_settings.VBR.in.viscous.methods_list={'HK2003'};
    model_settings.VBR.in.anelastic.methods_list={'xfit_mxw'};

    priors.T_min_K = 0+273;
    priors.T_max_K = 1500+273;
    priors.T_dT_z = 700;
    priors.T_type = 'uniform';
    priors.phi_min = 0;
    priors.phi_max = 0.1;
    priors.phi_type = 'uniform';
    priors.dg_min_um = 0.0001 * 1e6;
    priors.dg_max_um = 0.05 * 1e6;
    priors.dg_type = 'uniform'; % log_uniform
    model_settings.priors = priors;
    model_settings.Qobs_std_frac = 0.01;
    model_settings.Vobs_std_frac = 0.01;

    model_settings = calculate_invariant_values(model_settings);

    % get the synthetic target to fit
    [VBR_target, Vs, Q] = get_synthetic_model(model_settings);
    [Vs_noisy, Q_noisy] = add_noise(Vs, Q, model_settings);
    Vs_std = model_settings.Vobs_std_frac * Vs;
    Q_std = model_settings.Qobs_std_frac * Q;

    % initialize model
    model = initialize_model(model_settings);
    % predictions from initial models
    [Vs_pred, Q_pred] = get_model_predictions(model, model_settings);
    model.Vs_pred = Vs_pred;
    model.Q_pred = Q_pred;
    model_settings.active_iz = 1;
    model_stats = get_model_stats(Vs_noisy, Vs_std, Q_noisy, Q_std, model, model_settings);

    p_trend = zeros(model_settings.MCMC_max_steps,1);
    p_active = zeros(model_settings.MCMC_active_z_memory, 1);
    i_active = 1;
    n_active = 0;
    i_save = 1;
    if model_settings.progress_plots == 1
        f = figure();
        subplot(5,1,1)
        plot(model_settings.init_conds.z, Vs_noisy,'k')
        subplot(5,1,2)
        semilogy(model_settings.init_conds.z, Q_noisy,'k')
    end

    for i = 1:model_settings.MCMC_max_steps

        % perturb and predict
        model_1 = perturb_model(model, model_settings);
        [Vs_pred, Q_pred] = get_model_predictions(model_1, model_settings);
        model_1.Vs_pred = Vs_pred;
        model_1.Q_pred = Q_pred;
        model_1_stats = get_model_stats(Vs_noisy, Vs_std, Q_noisy, Q_std, model_1, model_settings);

        % choose old or new model
        p_total = model_stats.total_posterior;
        p_1_total = model_1_stats.total_posterior;
%        disp([p_1_total, p_total])
        if p_1_total > p_total
            model = model_1;
            model_stats = model_1_stats;
            p_save = p_1_total;
        else
%            frac = rand(1);
%            if frac < p_1_total / p_total
%                model = model_1;
%                model_stats = model_1_stats;
%                p_save = p_1_total;
%            else
%                p_save = p_total;
%            end
            p_save = p_total;
        end

        % store current (will overwrite eventually)
        p_active(i_active) = p_save;
        i_active = i_active + 1;
        if i_active > model_settings.MCMC_active_z_memory
            i_active = 1 ;
            n_active = n_active +1;
        end
        pdiff = mean(abs(p_save-p_active))/mean(p_save);
%        disp([i_active, p_save, pdiff])

        if pdiff < 1e-10 && n_active > 0
            disp("move it!")
            update_plots(model_settings,model)
            % move to next iz
            if model_settings.active_iz < model_settings.nz
                model_settings.active_iz = model_settings.active_iz + 1;
                p_active(:) = 100000;
            end
        end

        if rem(i, model_settings.output_every_n) == 0
            disp(['step ', num2str(i), ': ', num2str(p_save)])
            p_trend(i_save) = p_save;
            i_save = i_save + 1;
            if model_settings.progress_plots == 1
                update_plots(model_settings, model)
            end
        end


    end

    % store output conveniently
    result = struct();
    result.synth.VBR_target = VBR_target;
    result.synth.Vs = Vs;
    result.synth.Q = Q;
    result.synth.Vs_noisy = Vs_noisy;
    result.synth.Q_noisy = Q_noisy;
    result.model_settings = model_settings;
    result.z = model_settings.init_conds.z;
    result.model = model;
    result.p_trend = p_trend;
end

function update_plots(model_settings,model)
    active_iz = model_settings.active_iz;
    subplot(5,1,1)
    hold all
    plot(model_settings.init_conds.z, model.Vs_pred)
    plot(model_settings.init_conds.z(active_iz), model.Vs_pred(active_iz),'.')
    subplot(5,1,2)
    hold all
    semilogy(model_settings.init_conds.z, model.Q_pred)
    subplot(5,1,3)
    hold all
    plot(model_settings.init_conds.z, model.T_K)
    subplot(5,1,4)
    hold all
    plot(model_settings.init_conds.z, model.phi)
    subplot(5,1,5)
    hold all
    semilogy(model_settings.init_conds.z, model.dg_um/1e6)
    pause(0.01)
end

function model_settings = calculate_invariant_values(model_settings);

    % calculations and arrays needing initialization only once
    init_conds = struct();
    m_s = model_settings;
    z = linspace(0,m_s.zmax,m_s.nz);

    model_size = size(z);
    rho_o = ones(model_size) * m_s.rho_crust;
    rho_o(z>m_s.zcrust) = m_s.rho_mantle;

    chi = zeros(model_size);
    chi(z>m_s.zcrust) = 1.0;

    model_settings.init_conds.z = z;
    model_settings.init_conds.z_m = z*1e3;

    model_settings.init_conds.chi = chi;
    model_settings.init_conds.rho_o = rho_o;
    model_settings.init_conds.sig_MPa = m_s.sig_MPa * ones(model_size);

    p = model_settings.priors;

    T_mean_K_z = (p.T_max_K - p.T_min_K) * z / m_s.zmax + p.T_min_K;
    p.T_min_K_z = T_mean_K_z - p.T_dT_z / 2.0;
    p.T_min_K_z(p.T_min_K_z<273) = 273;
    p.T_max_K_z = T_mean_K_z + p.T_dT_z / 2.0;

    p.dg_min_um_log = log(p.dg_min_um);
    p.dg_max_um_log = log(p.dg_max_um);

    model_settings.priors = p;

end


function [VBR_out, Vs, Q] = get_synthetic_model(model_settings)
    % create a synthetic test to perturb and retrieve
    m_s = model_settings;
    z = m_s.init_conds.z;

    T_test = z/m_s.zLAB * m_s.Tpot;
    T_test(z>=m_s.zLAB) = m_s.Tpot;
    T_test = T_test + z * m_s.dTdz_ad;

    phi = zeros(size(T_test));
    phi(z>=m_s.zLAB) = m_s.phi0;
    phi(z>=m_s.zLAB+m_s.melt_layer_thickness) = 0.0;

    dg_m = ones(size(T_test)) * m_s.dg_lith;
    dg_m(z>=m_s.zLAB) = m_s.dg_asth;

    SV.z = z;
    SV.dg_um = dg_m * 1e6;
    SV.phi = phi;
    SV.T_K = T_test + 273;
    SV.f = m_s.freq;

    VBR_out = forward_calculation(SV, model_settings);

    meth = model_settings.VBR.in.anelastic.methods_list{1};
    Vs = VBR_out.out.anelastic.(meth).V;
    Q = VBR_out.out.anelastic.(meth).Q;
end

function VBR = forward_calculation(SV, model_settings)

    % copy over the present model state variables
    VBR = struct();
    VBR.in.SV = SV;

    % calculate density and pressure for this model
    T_K = SV.T_K;
    rho_o = model_settings.init_conds.rho_o;
    rho = Density_Thermal_Expansion(rho_o, T_K, 0.9);
    P0 = model_settings.P_surf_Pa;
    [rho, P_Pa] = Density_Adiabatic_Compression(rho_o, model_settings.init_conds.z_m,P0);
    VBR.in.SV.P_GPa = P_Pa/1e9; % pressure [GPa]
    VBR.in.SV.rho = rho; % density [kg m^-3]

    % copy over remaining VBR settings
    VBR.in.SV.sig_MPa = model_settings.init_conds.sig_MPa;
    VBR.in.elastic.methods_list = model_settings.VBR.in.elastic.methods_list;
    VBR.in.viscous.methods_list = model_settings.VBR.in.viscous.methods_list;
    VBR.in.anelastic.methods_list = model_settings.VBR.in.anelastic.methods_list;

    % finally do the calculation
    [VBR] = VBR_spine(VBR);

end

function [Vs_noisy, Q_noisy] = add_noise(Vs, Q, model_settings);
    Vs_noisy = 2 * (rand(size(Vs)) - 0.5) * model_settings.noise_magnitude;
    Vs_noisy = Vs .* (1 + Vs_noisy);

    Q_noisy = 2 * (rand(size(Q)) - 0.5) * model_settings.noise_magnitude;
    Q_noisy = Q .* (1 + Q_noisy);
end

function Tval = draw_T(model_settings, iz)

    frac = rand(1); % random between 0,1

    Tmin = model_settings.priors.T_min_K_z(iz);
    Tmax = model_settings.priors.T_max_K_z(iz);

    Tval =  frac * (Tmax - Tmin) + Tmin;
end

function dg = draw_dg(model_settings)
    frac = rand(1); % random between 0,1
    dg_max_exp = model_settings.priors.dg_max_um_log;
    dg_min_exp = model_settings.priors.dg_min_um_log;

    dg_exp = frac * (dg_max_exp - dg_min_exp) + dg_min_exp;
    dg = exp(dg_exp);
end

function phi = draw_phi(model_settings)
    frac = rand(1); % random between 0,1
    phi_min = model_settings.priors.phi_min;
    phi_max = model_settings.priors.phi_max;

    phi = frac * (phi_max - phi_min) + phi_min;
end


function model = initialize_model(model_settings);

    zsize = size(model_settings.init_conds.z);

    T = zeros(zsize);
    dg = zeros(zsize);
    phi = zeros(zsize);
    for iz = 1:numel(T)
        T(iz) = draw_T(model_settings, iz);
        dg(iz) = draw_dg(model_settings);
        phi(iz) = draw_phi(model_settings);
    end

    model.T_K = T;
    model.dg_um = dg;
    model.phi = phi;
end

function [Vs_pred, Q_pred] = get_model_predictions(model, model_settings)
    model.f = model_settings.freq;
    SV = model;
    VBR = forward_calculation(SV, model_settings);

    meth = model_settings.VBR.in.anelastic.methods_list{1};
    Vs_pred = VBR.out.anelastic.(meth).V;
    Q_pred = VBR.out.anelastic.(meth).Q;
end

function model_stats = get_model_stats(Vs_obs, Vs_std, Q_obs, Q_std, model, model_settings)
    zsize = size(model_settings.init_conds.z);

    % priors
    p_T = zeros(zsize);
    p_dg = zeros(zsize);
    p_phi = zeros(zsize);
    joint_prior = zeros(zsize);

    likelihood = zeros(zsize);
    posterior = zeros(zsize);

    total_posterior = 0;
    if model_settings.active_iz > 0
        upper_zi = model_settings.active_iz;
    else
        upper_zi = numel(p_T);
    end
    for iz = upper_zi:upper_zi

        Tmin = model_settings.priors.T_min_K_z(iz);
        Tmax = model_settings.priors.T_max_K_z(iz);
        T_x = model.T_K(iz);
        p_type = model_settings.priors.T_type;
        p_T(iz) = probability_distributions(p_type, T_x, Tmin, Tmax);

        phi_min = model_settings.priors.phi_min;
        phi_max = model_settings.priors.phi_max;
        phi_x = model.phi(iz);
        p_type = model_settings.priors.phi_type;
        p_phi(iz) = probability_distributions(p_type, phi_x, phi_min, phi_max);

        dmin = model_settings.priors.dg_min_um_log;
        dmax = model_settings.priors.dg_max_um_log;
        dg_x = log(model.dg_um(iz));
        p_type = model_settings.priors.dg_type;
        p_dg(iz) = probability_distributions(p_type, dg_x, dmin, dmax);

        joint_prior(iz) = p_dg(iz) * p_phi(iz) * p_T(iz);

        Q_pred = model.Q_pred(iz);
        Vs_pred = model.Vs_pred(iz);

        Vs_obs_i = Vs_obs(iz);
        Vs_std_i = Vs_std(iz);
        Q_obs_i = Q_obs(iz);
        Q_std_i = Q_std(iz);

        prob_type = 'likelihood from residuals';
%        disp([Q_obs_i, Q_std_i, Q_pred])
        Q_likeli = probability_distributions(prob_type, log10(Q_obs_i), log10(Q_std_i), log10(Q_pred));
%        disp('likeli')
%        disp([Q_obs_i, Q_std_i, Q_pred, Q_likeli])
        Vs_likeli = probability_distributions(prob_type, Vs_obs_i, Vs_std_i, Vs_pred);
%        disp([Vs_obs_i, Vs_std_i, Vs_pred, Vs_likeli])
        likelihood(iz) = Q_likeli * Vs_likeli;

        posterior(iz) = likelihood(iz) * joint_prior(iz);

        total_posterior = total_posterior + posterior(iz);
    end

    model_stats.p_T = p_T;
    model_stats.p_dg = p_dg;
    model_stats.p_phi = p_phi;
    model_stats.joint_prior = joint_prior;
    model_stats.likelihood = likelihood;
    model_stats.posterior = posterior;
    model_stats.total_posterior = total_posterior;

end


function model_1 = perturb_model(model, model_settings);

    % randomly select a variable and depth to perturb
    var = randi([1,3]);
    if model_settings.active_iz > 0
        z_i = model_settings.active_iz;
    else
        z_i = randi([1, numel(model.T_K)]);
    end
    model_1.T_K = model.T_K;
    model_1.phi = model.phi;
    model_1.dg_um = model.dg_um;

    if var == 1
        model_1.T_K(z_i) = draw_T(model_settings, z_i);
%        disp('change T')
%        disp([model_1.T_K(z_i), model.T_K(z_i)])
    elseif var == 2
        model_1.phi(z_i) = draw_phi(model_settings);
%        disp('change phi')
%        disp([model_1.phi(z_i), model.phi(z_i)])
    else
        model_1.dg_um(z_i) = draw_dg(model_settings);
%        disp('change dg')
%        disp([model_1.dg_um(z_i), model.dg_um(z_i)])
    end

end