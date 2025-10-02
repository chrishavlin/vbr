function [Prior_mod, sigmaPreds] = priorModelProbs( ...
    states, states_fields)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [Prior_mod, sigmaPreds] = priodModelProbs( ...
%       States, states_fields, ifnormal)
%
% Loop over the fields in States and calculate probabilities for each field
% to get total prior model pdf.  Each field (listed in states_fields) is
% a state variable name that we are varying.
%
% Assuming that the state variables are all independent of each other
%   p(var1, var2, ...) = p(var1) * p(var2) * ...
%
% Parameters
% ----------
% states: structure
%   A structure with the following fields for each state in the
%   list states_fields:
%    [field] : matrix
%       matrix of size (n_var1, n_var2, n_var3 ...) for all values of that
%       state variable
%    [field]_mean : scalar
%       mean (expected value) for that variable
%    [field]_std : scalar
%       standard deviation for that variable
%    [field]_pdf_type : string (optional)
%       If set, this is the PDF type to use (optional), must be one of
%       'normal', 'lognormal', 'uniformlog' or 'uniform'. This
%        option is overridden by the following field if it exists.
%    [field]_pdf : matrix (optional)
%       If there is a field [var_name]_pdf, then use that probability
%       as the prior for that variable instead of calculating it
%
% states_fields: cell array
%    names of all of the state variables we are varying. Each field name
%    should exist in the states structure.
%
% Returns
% -------
% [Prior_mod, sigmaPreds]
%   Prior_mod
%       joint probability of all combinations of the state variables
%   sigmaPreds
%       joint standard deviation for all combinations of the state variables
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  sigmaPreds = 1;
  marginals{numel(states_fields)} = 1;

  for i_field = 1:numel(states_fields)
    this_field=states_fields{i_field};
    std_field=[this_field,'_std']; % e.g., Tpot_std
    mn_field=[this_field,'_mean']; % e.g., Tpot_mean

    if isfield(states,[this_field,'_pdf_type'])
        pdf_type = states.([this_field,'_pdf_type']);
    else
        pdf_type = 'uniform';
    end

    % always override pdf_type if we have an input pdf
    if isfield(states,[this_field,'_pdf'])
        pdf_type = 'input';
    end

    switch pdf_type
        case 'input'
            marginals{i_field} = states.([this_field, '_pdf']);
            sigma = states.(std_field);
        case 'normal'
            % assume a normal distribution
            sigma =  states.(std_field); % standard deviation
            mu     = states.(mn_field); % mean value
            x      = states.(this_field); % measurements
            marginals{i_field} = probability_distributions(...
                'normal', x, mu, sigma);
        case 'lognormal'
            % assume a log normal distribution (assumes everything is properly normalized)
            sigma =  states.(std_field); % standard deviation
            mu     = states.(mn_field); % mean value
            x      = states.(this_field); % measurements
            marginals{i_field} = probability_distributions(...
                'lognormal', x, mu, sigma);
        case 'uniformlog'
            % uniform probability over natural log space
            x = log(states.(this_field));
            minv = min(x(:));
            maxv = max(x(:));
            sigma = 1;
            marginals{i_field} = probability_distributions(...
                'uniform', x, minv, maxv);
        otherwise
            % uniform PDF over total range
            sigma = 1;
            minv = min(states.(this_field)(:));
            maxv = max(states.(this_field)(:));
            x = states.(this_field); % measurements
            marginals{i_field} = probability_distributions(...
                'uniform', x, minv, maxv);
    end

    % Propagation of uncertainty for product of two real variables,
    %       f = A * B
    % sigma_f = |f| * sqrt((sigma_A/A)^2) + (sigma_B/B)^2 + 2(cov_AB/A/B))
    % from https://en.wikipedia.org/wiki/Propagation_of_uncertainty
    % cov_AB is the covariance of A and B, which is zero for our
    % (assumed) independent state variables
    sigmaPreds = sigmaPreds .* sigma;
  end

  Prior_mod = probability_distributions('joint independent', marginals);
end
