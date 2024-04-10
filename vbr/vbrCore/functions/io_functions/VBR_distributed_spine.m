function VBR_distributed_spine(VBR)

    % these will be user-inputs
    VBR.in.lut_dims.phi = linspace(0, 1, 1000);
    VBR.in.lut_dims.T_K = linspace(300, 2000, 800);
    VBR.in.lut_dims.ordering = {"phi"; "T_K"};
    VBR.in.lut_dims.include_ghosts = 0;

    n_dim = numel(VBR.in.lut_dims.ordering);

    % check for each ordering val in lut_dims



end


function phi = phi_vals(ival)
    phi = (ival - 0)./(1. - .0)
end