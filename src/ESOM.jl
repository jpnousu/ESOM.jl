module ESOM

using StaticArrays
using Interpolations

include("types.jl")


function create_interpolation_functions()
    return IntFunctions(
        LinearInterpolation([-40.0, -5.0, -1.0, 25.0, 35.0, 60.0],
                            [0.0, 0.0, 0.2, 1.53, 1.53, 0.0]),
        LinearInterpolation([-40.0, -3.0, 0.0, 7.0, 60.0],
                            [0.0, 0.0, 1.3, 1.3, 0.0]),
        LinearInterpolation([-40.0, -5.0, 1.0, 20.0, 40.0, 80.0],
                            [0.0, 0.0, 0.2, 1.0, 1.0, 0.0]),
        LinearInterpolation([-40.0, -5.0, 1.0, 13.0, 25.0, 50.0],
                            [0.0, 0.0, 0.2, 1.0, 1.0, 0.0]),
        LinearInterpolation([-40.0, -5.0, 1.0, 27.5, 35.0, 60.0],
                            [0.0, 0.0, 0.2, 1.95, 1.95, 0.0]),
        LinearInterpolation([-40.0, -30.0, -20.0, -10.0, 0.0, 10.0, 20.0, 30.0, 40.0, 50.0],
                            [0.03125, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 8.0]),
        LinearInterpolation([0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.417, 
                             1.333, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 4.0],
                            [0.0, 0.004, 0.026, 0.074, 0.154, 0.271, 0.432, 0.64, 
                             0.899, 1.0, 1.0, 0.844, 0.508, 0.305, 0.184, 0.111, 
                             0.067, 0.04, 0.024, 0.0]),
        LinearInterpolation([0.0, 0.133, 1.333, 2.333, 4.0],
                            [0.0, 1.0, 1.0, 0.0, 0.0]),
        LinearInterpolation([0.0, 0.067, 0.5, 2.333, 4.0, 10.0],
                            [0.0, 0.0, 1.0, 1.0, 0.0, 0.0])
    )
end

function create_parameters()
    h_humus = 0.06
    rho_humus = 100.0
    frac_L = 0.1
    frac_F = 0.2
    frac_H = 0.7
    frac_leaf = 0.5
    frac_woody = 0.5
    nitrogen = 0.7
    lignin = 25.0
    adjust = 1.0
    litterN = 1.2
    pH = 3.4
    ash = 5.0
    enable_peatbottom = 1
    enable_peatmiddle = 1
    enable_peattop = 1
    C_modif = 1.05
    contpara = 100
    mu_k1 = 0.092 * (lignin / nitrogen)^-0.7396 * adjust
    mu_k2 = 0.0027 * (lignin / nitrogen)^-0.3917 * adjust
    mu_k3 = 0.062 * (lignin / nitrogen)^-0.3972 * adjust
    nu = clamp(0.701 * pH - 1.6018 - 0.038 * pH^2, 0, 1)

    return Parameters(h_humus, rho_humus, frac_L, frac_F, frac_H, frac_leaf, frac_woody,
                      nitrogen, lignin, adjust, litterN, pH, ash, enable_peatbottom,
                      enable_peatmiddle, enable_peattop, C_modif, contpara, mu_k1, mu_k2,
                      mu_k3, nu)
end

function initialize_storages(params::Parameters)
    """
    Compute initial storages based on parameters.

    Returns a tuple of initial masses for various compartments.
    """
    LL_mass = params.h_humus * params.frac_L * params.rho_humus * params.frac_leaf * 100 / 100.
    LW_mass = params.h_humus * params.frac_L * params.rho_humus * params.frac_woody * 100 / 100.
    FL_mass = params.h_humus * params.frac_F * params.rho_humus * params.frac_leaf * 100 / 100.
    FW_mass = params.h_humus * params.frac_F * params.rho_humus * params.frac_woody * 100 / 100.
    H_mass = params.h_humus * params.frac_H * params.rho_humus * 100 / 100.
    P1_mass = 93.5 # kg m-2
    P2_mass = 93.5
    P3_mass = 93.5

    return LL_mass, LW_mass, FL_mass, FW_mass, H_mass, P1_mass, P2_mass, P3_mass
end

function create_initial_B_matrix(initial_masses::NTuple{8, Float64})
    """
    Create the initial decomposition storage vector (B).

    Parameters:
        initial_masses: Tuple of initial masses returned by `initialize_storages`.

    Returns an 11-element static vector.
    """
    LL_mass, LW_mass, FL_mass, FW_mass, H_mass, P1_mass, P2_mass, P3_mass = initial_masses

    B = @SVector [
        0.0,        # Fresh litter as input: needles and fine roots, kg m-2
        0.0,        # Woody debris as input: branches and coarse roots, kg m-2
        LL_mass,    # Leaf litter storage, kg m-2
        LW_mass,    # Woody litter storage, kg m-2
        FL_mass,    # F material storage (leaf & fine roots), kg m-2
        FW_mass,    # Woody F material (branches and coarse roots), kg m-2
        H_mass,     # Humus material storage, kg m-2
        P1_mass,    # Peat 1 storage, kg m-2
        P2_mass,    # Peat 2 storage, kg m-2
        P3_mass,    # Peat 3 storage, kg m-2
        0.0         # Cumulative output in organic material, kg m-2
    ]
    return B
end

# Modify the calculate_rates function to return a Rates object
function calculate_rates(forc::Forcing, params::Parameters, f::IntFunctions)
    k1 = (0.002 + 0.00009 * params.ash + 0.003 * params.litterN) * min(0.1754 * exp(0.0871 * forc.tair), 1.0) * f.Phi1236(forc.wn) * params.nu
    k2 = clamp((0.00114 - 0.00028 * params.litterN) * f.T2(forc.tair) * f.Phi1236(forc.wn) * params.nu, 0.0, 1.0)
    k3 = clamp((0.04 - 0.003 * params.litterN) * f.T3(forc.tair) * f.Phi1236(forc.wn), 0.0, 1.0)
    k4 = 0.005 * params.litterN * f.T4(forc.tair) * f.Phi4(forc.wn)
    k5 = 0.007 * f.T5(forc.tair) * f.Phi5(forc.wn)
    k6 = 0.0006 * f.T6(forc.tp_top) * f.Phi1236(forc.wn)
    k7 = 0.00045 * f.T7(forc.tp_top) * forc.peat_w1 * params.enable_peattop
    k8 = 0.0001 * f.T7(forc.tp_middle) * forc.peat_w2 * params.enable_peatmiddle
    k9 = 0.0001 * f.T7(forc.tp_bottom) * forc.peat_w3 * params.enable_peatbottom * 0.1
    
    # Return a Rates struct with the calculated values
    return Rates(k1, k2, k3, k4, k5, k6, k7, k8, k9)
end;

function construct_A_matrix(params, rates)
    # Initialize static 11x11 matrix with zeros
    A = @SMatrix zeros(Float64, 11, 11)

    # Construct the diagonal matrix using rates and params
    k_diag = zeros(Float64, 11)
    k_diag[3] = 1 - (params.C_modif * rates.k1 + rates.k3)                  
    k_diag[4] = 1 - (params.C_modif * rates.k1 * params.mu_k1 + rates.k3 * params.mu_k3)                 
    k_diag[5] = 1 - (params.C_modif * rates.k2 + rates.k4 + rates.k5)                          
    k_diag[6] = 1 - (params.C_modif * rates.k2 * params.mu_k2 + rates.k4 + rates.k5)                           
    k_diag[7] = 1 - (params.C_modif * rates.k6)                          
    k_diag[8] = 1 - (params.C_modif * rates.k7)                        
    k_diag[9] = 1 - (params.C_modif * rates.k8)                         
    k_diag[10] = 1 - (params.C_modif * rates.k9)                        
    k_diag[11] = 1  

    # Subdiagonal vectors
    k_low0 = zeros(Float64, 11)
    k_low0[7] = rates.k4 + rates.k5                               
    k_low0[11] = params.C_modif * rates.k9

    k_low1 = zeros(Float64, 11)
    k_low1[3] = 1                                      
    k_low1[4] = 1                                      
    k_low1[5] = rates.k3                                    
    k_low1[6] = rates.k3 * params.mu_k3                        
    k_low1[7] = rates.k4 + rates.k5                               
    k_low1[11] = params.C_modif * rates.k8 

    # Higher subdiagonals
    k_low2 = zeros(Float64, 11); k_low2[11] = params.C_modif * rates.k7
    k_low3 = zeros(Float64, 11); k_low3[11] = params.C_modif * rates.k6
    k_low4 = zeros(Float64, 11); k_low4[11] = params.C_modif * rates.k2 * params.mu_k2
    k_low5 = zeros(Float64, 11); k_low5[11] = params.C_modif * rates.k2
    k_low6 = zeros(Float64, 11); k_low6[11] = params.C_modif * rates.k1 * params.mu_k1
    k_low7 = zeros(Float64, 11); k_low7[11] = params.C_modif * rates.k1

    # Fill the matrix diagonals
    for i in 1:11
        A = setindex(A, k_diag[i], i, i)
    end

    # Fill the subdiagonals
    for i in 2:11
        A = setindex(A, k_low0[i], i, i-1)
    end
    for i in 3:11
        A = setindex(A, k_low1[i], i, i-2)
    end
    for i in 4:11
        A = setindex(A, k_low2[i], i, i-3)
    end
    for i in 5:11
        A = setindex(A, k_low3[i], i, i-4)
    end
    for i in 6:11
        A = setindex(A, k_low4[i], i, i-5)
    end
    for i in 7:11
        A = setindex(A, k_low5[i], i, i-6)
    end
    for i in 8:11
        A = setindex(A, k_low6[i], i, i-7)
    end
    for i in 9:11
        A = setindex(A, k_low7[i], i, i-8)
    end

    return A
end;

function compute_A_x_B(A, B)
    return A * B
end;

end # module ESOM
