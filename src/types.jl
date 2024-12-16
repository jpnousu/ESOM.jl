# type structures

struct IntFunctions
    T2
    T3
    T4
    T5
    T6
    T7
    Phi1236
    Phi4
    Phi5
end

struct Parameters
    h_humus::Float64 # Mor humus thickness in (m)
    rho_humus::Float64 # Bulk density of the mor (kg m3)
    frac_L::Float64 # Share of undecomposed litter (L) from the mor thickness (fraction 0...1)
    frac_F::Float64 # Share of partly decomposed F material from the mor thickness (fraction 0...1)
    frac_H::Float64 # Share of humified H material from the mor thickness (fraction 0...1)
    frac_leaf::Float64 # Share of non-woody material from L and F material (fraction 0...1)
    frac_woody::Float64 # Share of non-woody material from L and F material (fraction 0...1)
    nitrogen::Float64 # content in gravimetric %
    lignin::Float64 # content in gravimetric %
    adjust::Float64 # Adjustment factor
    litterN::Float64 # Litter nitrogen content
    pH::Float64 # pH value
    ash::Float64 # Ash content in %
    enable_peatbottom::Int64 # Enable peat bottom (1: enabled, 0: disabled)
    enable_peatmiddle::Int64 # Enable peat middle (1: enabled, 0: disabled)
    enable_peattop::Int64 # Enable peat top (1: enabled, 0: disabled)
    C_modif::Float64 # C modification factor
    contpara::Int64 # Continuous parameter
    mu_k1::Float64 # Rate constant k1
    mu_k2::Float64 # Rate constant k2
    mu_k3::Float64 # Rate constant k3
    nu::Float64 # nu parameter from pH
end

mutable struct Forcing
    L0w::Float64
    L0nw::Float64
    tair::Float64
    tp_top::Float64
    tp_middle::Float64
    tp_bottom::Float64
    wn::Float64
    peat_w1::Float64
    peat_w2::Float64
    peat_w3::Float64
end

# Define a mutable struct to hold the rates
mutable struct Rates
    k1::Float64
    k2::Float64
    k3::Float64
    k4::Float64
    k5::Float64
    k6::Float64
    k7::Float64
    k8::Float64
    k9::Float64
end