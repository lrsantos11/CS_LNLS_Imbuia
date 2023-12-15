using StatsBase
using FFTW
using LinearAlgebra
using Plots
include("../tests/tests_utils.jl")
include("../src/data_utils.jl")
include("../src/proximal_subgradient_method.jl")
include("../plots/plots_util.jl")

function get_s(methods:: Vector{Symbol}, iDFT:: Array{<:Number}, interferogram:: Vector{<:Number}; λ=1.0)
    # Your method here
    spectrum = 0

    for mtd ∈ methods
        func     = eval(mtd)
        spectrum = func(iDFT, interferogram; λ=λ)
    end

    return spectrum
end

function get_spectogram_from_sample(interferogram::Vector{Number}, percent:: Number)
    n = length(interferogram)
    sampleₙ, interferogramₙ, iDFTₙ = get_normally_dist_data(collect(1:n), interferogram, n)

    # m sub-samples
    rows_id        = sample(1:n, round(Int, percent*n))
    sampleₘ        = sampleₙ[rows_id]
    interferogramₘ = interferogramₙ[rows_id]
    iDFTₘ          = iDFTₙ[rows_id, :]

    return sampleₘ, interferogramₘ, iDFTₘ
end

nist_csv_file_re = "63148-62-9-IR_Re_Silicone oil.csv"
df_nist_data_re  = CSV.read("data/nist_data/" * nist_csv_file_re, DataFrame)
nist_csv_file_im = "63148-62-9-IR_Im_Silicone oil.csv"
df_nist_data_im  = CSV.read("data/nist_data/" * nist_csv_file_im, DataFrame)
frequency_range  = df_nist_data_re[:, :x]
spectrum_real    = df_nist_data_re[:, :y]
spectrum_imag    = df_nist_data_im[:, :y]

interferogram = ifft(spectrum_imag)

# If you want to try the full signal... 
#spectrum_full = spectrum_real+im*spectrum_imag
#interferogram_full = ifft(spectrum_full)
#get_spectogram_from_sample(methods, interferogram_full, percent)

p₁=plot(frequency_range, real.(fft(interferogram)), label="Data")
plot!(frequency_range, spectrum_imag, label="Reconstructed from double FFT")
plot!(title="Interferogram of silicone")

percent=0.5

sampleₘ, interferogramₘ, iDFTₘ, spectrumₘ = get_spectogram_from_sample(interferogram, percent)

p₂, p₃ = plot_all(frequency_range, spectrumₘ, spectrum_imag, sampleₘ, [i for i=1:length(interferogram)], iDFTₘ, real.(interferogram))
plot(p₂, p₃)