function reconstruct_interferogram(t_data:: Vector{Number}, t_sample:: Vector{Number}, spectogram_sample:: Vector{Number}, iDFT:: Array{Number}, interferogram_data:: Vector{Number})
    interferogram = iDFT*spectogram_sample

    plot(t_data, interferogram_data, color=:blue, label="Interferogram from iDFT")
    scatter!(t_sample, real(interferogram), color=:green, label="Interferogram from samples")
    plot!(title="Interferogram from iDFT of data vs. reconstructed interferogram from samples", xlabel="Optical path [cm]", ylabel="Volts [V]")
end

function plot_spectogram(hz_range:: Vector{Number}, spectogram_sample:: Vector{Number}, spectogram_data:: Vector{Number})
    plot(hz_range, real.(spectogram_data), color=:red, lw =2, alpha = .5 , label="Data")
    plot!(hz_range, real.(spectogram_sample), color=:green, lw =2, alpha = .5, linestyle=:dash, label="Spectogram from samples (real)")
    plot!(hz_range, imag.(spectogram_sample), color=:yellow3, linestyle=:dash, label="Spectogram from samples (imaginary)")
    plot!(title="Spectogram data vs. inverse problem", xlabel="Ratio of "*L"\lambda", ylabel="Fourier frequency term")
end

function plot_all(t_data:: Vector{Number}, t_sample:: Vector{Number}, spectogram_sample:: Vector{Number}, iDFT:: Array{Number}, interferogram_data:: Vector{Number}, hz_range:: Vector{Number}, spectogram_data:: Vector{Number}) 
    p₁, p₂=plot_spectogram(hz_range, spectogram_sample, spectogram_data), reconstruct_interferogram(t_data, t_sample, spectogram_sample, iDFT, interferogram_data)
end