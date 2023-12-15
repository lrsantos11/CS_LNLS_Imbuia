function reconstruct_interferogram(t_data:: Vector{<:Number}, t_sample:: Vector{<:Number}, spectrum_sample:: Vector{<:Number}, iDFT:: Array{<:Number}, interferogram_data:: Vector{<:Number})
    interferogram = iDFT*spectrum_sample

    plot(t_data, interferogram_data, color=:blue, label="Interferogram from iDFT")
    scatter!(t_sample, real(interferogram), color=:green, label="Interferogram from samples")
    plot!(title="Interferogram from iDFT of data vs. reconstructed interferogram from samples", xlabel="Optical path [cm]", ylabel="Volts [V]")
end

function plot_spectrum(frequency_range:: Vector{<:Number}, spectrum_sample:: Vector{<:Number}, spectrum_data:: Vector{<:Number})
    plot(frequency_range, real.(spectrum_data), color=:red, lw =2, alpha = .5 , label="Data")
    plot!(frequency_range, real.(spectrum_sample), color=:blue, lw =2, alpha = .5, label="Spectrum from samples (real)")
    plot!(frequency_range, imag.(spectrum_sample), color=:green, label="Spectrum from samples (imaginary)")
    plot!(title="spectrum data vs. inverse problem", xlabel="Ratio of "*L"\lambda", ylabel="Fourier frequency term")
end

function plot_all(frequency_range:: Vector{<:Number}, spectrum_sample:: Vector{<:Number}, spectrum_data:: Vector{<:Number}, t_sample:: Vector{<:Number}, t_data:: Vector{<:Number}, iDFT:: Array{<:Number}, interferogram_data:: Vector{<:Number}) 
    p₁, p₂=plot_spectrum(frequency_range, spectrum_sample, spectrum_data), reconstruct_interferogram(t_data, t_sample, spectrum_sample, iDFT, interferogram_data)
end