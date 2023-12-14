using DataFrames
using CSV
using Plots


# Load data
nist_data_files = readdir("data/nist_data/")
filter!(data -> occursin(".csv", data), nist_data_files)
# for nist_csv_file in nist_data_files
# end
function plot_nist_data(nist_data_files)
    num_data = length(nist_data_files)
    plt = plot(layout=num_data, size=(800, 800))
    for (index, nist_csv_file) in enumerate(nist_data_files)
        df_nist_data = CSV.read("data/nist_data/" * nist_csv_file, DataFrame)
        nist_data_x = df_nist_data[:, :x]
        nist_data_y = df_nist_data[:, :y]
        plot!(plt, nist_data_x, nist_data_y, label=nist_csv_file[1:end-4], legend=:topright, xlabel="Wavenumbers (cm-1)", ylabel="(micromol/mol)-1m-1 (base 10)", subplot=index)
    end
    return plt
end

plot_nist_data(nist_data_files)
