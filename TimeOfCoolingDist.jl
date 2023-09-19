## --- Read QTQT output

using StatGeochem
# Make sure we're running in the directory where the script is located
cd(@__DIR__)
# Read data from file
filenamebase = "QTQt-OUTPUT-FILE"
ext = ".txt"

# print step is either "2" if there are no zero/values between t-T pairs, if value is given, set to "3"
qtqt_print_step = 2
filepath = joinpath(@__DIR__, filenamebase*ext)

lines = readlines(filepath)
limits = findall(x->contains(x,"CHAIN"),lines)
start, stop = limits[1]+1, limits[2]-1

# Parse into array-of-arrays
parsed = lines[start:stop] .|> x -> delim_string_parse(x, ' ', Float64)

## --- Calculate time of first cooling through a given isotherm
# enter time window to search for cooling signal
starttime = 750
endtime = 500
dt = -1 # Must be negative!
Tq = 120 # half-max isotherm temperature

time_interp = starttime:dt:endtime
coolingdist = Array{Float64}(undef, length(parsed))
for i=1:length(parsed)
    time = parsed[i][5:qtqt_print_step:end]
    temperature = parsed[i][6:qtqt_print_step:end]
    temperature_interp = linterp1s(time, temperature, time_interp)

    coolingdist[i] = NaN
    if temperature_interp[1] > Tq
        @inbounds for j = 2:length(time_interp)
            if temperature_interp[j] <= Tq
                coolingdist[i] = time_interp[j]
                break
            end
        end
    end
end

# Save results
exportdataset((;coolingdist=coolingdist), "$(filenamebase)-coolingdist-$(starttime)-$(endtime)Ma-$(Tq)C.csv", ',')
## --- Plot results

using Plots
h = histogram(coolingdist, normalize=true, label="", xticks=endtime:50:starttime, xminorticks=5, tick_dir=:out, framestyle=:box)
plot!(h, xlabel="Time of first cooling through half-max temperature (Ma)", ylabel="Probability Density", )
plot!(h, xlims=(endtime, starttime), ylims=(0, ylims(h)[2]))
savefig(h, "$(filenamebase)-coolingdist-$(starttime)-$(endtime)Ma-$(Tq)C.pdf")
display(h)

# Calculate histogram
binedges = endtime:1:starttime
N = histcounts(coolingdist, binedges)
bincenters = cntr(binedges)
peak_index = argmax(N)
peak = N[peak_index]
peak_time = bincenters[peak_index]

# Calculate FDHM: full duration at half-maximum
left_fdhm = linterp1s(N[1:peak_index], bincenters[1:peak_index], peak/2)
right_fdhm = linterp1s(N[peak_index:end], bincenters[peak_index:end], peak/2)
fdhm = right_fdhm - left_fdhm
flush(stdout)
println("Cooling peak: $(peak_time) + $(right_fdhm - peak_time) - $(peak_time - left_fdhm) Ma")

## --- End of File
