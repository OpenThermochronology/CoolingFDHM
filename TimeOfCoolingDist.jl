## --- Read QTQT output

using StatGeochem
# Make sure we're running in the directory where the script is located
cd(@__DIR__)
# Read data from file
filenamebase = "13LVA04_v5.8.0_v3"
ext = ".txt"
qtqt_print_step = 2
filepath = joinpath(@__DIR__, filenamebase*ext)

lines = readlines(filepath)
limits = findall(x->contains(x,"CHAIN"),lines)
start, stop = limits[1]+1, limits[2]-1

# Parse into array-of-arrays
parsed = lines[start:stop] .|> x -> delim_string_parse(x, ' ', Float64)

## --- Plot image with colormap, to make sure it has been read properly

#using Plots
#using StatsBase: fit, Histogram

#xresolution = 2000
#yresolution = 1000
#tmax = nanmean(parsed .|> x -> maximum(x[5:qtqt_print_step:end-1]))
#Tmax = maximum(parsed .|> x -> maximum(x[6:qtqt_print_step:end]))

# Resize the post-burnin part of the stationary distribution
#tTdist = Array{Float64}(undef, xresolution, length(parsed))
#xq = range(0,tmax,length=xresolution)
#for i=1:length(parsed)
#    tTdist[:,i] = linterp1s(parsed[i][5:qtqt_print_step:end-1],parsed[i][6:qtqt_print_step:end],xq)
#end

# Calculate composite image
#tTimage = zeros(yresolution,size(tTdist,1))
#yq = range(0,Tmax,length=yresolution+1)
#for i=1:size(tTdist,1)
#    tTimage[:,i] = fit(Histogram,tTdist[i,.!isnan.(tTdist[i,:])],yq,closed=:right).weights
#end

# Plot image with colorscale # ylcn ramp
#A = imsc(tTimage,ylcn,0,nanpctile(tTimage,98.5))
#h = plot(xlabel="Time (Ma)",ylabel="Temperature (°C)",yticks=0:50:Tmax,xticks=0:20:tmax,yminorticks=2,xminorticks=2,tick_dir=:out,framestyle=:box)
#plot!(h,xq,cntr(yq),A,yflip=false,xflip=false,legend=false,aspectratio=tmax/Tmax/1.5,xlims=(0,tmax),ylims=(0,Tmax))

## --- Calculate time of first cooling through a given isotherm

starttime = 750
endtime = 500
dt = -1 # Must be negative!
Tq = 120

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
