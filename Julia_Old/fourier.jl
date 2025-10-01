using FFTW
include("pulses.jl")

# fourier transform of a function handle (fun) in an interval of num_stds 
# first determine the std of the fourier transform via fft (mean and std of the fft)
function fft_mean_and_std(val_dict::Dict, pulse_param::Pulse_Param_Struct; num_points::Int=2^10, type::Symbol=:positive)
    # determine fft and freuencies 
    # type can be :symmetric (also accept negative frequencies) or :positive
    pulse_fun = val_dict["\\beta"]
    kappa = val_dict["\\kappa"]
    T = pulse_param.T
    ts = collect(range(0, stop=T, length=num_points))
    fs = fftshift(fftfreq(num_points, num_points/(T)))
    # absolute value of fs 
    if type == :positive
        fs = abs.(fs)
    end
    fft_pulse = fftshift(fft(pulse_fun.(ts)))
    # determine mean and std of the fft
    mean_fft = sum(fs .* abs2.(fft_pulse)) / sum(abs2.(fft_pulse))
    std_fft = sqrt(sum((fs .- mean_fft).^2 .* abs2.(fft_pulse)) / sum(abs2.(fft_pulse)))
    return mean_fft, std_fft#, fs, fft_pulse
end
## Test 
#mean_fft, std_fft = fft_mean_and_std(val_dict, pulse_param; num_points=2^6)

# function to determine the fourier transform of a function handle fun in an interval of num_stds
function fourier_in_interval(val_dict::Dict, pulse_param::Pulse_Param_Struct; num_stds=2.5, num_points=2^10, reltol::Float64=1e-3, type::Symbol=:symmetric)
    # check if type  is symmetric or positive
    if !(type in [:positive, :symmetric])
        error("type should be :positive or :symmetric")
    end
    mean_fft, std_fft = fft_mean_and_std(val_dict, pulse_param; num_points=num_points, type=:negative)
    if type == :positive
        min_fs, max_fs = 0, mean_fft + num_stds * std_fft
    elseif type == :symmetric
        max_fs = mean_fft + num_stds * std_fft
        min_fs = mean_fft - num_stds * std_fft
    end
    fs::Vector{Float64} = collect(range(min_fs, stop=max_fs, length=num_points))
    # calculate the fourier transform of the pulse in the interval of num_stds
    pulse_fun = val_dict["\\beta"]
    kappa = val_dict["\\kappa"]
    T = pulse_param.T
    # fourier formula: F_fun = 1/sqrt(2*pi) * int_{0}^{T} fun(t) * exp(-i * omega * t) dt
    # we use quadgk to calculate the integral
    F_fun(omega) = quadgk(t -> pulse_fun(t) * exp(-1im * omega * 2*pi * t), 0, T, rtol=reltol)[1] / sqrt(2*pi)
    F_funs = F_fun.(fs)
    return fs, F_funs
end