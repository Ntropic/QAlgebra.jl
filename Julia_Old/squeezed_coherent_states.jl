# In this file we define scripts to determine parameters of coherent or squeezed displaced StatsBase
# from operator expectation values. estimate errors of the approximation
# calculate overlaps (coming soon)

##### Coherent States ####################################################################################
##########################################################################################################

# take operators s, -s (for checks +-s) where s is a spin state  
function fit_coherent_state(s::ComplexF64, ms::ComplexF64)
    # - is given by -s/s 
    alpha = ms/s
    return alpha
end
function coherent_state_test(alpha::ComplexF64, val::ComplexF64, p::Int=1, q::Int=1; eps::Float64=1e-10)
    # value is of operator type +^p-^q (by default +-=n)
    expected_value = conj(alpha)^p*alpha^q
    relative_error = abs(expected_value-val)/abs(expected_value+eps)
    return relative_error
end

##### Displaced Squeezed States ##########################################################################
##########################################################################################################

## Displaced Squeezed states ( D(alpha)S(x)|0>=|x, alpha> )
# takes operators -, +- , -- to fit parameters 
function fit_displaced_squeezed_state(s::ComplexF64, ms::ComplexF64, mms::ComplexF64)
    # - is given by -s/s 
    m = ms/s
    pm = pms/s
    pm_angle = angle(pm)
    mm = mms/s
    alpha = m
    r = asinh(sqrt(real(pm - abs2(alpha))))
    theta = angle((mm-alpha^2)/(1/2*sinh(2*r)))
    return alpha, r, theta, pm_angle
end
function calculate_displaced_squeezed_state_expectation(alpha::ComplexF64, r::Float64, theta::Float64, p::Int=1, q::Int=2)
    # value is of operator type +^p-^q (by default +--)
    # a tranforms as a   ==> a cosh(r) + a' * exp(i*theta) sinh(r) + alpha
    # and a' does so via ==> a' cosh(r) + a * exp(-i*theta) sinh(r) + alpha^*
    if p+q > 4
        error("Currently only supports up to 4th order cavity operators")
    elseif p+q <= 2
        error("Cavity operators of order <= 2 are already used to determine the parameters of the discplaced squeezed state.")
    end
    cr = cosh(r)
    sr = sinh(r)
    etheta = exp(1im*theta)
    if p==0 && q==3 #---
        expected_value = alpha^3 + 3*alpha*etheta*cr*sr 
    elseif p==1 && q==2 #+--
        expected_value = conj(alpha)*alpha^2 + 2*alpha*sr^2 + conj(alpha) * etheta * cr*sr
    elseif p==0 && q==4 #----
        expected_value = alpha^4 + 6*alpha^2*cr * sr * etheta + 3*etheta^2*cr^2*sr^2
    elseif p==1 && q==3 #+---
        expected_value = conj(alpha)*alpha^3 + 3*alpha^2*sr^2 + 3*conj(alpha)*alpha * etheta * cr * sr + 3*etheta*sr^3*cr 
    elseif p==2 && q==2 #++--
        expected_value = conj(alpha)^2*alpha^2 + alpha^2*cr*sr*conj(etheta) + 4*conj(alpha)*alpha*sr^2 + conj(alpha)^2*etheta*cr*sr + cr^2*sr^2 + 2*sr^4
    end
    return expected_value
end
function displaced_squeezed_state_expectation(p::Int, q::Int)
    # value is of operator type +^p-^q (by default +--)
    # a transforms as a   ==> a cosh(r) + a' * exp(i*theta) sinh(r) + alpha
    # and a' does so via ==> a' cosh(r) + a * exp(-i*theta) sinh(r) + alpha^*
    if p+q > 4
        error("Currently only supports up to 4th order cavity operators")
    elseif p+q <= 2
        error("Cavity operators of order <= 2 are already used to determine the parameters of the displaced squeezed state.")
    end

    if p == 0 && q == 3 # ---
        return function (alpha::ComplexF64, r::Float64, theta::Float64)
            cr = cosh(r)
            sr = sinh(r)
            etheta = exp(1im * theta)
            return alpha^3 + 3 * alpha * etheta * cr * sr
        end
    elseif p == 1 && q == 2 # +--
        return function (alpha::ComplexF64, r::Float64, theta::Float64)
            cr = cosh(r)
            sr = sinh(r)
            etheta = exp(1im * theta)
            return conj(alpha) * alpha^2 + 2 * alpha * sr^2 + conj(alpha) * etheta * cr * sr
        end
    elseif p == 0 && q == 4 # ----
        return function (alpha::ComplexF64, r::Float64, theta::Float64)
            cr = cosh(r)
            sr = sinh(r)
            etheta = exp(1im * theta)
            return alpha^4 + 6 * alpha^2 * cr * sr * etheta + 3 * etheta^2 * cr^2 * sr^2
        end
    elseif p == 1 && q == 3 # +---
        return function (alpha::ComplexF64, r::Float64, theta::Float64)
            cr = cosh(r)
            sr = sinh(r)
            etheta = exp(1im * theta)
            conj(alpha) * alpha^3 + 3 * alpha^2 * sr^2 + 3 * conj(alpha) * alpha * etheta * cr * sr + 3 * etheta * sr^3 * cr
        end
    elseif p == 2 && q == 2 # ++--
        return function (alpha::ComplexF64, r::Float64, theta::Float64)
            cr = cosh(r)
            sr = sinh(r)
            etheta = exp(1im * theta)
            return conj(alpha)^2 * alpha^2 + alpha^2 * cr * sr * conj(etheta) + 4 * conj(alpha) * alpha * sr^2 + conj(alpha)^2 * etheta * cr * sr + cr^2 * sr^2 + 2 * sr^4
        end
    else
        error("Unsupported operator configuration.")
    end
end
function displaced_squeezed_state_test(alpha::ComplexF64, r::Float64, theta::Float64, val::ComplexF64, ; eps::Float64=1e-10)
    expected_value_fun = displaced_squeezed_state_expectation(p, q)
    expected_value = expected_value_fun(alpha, r, theta)
    relative_error = abs(expected_value-val)/abs(expected_value+eps)
    return relative_error
end
## Test 
#x = 1.0 + 1.0im
#r = abs(x)
#theta = angle(x)
#alpha = 2.0 + 2.0im
#s = 1.2 + 1.2im
#ms = alpha*s
#pms = (abs2(alpha) + sinh(r)^2)*s
#mms = (alpha^2+exp(im*theta)*sinh(2*r)/2)*s 
#alpha_fit, r_fit, theta_fit = fit_displaced_squeezed_state(s, ms, mms)

