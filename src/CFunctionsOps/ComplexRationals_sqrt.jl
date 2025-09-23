issquare(n::Int) = n ≥ 0 && isqrt(n)^2 == n

# √ of a rational if it’s a rational square; else nothing
function sqrt_rational(q::Rational{Int})
    n, d = numerator(q), denominator(q)
    n ≥ 0 || return nothing
    rn, rd = isqrt(n), isqrt(d)
    (rn*rn == n && rd*rd == d) ? rn // rd : nothing
end

# canonicalize: c > 0 and gcd(a,b,c) = 1
function normalize(z::ComplexRational)
    a, b, c = z.a, z.b, z.c
    c == 0 && throw(ArgumentError("ComplexRational: denominator c must be nonzero"))
    if c < 0
        a = -a; b = -b; c = -c
    end
    g = gcd(gcd(abs(a), abs(b)), c)
    g == 0 && return ComplexRational(0, 0, 1)  # covers a=b=0
    ComplexRational(div(a,g), div(b,g), div(c,g))
end

# principal rational sqrt (if it exists): Re ≥ 0; if Re == 0 then Im ≥ 0
function principal_sqrt(z::ComplexRational)
    z = normalize(z)
    a, b, c = z.a, z.b, z.c

    # |z| = sqrt((a^2+b^2)/c^2) must be rational ⇒ a^2+b^2 is a perfect square
    n2 = a^2 + b^2
    issquare(n2) || return nothing
    norm = isqrt(n2) // c
    re   = a // c

    s1 = (norm + re) // 2   # = x^2
    s2 = (norm - re) // 2   # = y^2

    x = sqrt_rational(s1); x === nothing && return nothing
    y0 = sqrt_rational(s2); y0 === nothing && return nothing

    # choose sign for y to match principal branch:
    # normally sign(y) = sign(b); but for negative reals (b==0,a<0) take y ≥ 0
    y = if b < 0
        -y0
    elseif b > 0
        y0
    else
        (a < 0 ? y0 : 0//1)  # pure imaginary for negative real; 0 for nonnegative real
    end

    # ensure principal (Re ≥ 0; if Re == 0 then Im ≥ 0)
    if x == 0//1 && y < 0
        y = -y
    end

    # pack into (A + iB) / C with common denominator and normalize
    dx, dy = denominator(x), denominator(y)
    C = lcm(dx, dy)
    A = numerator(x) * (C ÷ dx)
    B = numerator(y) * (C ÷ dy)
    normalize(ComplexRational(A, B, C))
end

has_rational_sqrt(z::ComplexRational) = principal_sqrt(z) !== nothing