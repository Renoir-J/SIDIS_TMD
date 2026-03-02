using QuadGK
using PolyLog
using SpecialFunctions
using StaticArrays
using .FastGK

include("../../Wrappers/NangaParbat-master/TMDPDF_MAP.jl")
tmdpdf_init(name=config_name)

# active collinear pdfs

# https://arxiv.org/pdf/1604.07869

const active_flavors = ("u","ub","d","db","s","sb","c","cb","b","bb")

# anti-quark rules

function anti_flavor(flavor)
    if flavor == "u"
        type = "ub"
    elseif flavor == "ub"
        type = "u"
    elseif flavor == "d"
        type = "db"
    elseif flavor == "db"
        type = "d"
    elseif flavor == "s"
        type = "sb"
    elseif flavor == "sb"
        type = "s"
    elseif flavor == "c"
        type = "cb"
    elseif flavor == "cb"
        type = "c"
    elseif flavor == "b"
        type = "bb"
    elseif flavor == "bb"
        type = "b"
    elseif flavor == "g"
        type = "g"
    else
        error("Unknown flavor: $flavor")
    end

    return type
end

# NLO

# qq

C1_PDF_qq_Reg(x) = (
    2*CF*(1-x)
)

C1_PDF_qq_Delta = (
    -CF*π^2/6
)

# qg

C1_PDF_qg_Reg(x) = (
    4*TF*x*(1-x) #Alexey: 4*TF*x*(1-x) #MAP: 8*TF*x*(1-x)
)

# Convolution

function F1_q_to_h(; x::Float64, μ::Float64, flavor::String)

    fq(y) = PDF_func(x=x/y, μ=μ, flavor=flavor)
    fg(y) = PDF_func(x=x/y, μ=μ, flavor="g")

    integrand_qq(y) = (1/y)*(
        C1_PDF_qq_Reg(y)*fq(y)
    )
    constant_qq = C1_PDF_qq_Delta*fq(1)

    integrand_qg(y) = (1/y)*(
        C1_PDF_qg_Reg(y)*fg(y)
    )

    integrand(y) = integrand_qq(y) + integrand_qg(y)
    constant = constant_qq

    #integrated, err = quadgk(integrand, x, 1.0, rtol=rtol_TMD)
    integrated = FastGK.integrate(f=integrand, a=x, b=1.0, rtol=rtol_TMD)

    return (integrated + constant)
end

# NNLO

# qq

C2_PDF_qq_Reg(x,nf) = (
    (4*(6*log(x)*(13 + 8*log1p(-x)*(-1 + x) - 3*x) + 6*log(x)^2*(3 - 2*x) + 6*nf*(-1 + x) + 13*π^2*(-1 + x) + 4*log(x)^3*(1 + x) + 6*(22 + (-22 + 5*log1p(-x))*x) + 
    12*(-1 + x)*reli(2, 1 - x)))/27
)

C2_PDF_qq_Delta(nf) = (
    (14568 - 656*nf - 1809*π^2 + 90*nf*π^2 + 60*π^4 + 252*(-33 + 2*nf)*z3)/243
)

C2_PDF_qq_D0(x,nf) = (
    (-2*(1 + x^2)*(3*log(x)*(3*log(x)*(21 + 6*log(x) - 2*nf) - 4*(-66 + 5*nf + 12*π^2)) - 4*(-606 + 28*nf + 603*z3) + 864*log(x)*reli(2, 1 - x) + 648*log(x)*reli(2, x) + 
    360*reli(3, 1 - x) + 144*reli(3, x)))/81
)

C2_PDF_qq_D0_x1(nf) = (
    (4*(-2424 + 112*nf + 2268*z3))/81
)

C2_PDF_qq_D1(x,nf) = (
    (-80*(1 + x^2)*(2*log(x)^2 - reli(2, 1 - x)))/9
)

C2_PDF_qq_D1_x1(nf) = 0

C2_PDF_qq_D2(x,nf) = (
    (32*log(x)*(1 + x^2))/9
)

C2_PDF_qq_D2_x1(nf) = 0

# qg

C2_PDF_qg_Reg(x,nf) = (
    -131/3 + 172/(9*x) + 93*x + (5*log1p(-x)*(3 - 4*x)*x)/3 + (6*log1p(-x)^2 - π^2)*(-1 + x)*x + (4*(-6*log1p(-x)*(log1p(-x) - 2*log(x)) + (3 + 2*log(x))*π^2)*(-1 + x)*x)/9 - 
    (730*x^2)/9 + log(x)^3*(1 + 2*x) + (2*log(x)*(24 + 2*π^2 + 3*(15 - 8*x)*x))/9 + (log(x)^2*(1 + 4*(3 - 2*x)*x))/3 - (2*log(x)^3*(1 - 2*x + 4*x^2))/9 + 
    (2*log(x)*(21 - 6*(5 + π^2)*x + 68*x^2))/3 - (log(x)^2*(3 + 4*x*(-3 + 11*x)))/2 + (4*(-2 + x*(3 + x*(-9 + 14*x)))*reli(2, 1 - x))/x - 6*x*(1 + x)*reli(2, 1 - x^2) + 
    (1 + 2*(-1 + x)*x)*(log1p(-x)^3 - 6*log1p(-x)*log(x)^2 - 9*z3 - 6*(log1p(-x) + log(x))*reli(2, 1 - x) + 6*reli(3, 1 - x) - 12*reli(3, x)) + 
    (2*(1 + 2*(-1 + x)*x)*((-2*log1p(-x)^3)/3 + 32*z3 + ((log1p(-x) - log(x))*(π^2 + 6*reli(2, 1 - x) - 6*reli(2, x)))/3 - 4*reli(3, 1 - x) - 4*reli(3, x)))/3 + 
    (1 + 2*x*(1 + x))*(log(x)*(3*log(1+x)*(-log(1+x) + log(x)) + π^2) - 3*z3 - 3*log(x)*reli(2, x^2) + 3*reli(3, x^2) + 6*reli(3, (1 + x)^(-1)) - 6*reli(3, x/(1 + x)))
)

# qq'

C2_PDF_qqp_Reg(x,nf) = (
    (2*(344 + 9*(-70 + log(x)*(28 + log(x)*(-3 + 2*log(x))))*x + 9*(62 + log(x)*(-40 + log(x)*(-3 + 2*log(x))))*x^2 - 8*(34 - 48*log(x) + 9*log(x)^2)*x^3 + 
    72*(-1 + x)*(2 + x*(-1 + 2*x))*reli(2, 1 - x)))/(81*x)
)

# qqbar

C2_PDF_qqb_Reg(x,nf) = (
    (-2*(-30*(-1 + x) + log(x)*(6 + 22*x) - 16*reli(2, 1 - x) + 4*(1 + x)*reli(2, 1 - x^2) - 
    (2*(1 + x^2)*(6*log(1+x)^2*log(x) - 6*log(1+x)*log(x)^2 + log(x)^3 + 6*z3 + 6*log(x)*reli(2, x^2) - 6*reli(3, x^2) - 12*reli(3, (1 + x)^(-1)) + 12*reli(3, x/(1 + x))))/
    (3*(1 + x))))/9
)

function F2_q_to_h(; x::Float64, μ::Float64, flavor::String)

    fq(y) = PDF_func(x=x/y, μ=μ, flavor=flavor)
    fq1 = PDF_func(x=x, μ=μ, flavor=flavor)
    fqb(y) = PDF_func(x=x/y, μ=μ, flavor=anti_flavor(flavor))
    fg(y) = PDF_func(x=x/y, μ=μ, flavor="g")

    nf = nf_func(μ)

    integrand_qq(y) = (
           1/y*C2_PDF_qq_Reg(y,nf)*fq(y)
        + (1/y*C2_PDF_qq_D0(y,nf)*fq(y) - C2_PDF_qq_D0_x1(nf)*fq1)/(1-y)
        + (1/y*C2_PDF_qq_D1(y,nf)*fq(y) - C2_PDF_qq_D1_x1(nf)*fq1)*log1p(-y)/(1-y)
        + (1/y*C2_PDF_qq_D2(y,nf)*fq(y) - C2_PDF_qq_D2_x1(nf)*fq1)*log1p(-y)^2/(1-y)
    )
    constant_qq = fq1 * (
          C2_PDF_qq_Delta(nf) 
        + C2_PDF_qq_D0_x1(nf) * log1p(-x)
        + C2_PDF_qq_D1_x1(nf) * log1p(-x)^2/2
        + C2_PDF_qq_D2_x1(nf) * log1p(-x)^3/3
    ) 

    integrand_qg(y) = (1/y)*(
        C2_PDF_qg_Reg(y,nf)*fg(y)
    )
    
    fq_sum(y) = sum(PDF_func(x=x/y, μ=μ, flavor=f) for f in active_flavors)
    integrand_qqp(y) = (1/y)*(
        C2_PDF_qqp_Reg(y,nf)*fq_sum(y)
    )

    integrand_qqb(y) = (1/y)*(
        C2_PDF_qqb_Reg(y,nf)*fqb(y)
    )

    integrand(y) = (
          integrand_qq(y) 
        + integrand_qg(y)
        + integrand_qqp(y)
        + integrand_qqb(y)
    )
    constant = constant_qq

    #integrated, err =quadgk(integrand, x, 1.0, rtol=rtol_TMD)
    integrated = FastGK.integrate(f=integrand, a=x, b=1.0, rtol=rtol_TMD)

    return (integrated + constant)
end

# -----------------------------------------------------------------------

function TMDPDF_raw_func(; b::Float64, x::Float64, order::Int64=2)

    #bstar = Main.bstar_func(b=b)

    μ_safe = b0/b#star

    as = αs_func(μ_safe)/(4π)

    flavors = ["u", "ub", "d", "db", "s", "sb", "c", "cb", "b", "bb"]

    f = Dict{String, Float64}()
    f0 = Dict{String, Float64}()
    f1 = Dict{String, Float64}()
    f2 = Dict{String,Float64}()

    for name in flavors

        f0[name] = PDF_func(x = x, μ = μ_safe, flavor = name)
        f1[name] = f0[name] + as*F1_q_to_h(x=x, μ=μ_safe, flavor=name)
        f2[name] = f1[name] + as^2*F2_q_to_h(x=x, μ=μ_safe, flavor=name)

    end
        
    if order == 0
        f = f0
    elseif order == 1
        f = f1
    elseif order == 2
        f = f2
    else
        error("Unsupported Order")
    end   

    fu = f["u"]
    fub = f["ub"]
    fd = f["d"]
    fdb = f["db"]
    fs = f["s"]
    fsb = f["sb"]
    fc = f["c"]
    fcb = f["cb"]
    fb = f["b"]
    fbb = f["bb"]

    return fu, fub, fd, fdb, fs, fsb, fc, fcb, fb, fbb
end

function TMDPDF_func(; b::Float64, x::Float64, Q::Float64)

    bstar = Main.bstar_func(b=b,Q=Q)
    μi = b0/bstar
    μf = Q
    ζi = (b0/bstar)^2
    ζf = Q^2

    fu, fub, fd, fdb, fs, fsb, fc, fcb, fb, fbb = TMDPDF_raw_func(b=bstar, x=x, order=XLO)
    SP = PertSudakov_numerical(μi=μi, ζi=ζi, μf=μf, ζf=ζf)
    #SP = PertSudakov_numerical(b=bstar, μi=μi, ζi=ζi, μf=μf, ζf=ζf, order=XLL)
    #SNP = NP_f_func_CFR25(x=x, b=b, Q=Q)

    SNP = NP_f_func_CFR25(x=x, b=b, Q=Q)

    return (fu, fub, fd, fdb, fs, fsb, fc, fcb, fb, fbb) .* (SP*SNP)
end

function xTMDPDF_raw_func(x, b)

    if if_grid == true
        return xTMDPDF_raw_grid(x,b)
    else
        return x .* TMDPDF_raw_func(b=b, x=x) 
    end
end

function μ_evolution_integrand(; μ::Float64, ζf::Float64, order::Int64)

    Γ = Γ_func(μ=μ, order=order)
    γV = γV_func(μ=μ, order=order)

    return (Γ*log(ζf/μ^2) + γV)/μ
end

function PertSudakov_numerical(; μi::Float64, ζi::Float64, μf::Float64, ζf::Float64)

    integrand(μ) = μ_evolution_integrand(μ=μ, ζf=ζf, order=XLL)
    integral, err = quadgk(integrand, μi, μf, rtol=rtol_Sudakov)

    CS = CS_func(μ = μi, order = XLL)

    return exp(-integral) * (ζf/ζi)^(-CS)
end

function PertSudakov_func(b, Q)
    if if_grid == true
        return PertSudakov_grid(b,Q)
    else
        μi = b0/b
        ζi = (b0/b)^2
        μf = Q
        ζf = Q^2
        return PertSudakov_numerical(μi=μi, ζi=ζi, μf=μf, ζf=ζf) # C wrapper
    end
end

function TMDPDF_pert_func(; b::Float64, x::Float64, Q::Float64)

    bstar = bstar_func(b=b,Q=Q) 

    TMDPDF_raw = xTMDPDF_raw_func(x, bstar)

    SP = PertSudakov_func(bstar, Q)

    return TMDPDF_raw .* (SP/x)
end

function TMD_per_nucleon_func(fu, fub, fd, fdb, isoscalarity)  

    if isoscalarity == 1.0
        return fu, fub, fd, fdb
    else
        ZdA = abs(isoscalarity)
        NdA = 1 - ZdA

        fu_mix = ZdA*fu+NdA*fd
        fub_mix = ZdA*fub+NdA*fdb
        fd_mix = ZdA*fd+NdA*fu
        fdb_mix = ZdA*fdb+NdA*fub               

        if isoscalarity > 0.0
            return fu_mix, fub_mix, fd_mix, fdb_mix
        elseif isoscalarity < 0.0
            return fub_mix, fu_mix, fdb_mix, fd_mix
        else
            return error("For neutron/anti-neutron pass a very small non-zero isoscalarity, like eps = 1e-6/-1e-6")
        end
    end 
end

function TMDPDF_func(; b::Float64, x::Float64, Q::Float64, isoscalarity::Float64=1.0)

    fu, fub, fd, fdb, fs, fsb, fc, fcb, fb, fbb = TMDPDF_pert_func(b=b, x=x, Q=Q)

    if if_grid == true
        SNP = NP_f_grid(b, Q)
    else
        SNP = NP_f_func(b, Q)
    end

    Fu, Fub, Fd, Fdb = 
    TMD_per_nucleon_func(fu, fub, fd, fdb, isoscalarity)

    (Fs, Fsb, Fc, Fcb, Fb, Fbb) = (fs, fsb, fc, fcb, fb, fbb) 

    return (Fu, Fub, Fd, Fdb, Fs, Fsb, Fc, Fcb, Fb, Fbb) .* SNP 
end

function TMDPDF_kt_func(; kt::Float64, x::Float64, Q::Float64)
    integrand(b) = begin
        vals = TMDPDF_func(b=b, x=x, Q=Q)  
        s = b * besselj0(b*kt)/(2*π)
        s .* SVector{10,Float64}(vals)      
    end
    out, _ = quadgk(integrand, 0.001, 30.0, rtol=1e-3)
    return out
end