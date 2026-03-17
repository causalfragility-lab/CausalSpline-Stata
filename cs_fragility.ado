*! cs_fragility.ado  v1.0.1  2026-03-17  Stata 14.1 compatible ASCII only
*! Geometric fragility curve for causalspline - with dual panel plot
*!
*! Syntax:
*!   cs_fragility [, type(curvature_ratio|inverse_slope)
*!                   savefragility(filename) saving(filename) noplot]

program define cs_fragility, rclass
    version 14.0

    syntax [, Type(string) SAVEFragility(string) SAVing(string) NOPlot ]

    if "`type'" == "" local type "curvature_ratio"
    if !inlist("`type'", "curvature_ratio", "inverse_slope") {
        di as error "type() must be: curvature_ratio  inverse_slope"
        exit 198
    }
    if "`e(cmd)'" != "causalspline" {
        di as error "cs_fragility requires causalspline to be run first"
        exit 301
    }

    // Get derivatives
    qui cs_gradient
    tempname ct ce cse d1 d2
    mat `ct'  = r(grad_t)
    mat `ce'  = r(grad_mu)
    mat `cse' = r(grad_se)
    mat `d1'  = r(grad_d1)
    mat `d2'  = r(grad_d2)

    local ng = e(evalgrid)

    // Adaptive eps = 0.05 * median|d1|
    mata: _cs_frag_eps("`d1'", `ng')
    local eps = r(eps)

    // Compute fragility
    tempname frag frag_norm hiflag zone_num
    mat `frag'      = J(`ng', 1, .)
    mat `frag_norm' = J(`ng', 1, .)

    forval j = 1/`ng' {
        local d1j = `d1'[`j',1]
        local d2j = `d2'[`j',1]
        local sej = `cse'[`j',1]

        if `d1j' != . & `d2j' != . {
            local ad1   = abs(`d1j')
            local ad2   = abs(`d2j')
            local denom = `ad1' + `eps'

            if "`type'" == "curvature_ratio" {
                mat `frag'[`j',1] = `ad2' / `denom'
            }
            else {
                mat `frag'[`j',1] = 1 / `denom'
            }
            if `sej' != . & `sej' > 0 {
                mat `frag_norm'[`j',1] = `frag'[`j',1] / `sej'
            }
        }
    }

    // Quantile thresholds
    mata: _cs_frag_quantiles("`frag'", `ng')
    local q50 = r(q50)
    local q75 = r(q75)

    // Assign zones
    mat `zone_num' = J(`ng', 1, .)
    mat `hiflag'   = J(`ng', 1, 0)

    forval j = 1/`ng' {
        local fv = `frag'[`j',1]
        if `fv' != . {
            if `fv' > `q75' {
                mat `zone_num'[`j',1] = 3
                mat `hiflag'[`j',1]   = 1
            }
            else if `fv' > `q50' {
                mat `zone_num'[`j',1] = 2
            }
            else {
                mat `zone_num'[`j',1] = 1
            }
        }
    }

    // Display table
    di as text " "
    di as text "{hline 75}"
    di as text " Fragility Curve   Type: " as result "`type'"
    di as text "{hline 75}"
    di as text %10s "t" %13s "E[Y(t)]" %13s "Fragility" ///
               %10s "Frag_norm" %10s "Zone"
    di as text "{hline 75}"

    foreach j in 1 11 26 51 76 91 `ng' {
        if `j' >= 1 & `j' <= `ng' {
            local fv = `frag'[`j',1]
            local zv = `zone_num'[`j',1]
            local zlab = cond(`zv'==3,"high",cond(`zv'==2,"moderate","low"))
            if `fv' != . {
                di as result %10.3f `ct'[`j',1] ///
                             %13.4f `ce'[`j',1] ///
                             %13.5f `fv'          ///
                             %10.5f `frag_norm'[`j',1] ///
                             %10s   "`zlab'"
            }
        }
    }
    di as text "{hline 75}"
    di as text "q50 = " %7.5f `q50' "  q75 = " %7.5f `q75'
    di as text "Zones: low < q50 <= moderate < q75 <= high"

    // Save fragility dataset if requested
    if "`savefragility'" != "" {
        preserve
            qui drop _all
            qui set obs `ng'
            qui gen double t              = .
            qui gen double estimate       = .
            qui gen double se             = .
            qui gen double d1             = .
            qui gen double d2             = .
            qui gen double fragility      = .
            qui gen double fragility_norm = .
            qui gen byte   high_fragility = .
            qui gen byte   zone           = .
            forval j = 1/`ng' {
                qui replace t              = `ct'[`j',1]        in `j'
                qui replace estimate       = `ce'[`j',1]        in `j'
                qui replace se             = `cse'[`j',1]       in `j'
                qui replace d1             = `d1'[`j',1]        in `j'
                qui replace d2             = `d2'[`j',1]        in `j'
                qui replace fragility      = `frag'[`j',1]      in `j'
                qui replace fragility_norm = `frag_norm'[`j',1] in `j'
                qui replace high_fragility = `hiflag'[`j',1]    in `j'
                qui replace zone           = `zone_num'[`j',1]  in `j'
            }
            label define zonelbl 1 "low" 2 "moderate" 3 "high"
            label values zone zonelbl
            qui save "`savefragility'", replace
        restore
        di as text "  Fragility saved: " as result "`savefragility'"
    }

    // Dual panel plot (mirrors R plot.fragility_curve)
    if "`noplot'" == "" {
        preserve
            qui drop _all
            qui set obs `ng'
            qui gen double t         = .
            qui gen double estimate  = .
            qui gen double ci_lo     = .
            qui gen double ci_hi     = .
            qui gen double fragility = .
            qui gen byte   hi_flag   = .

            tempname clo chi
            mat `clo' = e(curve_lo)
            mat `chi' = e(curve_hi)

            forval j = 1/`ng' {
                qui replace t         = `ct'[`j',1]   in `j'
                qui replace estimate  = `ce'[`j',1]   in `j'
                qui replace ci_lo     = `clo'[`j',1]  in `j'
                qui replace ci_hi     = `chi'[`j',1]  in `j'
                qui replace fragility = `frag'[`j',1] in `j'
                qui replace hi_flag   = `hiflag'[`j',1] in `j'
            }

            // Find high fragility x-range for shading
            qui sum t if hi_flag == 1
            local hf_n = r(N)

            // Top panel: dose-response curve
            if `hf_n' > 0 {
                local hf_min = r(min)
                local hf_max = r(max)
                twoway (rarea ci_lo ci_hi t, lwidth(none) color(ltblue)) ///
                       (line estimate t, lcolor(navy) lwidth(medthick)), ///
                    xline(`hf_min' `hf_max', lpattern(dash) lcolor(red) lwidth(thin)) ///
                    xtitle("") ytitle("E[Y(t)]") ///
                    title("Dose-response with fragility regions") ///
                    subtitle("Shaded high fragility (top 25%) | Type: `type'") ///
                    legend(off) name(cs_top, replace) nodraw
            }
            else {
                twoway (rarea ci_lo ci_hi t, lwidth(none) color(ltblue)) ///
                       (line estimate t, lcolor(navy) lwidth(medthick)), ///
                    xtitle("") ytitle("E[Y(t)]") ///
                    title("Dose-response with fragility regions") ///
                    subtitle("Type: `type'") ///
                    legend(off) name(cs_top, replace) nodraw
            }

            // Bottom panel: fragility curve
            twoway (line fragility t, lcolor(navy) lwidth(medthick)), ///
                yline(`q75', lpattern(dash) lcolor(red) lwidth(thin)) ///
                yline(`q50', lpattern(dot)  lcolor(gs8) lwidth(thin)) ///
                xtitle("Treatment (T)") ytitle("Fragility") ///
                note("Red dashed = 75th pct (high) | Grey dotted = 50th pct (moderate)") ///
                legend(off) name(cs_bot, replace) nodraw

            // Combine panels
            graph combine cs_top cs_bot, cols(1) ///
                title("CausalSpline Fragility Diagnostics") ///
                scheme(s2color)

            if "`saving'" != "" {
                graph save "`saving'", replace
                di as text "  Plot saved: " as result "`saving'"
            }
        restore
    }

    // rreturn
    return scalar q50  = `q50'
    return scalar q75  = `q75'
    return scalar eps  = `eps'
    return local  type = "`type'"
    return mat frag_t    = `ct'
    return mat frag_mu   = `ce'
    return mat frag_d1   = `d1'
    return mat frag_d2   = `d2'
    return mat fragility = `frag'
    return mat frag_norm = `frag_norm'
    return mat hiflag    = `hiflag'
    return mat zone      = `zone_num'
end


mata:

void _cs_frag_eps(string scalar d1name, real scalar ng)
{
    real matrix  ad
    real scalar  med_slope, eps

    ad = st_matrix(d1name)
    ad = abs(ad)
    ad = select(ad, ad :!= .)
    if (rows(ad) > 0) {
        ad        = sort(ad, 1)
        med_slope = ad[ceil(rows(ad)/2), 1]
    }
    else {
        med_slope = 1
    }
    eps = 0.05 * max((med_slope, 1e-8))
    st_numscalar("r(eps)", eps)
}


void _cs_frag_quantiles(string scalar fragname, real scalar ng)
{
    real matrix  fv
    real scalar  n, q50, q75

    fv = st_matrix(fragname)
    fv = select(fv, fv :!= .)
    fv = sort(fv, 1)
    n  = rows(fv)
    if (n > 0) {
        q50 = fv[ceil(0.50 * n), 1]
        q75 = fv[ceil(0.75 * n), 1]
    }
    else {
        q50 = 0
        q75 = 0
    }
    st_numscalar("r(q50)", q50)
    st_numscalar("r(q75)", q75)
}

end
