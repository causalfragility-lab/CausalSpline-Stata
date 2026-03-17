*! cs_plot.ado  v1.0.0  2026-03-17
*! Plot dose-response curve from causalspline
*! Stata 14.1 compatible
*!
*! Syntax:
*!   cs_plot [, title(string) saving(filename)]

program define cs_plot
    version 14.0

    syntax [, title(string) SAVing(string) ]

    if "`e(cmd)'" != "causalspline" {
        di as error "cs_plot requires causalspline to be run first"
        exit 301
    }

    local ng     = e(evalgrid)
    local method = upper(e(method))
    local tvar   = e(treatment)
    if "`title'" == "" local title "Causal Dose-Response [`method']"

    preserve
        qui drop _all
        qui set obs `ng'
        qui gen double t     = .
        qui gen double est   = .
        qui gen double ci_lo = .
        qui gen double ci_hi = .

        tempname ct ce clo chi
        mat `ct'  = e(curve_t)
        mat `ce'  = e(curve_est)
        mat `clo' = e(curve_lo)
        mat `chi' = e(curve_hi)

        forval j = 1/`ng' {
            qui replace t     = `ct'[`j',1]  in `j'
            qui replace est   = `ce'[`j',1]  in `j'
            qui replace ci_lo = `clo'[`j',1] in `j'
            qui replace ci_hi = `chi'[`j',1] in `j'
        }

        local ymean = e(t_mean)

        twoway  (rarea ci_lo ci_hi t, lwidth(none) color(ltblue))  ///
                (line est t, lcolor(navy) lwidth(medthick)),        ///
            yline(`ymean', lpattern(dash) lcolor(gs10))             ///
            xtitle("Treatment (`tvar')")                            ///
            ytitle("E[Y(t)]")                                       ///
            title("`title'")                                        ///
            subtitle("95% CI shaded  |  Dashed = marginal mean")   ///
            legend(order(2 "E[Y(t)]" 1 "95% CI")                   ///
                   position(6) cols(2) size(small))                 ///
            scheme(s2color)

        if "`saving'" != "" {
            graph save "`saving'", replace
            di as text "  Plot saved: " as result "`saving'"
        }
    restore
end
