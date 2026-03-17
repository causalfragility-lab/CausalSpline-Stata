*! cs_gradient.ado  v1.0.0  2026-03-17
*! Compute first and second derivatives of causalspline curve
*!
*! Syntax:
*!   cs_gradient [, savegradient(filename) verbose]
*!
*! Uses e() results from causalspline

program define cs_gradient, rclass
    version 14.0

    syntax [, SAVEgradient(string) VERBose ]

    // -- Check that causalspline was run --------------------------------------
    if "`e(cmd)'" != "causalspline" {
        di as error "cs_gradient requires causalspline to be run first"
        exit 301
    }

    // -- Retrieve stored curve ------------------------------------------------
    tempname ct ce cse
    mat `ct'  = e(curve_t)
    mat `ce'  = e(curve_est)
    mat `cse' = e(curve_se)

    local ng = e(evalgrid)

    // -- Compute numerical derivatives ----------------------------------------
    tempname d1 d2
    mat `d1' = J(`ng', 1, .)
    mat `d2' = J(`ng', 1, .)

    // Central differences for interior; forward/backward at boundaries
    forval j = 1/`ng' {
        if `j' == 1 {
            // Forward difference
            if `ng' >= 2 {
                local h = `ct'[2,1] - `ct'[1,1]
                mat `d1'[`j',1] = (`ce'[2,1] - `ce'[1,1]) / `h'
                mat `d2'[`j',1] = .
            }
        }
        else if `j' == `ng' {
            // Backward difference
            local h = `ct'[`ng',1] - `ct'[`=`ng'-1',1]
            mat `d1'[`j',1] = (`ce'[`ng',1] - `ce'[`=`ng'-1',1]) / `h'
            mat `d2'[`j',1] = .
        }
        else {
            // Central difference
            local jm1 = `j' - 1
            local jp1 = `j' + 1
            local hl = `ct'[`j',1]  - `ct'[`jm1',1]
            local hr = `ct'[`jp1',1] - `ct'[`j',1]

            // First derivative (central)
            mat `d1'[`j',1] = (`ce'[`jp1',1] - `ce'[`jm1',1]) / (`hl' + `hr')

            // Second derivative
            local num = (`ce'[`jm1',1]/(`hl'*(`hl'+`hr'))) ///
                      - (`ce'[`j',1] /(`hl'*`hr'))          ///
                      + (`ce'[`jp1',1]/(`hr'*(`hl'+`hr')))
            mat `d2'[`j',1] = 2 * `num'
        }
    }

    // -- Display ---------------------------------------------------------------
    di as text _n "{hline 75}"
    di as text " Gradient Curve (Derivatives of Dose-Response Function)"
    di as text "{hline 75}"
    di as text %10s "t" %13s "E[Y(t)]" %11s "SE" ///
               %13s "dE/dt" %13s "d2E/dt2"
    di as text "{hline 75}"

    foreach j in 1 11 26 51 76 91 `ng' {
        if `j' >= 1 & `j' <= `ng' {
            local d2val = `d2'[`j',1]
            if `d2val' == . {
                di as result %10.3f `ct'[`j',1]  ///
                             %13.4f `ce'[`j',1]  ///
                             %11.4f `cse'[`j',1] ///
                             %13.4f `d1'[`j',1]  ///
                             %13s   "."
            }
            else {
                di as result %10.3f `ct'[`j',1]  ///
                             %13.4f `ce'[`j',1]  ///
                             %11.4f `cse'[`j',1] ///
                             %13.4f `d1'[`j',1]  ///
                             %13.6f `d2'[`j',1]
            }
        }
    }
    di as text "{hline 75}"
    di as text "First derivative  = marginal causal effect dE[Y(t)]/dt"
    di as text "Second derivative = curvature (acceleration of effect)"

    // -- Save if requested -----------------------------------------------------
    if "`savegradient'" != "" {
        preserve
            qui drop _all
            qui set obs `ng'
            qui gen double t                 = .
            qui gen double estimate          = .
            qui gen double se                = .
            qui gen double derivative        = .
            qui gen double second_derivative = .
            forval j = 1/`ng' {
                qui replace t                 = `ct'[`j',1]  in `j'
                qui replace estimate          = `ce'[`j',1]  in `j'
                qui replace se                = `cse'[`j',1] in `j'
                qui replace derivative        = `d1'[`j',1]  in `j'
                qui replace second_derivative = `d2'[`j',1]  in `j'
            }
            label var t                 "Treatment"
            label var estimate          "E[Y(t)]"
            label var se                "SE of E[Y(t)]"
            label var derivative        "First derivative dE/dt"
            label var second_derivative "Second derivative d2E/dt2"
            qui save "`savegradient'", replace
        restore
        di as text "  Gradient saved: " as result "`savegradient'"
    }

    // -- rreturn ---------------------------------------------------------------
    return mat grad_t  = `ct'
    return mat grad_mu = `ce'
    return mat grad_se = `cse'
    return mat grad_d1 = `d1'
    return mat grad_d2 = `d2'

end
