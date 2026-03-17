*! cs_region.ado  v1.0.0  2026-03-17
*! Regional fragility summary for causalspline
*!
*! Syntax:
*!   cs_region, a(#) b(#) [type(curvature_ratio|inverse_slope)]

program define cs_region, rclass
    version 14.0

    syntax , A(real) B(real) [ Type(string) ]

    if "`type'" == "" local type "curvature_ratio"
    if !inlist("`type'", "curvature_ratio", "inverse_slope") {
        di as error "type() must be: curvature_ratio  inverse_slope"
        exit 198
    }
    if `a' >= `b' {
        di as error "a() must be less than b()"
        exit 198
    }
    if "`e(cmd)'" != "causalspline" {
        di as error "cs_region requires causalspline to be run first"
        exit 301
    }

    // -- Clamp interval to observed treatment support --------------------------
    local t_min = e(t_min)
    local t_max = e(t_max)
    local a0 = max(`a', `t_min')
    local b0 = min(`b', `t_max')

    if `a0' >= `b0' {
        di as error "Requested interval does not overlap observed support " ///
            "[" `t_min' ", " `t_max' "]"
        exit 198
    }

    // -- Get fragility over full grid ------------------------------------------
    qui cs_fragility, type(`type')
    tempname ft frag
    mat `ft'   = r(frag_t)
    mat `frag' = r(fragility)
    local ng = e(evalgrid)

    // -- Extract points in [a0, b0] --------------------------------------------
    // Collect into Mata for integration
    mata: _cs_region_integrate("`ft'", "`frag'", `ng', `a0', `b0')
    local integral = r(integral)
    local average  = r(average)
    local npoints  = r(npoints)

    // -- Display ---------------------------------------------------------------
    di as text _n "{hline 55}"
    di as text " Regional Fragility Summary"
    di as text "{hline 55}"
    di as text "  Interval       : [" %6.3f `a0' ", " %6.3f `b0' "]"
    di as text "  Type           : " as result "`type'"
    di as text "  Grid points    : " as result `npoints'
    di as text "{hline 55}"
    di as text "  Integral fragility : " as result %10.5f `integral'
    di as text "  Average fragility  : " as result %10.5f `average'
    di as text "{hline 55}"

    // -- Interpretation guide -------------------------------------------------
    di as text _n "  Interpretation (curvature_ratio):"
    di as text "  < 0.5  : stable region (low curvature relative to slope)"
    di as text "  0.5-2  : moderate structural change"
    di as text "  > 2    : high fragility (threshold / turning point)"

    // -- rreturn ---------------------------------------------------------------
    return scalar a                 = `a0'
    return scalar b                 = `b0'
    return scalar integral_fragility = `integral'
    return scalar average_fragility  = `average'
    return scalar npoints            = `npoints'
    return local  type               = "`type'"

end


mata:

void _cs_region_integrate(string scalar ftname, string scalar fragname,
                           real scalar ng, real scalar a0, real scalar b0)
{
    real matrix ft, fv
    ft = st_matrix(ftname)
    fv = st_matrix(fragname)

    // Select points in [a0, b0] with finite values
    real colvector x, y, keep
    keep = (ft :>= a0) :& (ft :<= b0) :& (fv :!= .)
    x    = select(ft, keep)
    y    = select(fv, keep)

    real scalar npts, integral, avg
    npts = rows(x)

    if (npts < 2) {
        integral = .
        avg      = .
    }
    else {
        // Trapezoidal rule
        real colvector dx, ymid
        dx     = x[2::npts] - x[1::npts-1]
        ymid   = (y[1::npts-1] + y[2::npts]) / 2
        integral = sum(dx :* ymid)
        avg      = integral / (x[npts] - x[1])
    }

    st_numscalar("r(integral)", integral)
    st_numscalar("r(average)",  avg)
    st_numscalar("r(npoints)",  npts)
}

end
