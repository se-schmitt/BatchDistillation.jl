# Reboiler - bottom stage. Receives reflux liquid from above, sends vapor up.
# Heat duty Q_reb drives vaporization. Holdup depletes as distillate is removed.
function Reboiler(; name, Q_reb=500.0, τ_V=1e-3)
    @named liq_in  = MaterialStream()
    @named vap_out = MaterialStream()

    sts = @variables begin
        n1(t); n2(t); n3(t)
        H_stage(t)
        T(t); β(t)
        x1(t); x2(t); x3(t)
        y1(t); y2(t); y3(t)
        n_tot(t); h_spec(t); HL(t); HV(t)
        hL(t); hV(t)
        V_out(t)
    end
    pars = @parameters begin
        Q_reb_p = Q_reb
        τ_V_p = τ_V
    end

    eqs = [
        D(n1) ~ liq_in.F*liq_in.z1 - V_out*y1,
        D(n2) ~ liq_in.F*liq_in.z2 - V_out*y2,
        D(n3) ~ liq_in.F*liq_in.z3 - V_out*y3,
        D(H_stage) ~ liq_in.F*liq_in.h - V_out*hV + Q_reb_p,

        n_tot  ~ n1 + n2 + n3,
        h_spec ~ H_stage / n_tot,
        HL     ~ (1 - β) * n_tot,
        HV     ~ β * n_tot,

        T  ~ flash_T(p_set,  H_stage, n1, n2, n3),
        β  ~ flash_β(p_set,  H_stage, n1, n2, n3),
        x1 ~ flash_x1(p_set, H_stage, n1, n2, n3),
        x2 ~ flash_x2(p_set, H_stage, n1, n2, n3),
        x3 ~ flash_x3(p_set, H_stage, n1, n2, n3),
        y1 ~ flash_y1(p_set, H_stage, n1, n2, n3),
        y2 ~ flash_y2(p_set, H_stage, n1, n2, n3),
        y3 ~ flash_y3(p_set, H_stage, n1, n2, n3),

        hL ~ h_liq(p_set, T, x1, x2, x3),
        hV ~ h_vap(p_set, T, y1, y2, y3),

        V_out ~ HV / τ_V_p,

        vap_out.F  ~ V_out,
        vap_out.z1 ~ y1, vap_out.z2 ~ y2, vap_out.z3 ~ y3,
        vap_out.T  ~ T,  vap_out.p ~ p_set,  vap_out.h ~ hV,
    ]
    compose(ODESystem(eqs, t, sts, pars; name), [liq_in, vap_out])
end