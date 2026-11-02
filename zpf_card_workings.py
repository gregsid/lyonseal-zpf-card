# -*- coding: utf-8 -*-
"""
zpf_card_workings.py -- the numbers behind "Your six steps, filled in for the ZPF Array".

Plain Python 3 (standard library only). Change the inputs in INPUTS below to your own design
values and run it:   python zpf_card_workings.py

Everything is an estimate from textbook physics. Each section names the formula it uses, so
anyone can check it by hand. The script first tests itself against cases with known answers
(ideal mirrors, perfect conductors) and stops if any test fails.
"""
import math, cmath

# ------------------------------------------------------------------ INPUTS (edit these)
INPUTS = dict(
    f_drive=1.0e9,          # wall drive frequency, Hz
    stroke=1.0e-9,          # wall motion AMPLITUDE, m (peak-to-peak travel is twice this)
    gap=200e-9,             # gap between the walls, m
    wavelength=1.55e-6,     # probe light, m
    slot_length=100e-6,     # length of the slot the light runs along, m
    eps_metal=complex(-129, 3.3),   # wall permittivity at the probe wavelength (silver, Johnson & Christy);
                                    # replace with your film's measured value
    P_in=10.0,              # electrical input power, W
    claimed_thrust=(3e-3, 15e-3),   # N (pre-seed page)
    P_detected=1e-4,        # probe power reaching the detector, W
    t_integrate=3 * 3600.0, # integration time, s
    stroke_cal=0.01,        # how well the stroke is known (fraction); 0.01 = 1 %
    dK_example=0.05,        # a vacuum change to scale against (5 %)
    v_cap=50.0,             # material limit on wall speed, m/s (Dodonov 2010, sec. 3.1)
    wilson_v_over_c=0.05,   # Wilson et al. 2011's effective boundary speed
    T_fridge=0.010,         # K, for the cold-circuit option
    wall=dict(length=100e-6, height=1e-6, thickness=100e-9, density=19300.0),  # one gold wall, for its mass
    Q_mech=1000.0,          # mechanical quality factor of the driven wall (a guess; measure it)
    actuator_area=1e-8,     # area of thickness-mode actuator holding the stroke, m^2 (100 um x 100 um)
    coupler_dB=(1.4, 5.0),  # light into (or out of) a metal slot, per transition, measured devices
    fiber_dB=(1.4, 2.8),    # each fiber grating coupler
)

# Sound in the wall materials at room temperature: (density kg/m^3, longitudinal speed, shear speed, m/s).
# Gold: CRC Handbook (hard-drawn); silicon [100]: Ioffe database; AlN c-axis: c33 = 395 GPa, rho = 3260.
MATERIALS = {"gold": (19300.0, 3240.0, 1200.0), "silicon": (2329.0, 8433.0, 5843.0), "AlN": (3260.0, 11000.0, 6200.0)}

# ------------------------------------------------------------------ constants (CODATA)
c = 299792458.0
hbar = 1.054571817e-34
h = 2 * math.pi * hbar
kB = 1.380649e-23
eV = 1.602176634e-19
alpha = 7.2973525693e-3
lambda_C_reduced = 3.8615926796e-13  # hbar/(m_e c), m


# ------------------------------------------------------------------ the light the gap holds
def gap_mode(eps_m, d, lam, eps_d=1.0, mu_d=1.0):
    """Fundamental (symmetric TM) mode between two metal half-spaces with a gap of width d.
    Dispersion relation: eps_m * kd * tanh(kd*d/2) + eps_d * km = 0,
    kd^2 = beta^2 - k0^2 eps_d mu_d,  km^2 = beta^2 - k0^2 eps_m.  Returns beta (complex, 1/m)."""
    k0 = 2 * math.pi / lam

    def f(beta):
        kd = cmath.sqrt(beta * beta - k0 * k0 * eps_d * mu_d)
        km = cmath.sqrt(beta * beta - k0 * k0 * eps_m)
        if km.real < 0:
            km = -km
        return eps_m * kd * cmath.tanh(kd * d / 2) + eps_d * km

    delta = 1 / (k0 * cmath.sqrt(-eps_m))            # field penetration into the metal
    b0 = k0 * cmath.sqrt(eps_d * mu_d) * (1 + delta / d)
    b1 = b0 * (1 + 1e-4)
    f0, f1 = f(b0), f(b1)
    for _ in range(200):                              # secant iteration
        if f1 == f0:
            break
        b2 = b1 - f1 * (b1 - b0) / (f1 - f0)
        b0, f0, b1, f1 = b1, f1, b2, f(b2)
        if abs(b1 - b0) < 1e-14 * abs(b1):
            break
    if abs(f1) > 1e-6 * abs(eps_m) * k0:
        raise RuntimeError("mode solver did not converge")
    return b1


def mode_numbers(eps_m, d, lam):
    k0 = 2 * math.pi / lam
    beta = gap_mode(eps_m, d, lam)
    n_eff = beta.real / k0
    L_power = 1 / (2 * beta.imag)                      # distance to 1/e power
    hK = 1e-6                                          # vacuum change in the gap: eps and mu together, by K
    S = (gap_mode(eps_m, d, lam, 1 + hK, 1 + hK).real - beta.real) / (k0 * hK)
    hd = 1e-12                                         # gap change ("breathing")
    dn_dd = (gap_mode(eps_m, d + hd, lam).real - beta.real) / (k0 * hd)   # per metre
    return dict(n_eff=n_eff, L_power=L_power, S=S, dn_dd=dn_dd)


# ------------------------------------------------------------------ Casimir pressure, real metals
def gauss_legendre(n):
    """Nodes and weights on [0, 1]."""
    xs, ws = [], []
    for i in range(1, n + 1):
        x = math.cos(math.pi * (i - 0.25) / (n + 0.5))
        for _ in range(100):
            p0, p1 = 1.0, x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
            dp = n * (x * p1 - p0) / (x * x - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < 1e-15:
                break
        xs.append(0.5 * (x + 1)); ws.append(1.0 / ((1 - x * x) * dp * dp))
    return xs, ws


def casimir_pressure(d, wp_eV=None, gamma_eV=0.0, nt=48, nx=600, xmax=60.0):
    """Pressure between two identical half-spaces at T = 0 (Lifshitz formula), in Pa.
    wp_eV None means ideal mirrors. Otherwise a Drude metal, eps(i xi) = 1 + wp^2 / (xi (xi + gamma)).
    P = (hbar c / 2 pi^2) (1/2d)^4  Int_0^inf dx x^3  Int_0^1 dt  Sum_TE,TM  r^2 e^-x / (1 - r^2 e^-x),
    with q = x/2d the vacuum normal wavenumber and xi = c q t the imaginary frequency."""
    ts, wts = gauss_legendre(nt)
    wp = None if wp_eV is None else wp_eV * eV / hbar
    gm = gamma_eV * eV / hbar
    hx = xmax / nx
    total = 0.0
    for j in range(nx + 1):
        x = j * hx
        if x == 0:
            continue
        q = x / (2 * d)
        e = math.exp(-x)
        inner = 0.0
        for t, w in zip(ts, wts):
            if wp is None:
                r2s = (1.0, 1.0)
            else:
                xi = c * q * t
                eps = 1 + wp * wp / (xi * (xi + gm))
                s = math.sqrt(1 + (eps - 1) * t * t)
                r2s = (((1 - s) / (1 + s)) ** 2, ((eps - s) / (eps + s)) ** 2)
            inner += w * sum(r2 * e / (1 - r2 * e) for r2 in r2s)
        wsimp = 1 if j == nx else (4 if j % 2 else 2)
        total += wsimp * x ** 3 * inner
    total *= hx / 3
    return hbar * c / (2 * math.pi ** 2) * (1 / (2 * d)) ** 4 * total


def casimir_ideal(d):
    return math.pi ** 2 * hbar * c / (240 * d ** 4)


# ------------------------------------------------------------------ self-tests
def self_tests():
    # 1. ideal mirrors reproduce pi^2 hbar c / 240 d^4
    for d in (10e-9, 200e-9):
        num, exact = casimir_pressure(d), casimir_ideal(d)
        assert abs(num / exact - 1) < 2e-3, ("Casimir ideal limit", d, num, exact)
    # 2. a near-perfect conductor gives the free-space mode: n_eff -> 1, S -> 1, no breathing
    m = mode_numbers(complex(-1e8, 1.0), 200e-9, 1.55e-6)
    assert abs(m["n_eff"] - 1) < 1e-3 and abs(m["S"] - 1) < 1e-3 and abs(m["dn_dd"]) * 1e-9 < 1e-6, ("PEC limit", m)
    k0, b = 2 * math.pi / 1.55e-6, gap_mode(complex(-1e8, 1.0), 200e-9, 1.55e-6).real
    S_opp = (gap_mode(complex(-1e8, 1.0), 200e-9, 1.55e-6, 1 + 1e-6, 1 / (1 + 1e-6)).real - b) / (k0 * 1e-6)
    assert abs(S_opp) < 1e-3, ("PEC: eps*mu fixed must leave the phase alone", S_opp)
    # 3. the small-penetration formula n_eff ~ 1 + delta/d holds for silver to a few per cent of (n_eff - 1)
    m = mode_numbers(complex(-129, 3.3), 200e-9, 1.55e-6)
    delta = 1.55e-6 / (2 * math.pi * math.sqrt(129))
    assert abs((m["n_eff"] - 1) / (delta / 200e-9) - 1) < 0.1, ("delta/d check", m, delta)
    # 3b. the thin-gap formula n_eff^2 ~ eps_d (mu_d + 2 delta/d) gives, for eps*mu fixed, dn/dK ~ delta / (d n_eff)
    k0, b = 2 * math.pi / 1.55e-6, gap_mode(complex(-129, 3.3), 200e-9, 1.55e-6).real
    S_opp = (gap_mode(complex(-129, 3.3), 200e-9, 1.55e-6, 1 + 1e-6, 1 / (1 + 1e-6)).real - b) / (k0 * 1e-6)
    assert abs(S_opp / (delta / (200e-9 * m["n_eff"])) - 1) < 0.15, ("eps*mu fixed, real walls", S_opp)
    # 4. a gold-like Drude metal gives LESS pressure than ideal mirrors, and the deficit grows as the gap shrinks
    r200 = casimir_pressure(200e-9, 9.0, 0.035) / casimir_ideal(200e-9)
    r10 = casimir_pressure(10e-9, 9.0, 0.035) / casimir_ideal(10e-9)
    assert 0 < r10 < r200 < 1, ("Drude ordering", r10, r200)
    return True


# ------------------------------------------------------------------ the report
def report(I):
    out = []
    p = out.append
    lam, L, d, A, f = I["wavelength"], I["slot_length"], I["gap"], I["stroke"], I["f_drive"]
    k0 = 2 * math.pi / lam

    p("== 2. The actuator")
    v = 2 * math.pi * f * A
    p("  peak wall speed v = 2 pi f A = %.3g m/s;  v/c = %.2g" % (v, v / c))
    p("  headroom to the material limit (%.0f m/s): x%.2g" % (I["v_cap"], I["v_cap"] / v))
    p("  drive / light frequency = f lambda / c = %.2g" % (f * lam / c))
    p("  photon yield vs Wilson, if it scales as (v/c)^2: %.1g" % ((v / c / I["wilson_v_over_c"]) ** 2))
    fc = f / 2
    Tq = h * fc / kB
    n_th = 1 / math.expm1(h * fc / (kB * I["T_fridge"]))
    p("  cold-circuit option: a mode at f/2 = %.3g GHz has h f / kB = %.1f mK; thermal photons at %.0f mK: %.2g"
      % (fc / 1e9, Tq * 1e3, I["T_fridge"] * 1e3, n_th))
    p("  can the wall make that stroke? a sound wave carrying speed v strains the material by v/v_sound and stresses it")
    p("  by rho * v_sound * v (estimate, plane wave):")
    for name, (rho, vs, vt) in MATERIALS.items():
        shear = 2 * (vt / vs) ** 2            # the part of the uniaxial-strain stress that drives yield, 2 mu S / (rho c v)
        p("    %-8s strain %.2g, stress %.2g Pa (shear part %.2g Pa); wavelength: compression %.3g um, shear %.3g um;"
          % (name, v / vs, rho * vs * v, shear * rho * vs * v, vs / f * 1e6, vt / f * 1e6))
        p("             length that moves as one when a wave runs along the wall (lambda/pi): compression %.2g um, shear %.2g um"
          % (vs / f / math.pi * 1e6, vt / f / math.pi * 1e6))
    W = I["wall"]
    m_wall = W["length"] * W["height"] * W["thickness"] * W["density"]
    w = 2 * math.pi * f
    c_ext_gold = 2030.0                         # gold, thin-rod (extensional) speed, CRC Handbook
    c_bend = math.sqrt(w * W["thickness"] * c_ext_gold / math.sqrt(12))   # thin-plate bending wave (rough when kh ~ 1)
    p("    a %.0f-nm gold wall's bending wave: speed %.0f m/s, lambda/pi = %.2g um moves as one (estimate)"
      % (W["thickness"] * 1e9, c_bend, c_bend / f / math.pi * 1e6))
    p("    10^13 cycles at f take %.1f hours" % (1e13 / f / 3600))
    s_in = W["density"] * w ** 2 * A * W["thickness"]
    g_rho, g_vl, g_vt = MATERIALS["gold"]
    p("    a thin gold layer riding on the actuator's free face carries only its inertia, rho w^2 A t: %.2g Pa (shear part"
      % s_in)
    p("    %.2g Pa) for %.0f nm; thin gold films start to yield near 1e8 Pa, and GHz-cycled aluminium films have failed"
      % (s_in * 2 * (g_vt / g_vl) ** 2, W["thickness"] * 1e9))
    p("    down to 1.4e7 Pa (Eberl et al. 2006)")
    P_mech = m_wall * w ** 3 * A ** 2 / (2 * I["Q_mech"])
    p("  mechanical loss of the gold wall alone, m w^3 A^2 / 2Q: m = %.3g pg, Q = %.0f -> %.2g W"
      % (m_wall * 1e15, I["Q_mech"], P_mech))
    a_rho, a_vl, a_vt = MATERIALS["AlN"]
    t_fbar = a_vl / (2 * f)                     # a half-wave thickness-mode plate
    E_fbar = 0.25 * a_rho * w ** 2 * A ** 2 * I["actuator_area"] * t_fbar
    P_fbar = w * E_fbar / I["Q_mech"]
    p("  a thickness-mode AlN actuator at f is %.2g um thick (half a wavelength); holding the stroke on %.3g mm^2 of it"
      % (t_fbar * 1e6, I["actuator_area"] * 1e6))
    p("  stores 1/4 rho w^2 A^2 (area) t = %.2g J and dissipates w E / Q = %.2g W at Q = %.0f (estimate)"
      % (E_fbar, P_fbar, I["Q_mech"]))
    T_warm = 0.1 * math.sqrt(P_fbar / 1e-3)     # a large fridge: ~1 mW of cooling at 100 mK, cooling ~ T^2 (estimate)
    p("  dilution fridges are rated at about 1e-5 W at 20 mK (none at 10 mK); that drive would hold the fridge near %.0f mK,"
      % (T_warm * 1e3))
    p("  where the f/2 mode holds %.1f thermal photons (it holds %.2f at 20 mK) (estimate)"
      % (1 / math.expm1(h * fc / (kB * T_warm)), 1 / math.expm1(h * fc / (kB * 0.020))))

    p("== 3. The observable")
    m = mode_numbers(I["eps_metal"], d, lam)
    p("  gap mode: n_eff = %.4f, distance to 1/e power = %.3g um" % (m["n_eff"], m["L_power"] * 1e6))
    p("  phase per unit vacuum change: dphi/dK = k0 L S, with S = dn_eff/dK = %.4f" % m["S"])
    p("  = %.4g rad per unit K over L = %.3g um" % (k0 * L * m["S"], L * 1e6))
    p("  a %.0f %% change would shift the phase by %.3g rad (if the whole length L moves in step)"
      % (I["dK_example"] * 100, k0 * L * m["S"] * I["dK_example"]))
    beta0 = gap_mode(I["eps_metal"], d, lam).real
    S_opp = (gap_mode(I["eps_metal"], d, lam, 1 + 1e-6, 1 / (1 + 1e-6)).real - beta0) / (k0 * 1e-6)
    p("  if eps and mu change oppositely (eps*mu fixed, only the impedance moves): S = %.3f, %.0f %% of the above;"
      % (S_opp, 100 * S_opp / m["S"]))
    p("  between perfect conductors it would be 0, and then the reflection carries the signal")
    dn_perp = 11 * math.pi ** 2 / 8100 * alpha ** 2 * (lambda_C_reduced / d) ** 4
    p("  standard QED (Scharnhorst): 0 along the plates; %.2g across a gap of %.0f nm" % (dn_perp, d * 1e9))

    p("== 4. The mimics (sizes as an equivalent vacuum change)")
    breath = m["dn_dd"] / m["S"]
    p("  breathing: dn_eff/dd = %.3g per nm  ->  %.2g of K per nm of gap change" % (m["dn_dd"] * 1e-9, abs(breath) * 1e-9))
    p("  at this stroke (amplitude %.3g nm): %.2g of K at the drive frequency" % (A * 1e9, abs(breath) * A))
    p("  gas in the gap (air, n - 1 = 2.7e-4): breathing changes it by %.2g" % (2.7e-4 * A / d))
    p("  round trip of the slot at n_eff: %.2g ps (the group index is larger)" % (2 * L * m["n_eff"] / c * 1e12))

    p("== 5. The ledger")
    F_max = I["P_in"] / c
    lo, hi = I["claimed_thrust"]
    p("  F <= P/c at %.0f W: %.3g N;  the claimed %.3g-%.3g N is x%.2g-x%.2g of it" % (I["P_in"], F_max, lo, hi, lo / F_max, hi / F_max))
    for dd in (d, 10e-9):
        Pid = casimir_ideal(dd)
        Pau = casimir_pressure(dd, 9.0, 0.035)
        p("  Casimir pressure at %.0f nm: ideal mirrors %.3g Pa (%.3g N/cm^2); gold-like Drude metal %.3g Pa (x%.2f of ideal)"
          % (dd * 1e9, Pid, Pid / 1e4, Pau, Pau / Pid))
    p("  (the Drude model keeps only gold's free electrons; with its interband transitions real gold gives about 0.6 of")
    p("   ideal at 200 nm, measured, and about 0.12 at 10 nm: Decca et al. 2007; Lambrecht & Reynaud 2000)")
    p("  pre-seed arithmetic: 0.6-9 N/cm^2 over 100 cm^2 = %.0f-%.0f N" % (0.6 * 100, 9 * 100))

    p("== 6. The bench floor")
    N = I["P_detected"] * I["t_integrate"] * lam / (h * c)
    dphi = 1 / math.sqrt(N)
    p("  photon noise: dphi = 1/sqrt(N) = %.2g rad  ->  dK = %.2g" % (dphi, dphi / (k0 * L * m["S"])))
    slot_dB = 10 * math.log10(math.e) * L / m["L_power"]
    lo_dB = slot_dB + 2 * I["coupler_dB"][0] + 2 * I["fiber_dB"][0]
    hi_dB = slot_dB + 2 * I["coupler_dB"][1] + 2 * I["fiber_dB"][1]
    p("  light budget: the slot itself costs %.1f dB; with the couplers, %.0f-%.0f dB in all, so %.2g W detected needs"
      % (slot_dB, lo_dB, hi_dB, I["P_detected"]))
    p("  %.2g-%.2g W launched (estimate)" % (I["P_detected"] * 10 ** (lo_dB / 10), I["P_detected"] * 10 ** (hi_dB / 10)))
    floor = abs(breath) * A * I["stroke_cal"]
    p("  breathing known to %.2g %%: floor dK = %.2g" % (I["stroke_cal"] * 100, floor))
    for target in (1e-7, 1e-6):
        p("  to reach dK = %.0e the stroke must be known to %.2g pm (1 part in %.0f)"
          % (target, target / abs(breath) * 1e12, abs(breath) * A / target))
    return "\n".join(out)


def table_b():
    rows = [("Gold, 1.5 um", complex(-106.9, 10.2), 1.5e-6), ("Silver, 1.55 um", complex(-129, 3.3), 1.55e-6),
            ("Aluminium, 1.5 um", complex(-226, 46), 1.5e-6), ("Gold, 350 nm", complex(-1.26, 5.58), 350e-9)]
    out = ["== Table B (gap 200 nm; permittivities: Johnson & Christy for gold and silver, Rakic for aluminium)"]
    for name, eps, lam in rows:
        m = mode_numbers(eps, 200e-9, lam)
        out.append("  %-18s eps = %-14s n_eff = %.3f   1/e power at %.3g um   S = %.3f"
                   % (name, "%g%+gi" % (eps.real, eps.imag), m["n_eff"], m["L_power"] * 1e6, m["S"]))
    return "\n".join(out)


if __name__ == "__main__":
    self_tests()
    print("self-tests passed\n")
    print(report(INPUTS))
    print()
    print(table_b())
