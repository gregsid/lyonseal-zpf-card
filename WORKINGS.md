# Your six steps, filled in for the ZPF Array

## What I assumed

| Input | Value | Where it comes from |
|---|---|---|
| Wall drive | about 1 GHz | your Apr 19 post: "GHz frequency dynamic boundary motion" |
| Wall stroke | 1 nm amplitude (2 nm peak to peak) | assumed |
| Gap width | 200 nm | assumed: the "illustrative gap" in your Oct 25 post (your site doesn't give the Array's) |
| Probe light | 1.55 µm, guided along the gap | my choice |
| Slot length | 100 µm, silver walls | my choice (gold loses most of the light within 23 µm) |
| Input power | 10 W | your pre-seed page, as read Nov 1, 2026 |

Numbers marked "estimate" are worked from these with standard physics and published data. The script, `zpf_card_workings.py`, redoes the page's calculations from these inputs (the actuator and its heat, the observable, breathing, the ledger, the light budget, the bench floor and Table B), checks itself against cases with known answers, and updates when you change the stroke, gap, drive frequency, probe wavelength, slot length, wall permittivity and the powers. Table B's four rows, the Casimir metal model and the wall materials' sound speeds are fixed in it. Figures taken from papers are in the source list at the end, under the step they support.

## 1. The degree of freedom

Your report defines K for the whole cosmos (§3.2) and says in §3.4 that the cosmic K can't automatically be applied to every local change. Your About page says a local ZPE rise changes permittivity "in balance with" the magnetic properties. Your Jan 20 post defines a local Km raised by resonance. These are different claims, and they predict different things. Your Oct 25 post says such questions "should drive the experiment before the first measurement is taken". The form below asks you to pick one and give it a number.

**Forecast form — fill in and date before the first run**

| Question | Your answer |
|---|---|
| What changes | e.g. ε and μ together, by K; or oppositely, keeping εμ fixed |
| By how much | ΔK = … per nm of stroke, at this gap and frequency |
| Where | the gap only, or the gap and the walls |
| When | steady, at the drive frequency Ω, or at 2Ω |
| Phase against the wall motion | with the displacement, or 90° from it (with the velocity) |
| Which light sees it | along the gap or across it; which polarization; which wavelengths |
| What it needs | resonance? one cavity, or the whole array? |
| How it scales | with stroke, drive frequency and gap width |
| A pass | a phase change of at least … |
| A fail | a null below … |
| What would convince a skeptic | … |

The first line decides what the bench should read, and reading both the phase and the reflection tells the two cases apart. If ε and μ change together, the phase moves and the reflection stays small (see the reflection note under Table B). If they change oppositely, keeping εμ fixed, the phase moves only about a tenth as much with real metal walls (not at all between perfect conductors), and a sharp edge of the changed region reflects about ten times more (estimates). Measuring reflection and transmission together is the standard way to separate ε from μ in a material (the Nicolson–Ross–Weir method).

Post the filled form with its date, or send it to someone outside the company, before the first run, so no one can say the result was fitted afterwards.

## 2. The actuator

- Walls moving 1 nm at 1 GHz reach about 6 m/s, which is 2×10⁻⁸ of light speed. The material itself limits any wall to about 50 m/s, whatever the frequency (Dodonov 2010, §3.1), so driving harder gains at most about 8×.
- Even 6 m/s is hard on the materials (estimates). Where the motion travels as a full sound wave, it stresses the wall by about ρ × (speed of sound) × v: about 0.4 GPa in gold, 0.2 GPa in AlN and 0.1 GPa in silicon (silver, which the optics prefers, isn't worked out here). In gold the shear part, which is what makes a metal yield, is about 0.1 GPa, and thin gold films start to yield near 0.1 GPa. A thin gold layer riding on the actuator's free face carries only its own inertia, about a fifth of that for 100 nm, but that is still above the 14 MPa at which aluminium films cycled at GHz have worn out, with no lower limit found, over more than 10¹⁴ cycles (Eberl et al. 2006); a few hours at 1 GHz is 10¹³ cycles. A reported GHz piezo transducer moves about 0.25 nm (at 1.06 GHz, 5 V). So 1 nm at 1 GHz needs a high-Q mechanical resonance driven by a low-loss piezoelectric such as AlN, and is optimistic for hours of running. PZT loses too much at 1 GHz, and off resonance 1 nm would take about 250 V.
- The drive is about 5 millionths of the frequency of 1.55-µm light. At that pace the light simply follows the walls, and standard physics predicts no new photons in it (Dodonov 2010; Wilson et al. 2011).
- Your Apr 19 post says η_b ≈ 0.08 was "derived from Wilson et al. dynamical Casimir photon production". Wilson's boundary was electrical, moved at a few percent of light speed ("approaching v/c ∼ 0.05", the paper's body), and was pumped at twice the frequency of the photons it made. If photon yield scales as (v/c)², as it does for moving mirrors (Lambrecht et al. 1996), your walls would give about 2×10⁻¹³ of Wilson's yield (estimate).
- Wilson's paper, citing Braggio et al., puts a mechanical mirror moving 1 nm at 2 GHz at about one photon a day. Those would be microwave photons from the wall motion itself, not light in the gap.
- If you want a dynamical-Casimir signal anyway, the place it could appear is an electrical circuit that uses the gap as its capacitor and resonates near half your drive frequency. It would test standard physics, not K. It needs about 10 mK: a 0.5-GHz mode holds about 0.1 thermal photons there, and 0.4 at 20 mK. Commercial dilution refrigerators are rated to remove about 10–30 µW at 20 mK; none is rated at 10 mK, where they manage a few µW at most. A thickness-mode drive holding 1 nm dissipates about 10 mW for each 100 µm × 100 µm of actuator at a quality factor of 1,000 (estimate), so a continuous drive would warm the refrigerator to a few hundred mK, where the mode holds a dozen or so thermal photons (estimate).

## 3. The observable

- If K changes by ΔK inside the gap, light guided along it shifts phase by about k₀ × L × ΔK. Here k₀ = 2π divided by the wavelength, and L is the length. Real metal walls change this by under 3% (estimate; 0.7% for silver at 1.55 µm).
- Standard quantum electrodynamics predicts no change for light along the plates (Scharnhorst 1990; Barton 1990), and about 10⁻²⁹ for light across a 200-nm gap. So the size has to come from your framework. The slot itself still answers the drive (§4, the mimics: about 0.2 radians per nanometre of gap change at this length), so what the bench looks for is an excess beyond those modelled responses, and the 10⁻⁷ in §6 is the calibration that excess would need, not a demonstrated floor.
- For scale: your Aug 4, 2025 post "Emergent Thermodynamics in Stochastic Electrodynamics: Implications of Variable Zero-Point Energy Density…" (zero-pointfieldtechnologies.com/post/emergent-thermodynamics-in-stochastic-electrodynamics-implications-of-variable-zero-point-energy-de) says c "varies ~5-10%" with the ZPE. If the Array did even 5% locally, light along a 100-µm silver slot at 1.55 µm would shift by about 20 radians — if the whole length moves in step (next point).
- At 1 GHz, sound in the walls has a wavelength of only a few µm: compression waves are about 3 µm long in gold, 8 µm in silicon and 11 µm in AlN. The waves that move a wall sideways as they run along it are shorter still: shear waves about 1.2 µm and bending waves about 0.6 µm in a thin gold wall. So a wall pushed from one point carries a wave, and at most a micron or two of a 100-µm slot moves in step. That shortens the effective L for any effect that follows the drive; a steady effect doesn't pay that penalty.
- To move the whole length in step, drive it evenly along its length. A thickness-mode (FBAR-type) actuator under the whole wall does this; at 1 GHz it is about 5 µm of AlN, and with a tuned edge frame its face moves as a nearly flat piston (Pensala & Ylilammi 2009). Its face moves along the plate's normal, so the moving wall must be the broad face of a plate many times wider than its 5-µm thickness: the slot becomes a stack of plates rather than two upright walls (estimate).

## 4. The mimics

Your Sep 8 post says a departure "must survive electromagnetic, thermal, mechanical, material and instrumental controls". Table A is that list for an optical readout of the gap, with the drive at frequency Ω. Sizes are estimates, given as the vacuum change that would give the same phase; where the size depends on your design, the column says so.

**Table A — the mimics**

| Mimic | Size, as an equivalent ΔK | Shows up at | Control |
|---|---|---|---|
| The drive changes the gap ("breathing") | 5×10⁻⁴ per nm of gap change | Ω, with the displacement | Scan the stroke and the relative phase of the two walls; scan gap width (breathing goes as 1/d²) |
| Heating from strain | depends on the design (about 1 K per 10⁻³ of volume strain) | Ω | Scan the drive frequency; model the heat |
| Ohmic and dielectric heating in the actuator | depends on the design | 0 and 2Ω (Ω only with a DC bias) | Scan the bias; switch the drive on and off slowly |
| Electrostatic force between the walls: drive voltage, DC offsets, patch potentials | depends on the design | 2Ω, and Ω with any DC offset | Put the drive on the actuator's buried electrode and hold the gap-facing face at RF ground; measure the contact potential between the walls and cancel it with a DC bias |
| Gas in the gap | about 10⁻⁶ (air's index is 3×10⁻⁴, and breathing changes it) | Ω | Run in vacuum |
| Water adsorbed on the walls | shifts the mode | steady, drifting | Bake out; run in vacuum |
| Heating by the probe light | depends on the design | steady | Scan the probe power |
| Electro-optic effects in nearby optics | set by the drive's stray field | Ω | Shield; move the optics |
| Electrical pickup in the detector | can fake any sideband | Ω | Light-off runs; heterodyne with an offset, so the optical signal appears at Ω ± δ while pickup stays at Ω |
| Laser frequency noise | grows with path mismatch | broadband | Balance the interferometer's paths |
| Motion of the output optics | small, with tilt | Ω | Rigid mounts; alignment scans |
| Light reflected back by the slot's exit | the reflected sideband can be 30–60% of the transmitted one | Ω | Matched or tapered exits (a round trip takes about 0.7 ps, too short to gate in time) |
| Sound reaching the couplers | depends on the design | Ω | Isolate the couplers; scan the drive frequency |
| Strain or heating in the walls' metal | depends on the design | Ω | Scan the gap width (this goes about as 1/d) and change the wall metal or its thickness; an undriven section of the slot tracks drift but can't separate this from K |

If your effect is steady or at 2Ω, heating and electrostatic forces become direct mimics. Moving both walls together removes breathing but may remove your effect too; scanning the relative phase of the walls separates the two. A gap-width scan separates three things at once: breathing goes about as 1/d², a change in the walls' metal about as 1/d, and a change of K in the gap hardly depends on d (estimate, from the script's solver). On a high-Q resonance the stroke changes with the drive frequency, so normalize every frequency scan to the measured stroke.

## 5. The ledger

- Your home page gives F ≤ P_rad/c for radiated light. It says anything above that "would require evidence for an additional, explicitly identified momentum-transport mechanism". At 10 W the limit is 33 nN. Name the mechanism, and where the reaction goes, before the thrust stand. The energy half too: where the 10 W ends up.
- An isolated system can't move its own centre of mass without coupling to something outside it (Puthoff, Little & Ibison 2002, §2.1). A static force such as the Casimir pull acts on both walls and gives no net thrust.
- For scale: between ideal mirrors, the whole static Casimir pressure is about 0.8 Pa across a 200-nm gap and about 13 N/cm² across a 10-nm gap. Real gold or silver gives about 60% of that at 200 nm (about 0.5 Pa, measured for gold) and only about 12% at 10 nm (about 1.5 N/cm², estimate), where the force is closer to the van der Waals 1/d³ law than the Casimir 1/d⁴.
- At the thrust stand, the standard checks are heat, forces from cables and magnetic fields, vibration and outgassing. To test the P/c limit, the stand needs a resolution of about 10 nN.

## 6. The experiment

1. **Drive off.** Measure the slot's effective index and loss, and compare with Table B. A mismatch means the slot isn't what you think (width, roughness, film quality).
2. **Drive on.** Read the transmitted phase at Ω with every control above, in both quadratures: breathing follows the walls' displacement, and anything that follows their velocity sits 90° away. Measure the stroke independently, for example with a second interferometer on the back of a wall: calibrating it by the breathing signal alone is circular if your effect also grows with stroke. For a steady effect, switch the drive on and off slowly and lock in to that, to beat thermal drift.
3. **Thrust stand last,** with the ledger written down.

How small a ΔK the bench can see (estimates):
- Photon noise allows about 10⁻¹² (0.1 mW reaching the detector through a 100-µm silver slot, a few hours). Getting 0.1 mW out takes roughly 1–10 mW launched (estimate): in measured devices, getting light into a metal slot and out again costs about 1.4–5 dB each way, plus 1.4–2.8 dB at each fiber coupler, and 100 µm of silver costs about 5 dB more. Most of the light lost in the slot heats the walls, so how much they can absorb is the open question (Table A, heating by the probe light).
- In practice the limit is how well breathing is calibrated out. Reaching 10⁻⁷ needs the stroke known to about 0.2 pm on 1 nm (1 part in 5,000); at 1% it is about 5×10⁻⁶.
- A steady drive-on/drive-off shift is limited by thermal drift instead, likely far worse.

Reading the outcome:
- A signal above your stated pass, surviving the scans, is worth replicating.
- A null below your stated fail is a publishable bound.
- Anything in between means improving the calibration first.

## Table B — the light a 200-nm gap holds

Estimates for an air gap between infinite plates, with bulk metal data; check with a mode solver for your geometry (the script solves the same equation for your metal data).

| Walls, wavelength | Effective index | Distance to 1/e power |
|---|---|---|
| Gold, 1.5 µm | 1.11 | 23 µm |
| Silver, 1.55 µm | 1.10 | 95 µm |
| Aluminium, 1.5 µm | 1.08 | 16 µm |
| Gold, 350 nm | 1.08 | 0.26 µm |

Only light polarized across the gap is guided. The next mode up exists only below about 450–500 nm, where the metals absorb strongly.

**Reflection.** In an ideal gap between perfect conductors, a change of ε and μ together reflects nothing. With real gold or silver walls, part of the light runs in the metal, so the change reflects about 6–15% as much as an ordinary index change giving the same phase. This is a different quantity from the 30–60% exit reflection in Table A, which is the slot's own end. A change that fades in smoothly along the slot reflects almost nothing either way, so a quiet reflection is weak evidence on its own.

## Sources

- Your report: Douglas Miller, ZPF Technologies LLC, *The Setterfield–Miller Zero-Point Field Framework*, Public Technical Report v1.1, October 2026 (§3.2, §3.4; step 1).
- Your pages, opened Oct 31, 2026: the home page, About and pre-seed (zero-pointfieldtechnologies.com); the posts "Hal Puthoff: The Consummate Salesman of the Physical Vacuum" (Sep 7, 2026, the six steps), "Whatever Happened to the Physical Vacuum?" (Sep 8, 2026), "The Tiny Cavities That Tried to Tap the Quantum Vacuum" (Oct 25, 2026), "Quantitative Vacuum Engineering: Predictions for Density Surge and Thrust…" (Apr 19, 2026), the Jan 20, 2026 post on the Miller Constant, and "Emergent Thermodynamics in Stochastic Electrodynamics: Implications of Variable Zero-Point Energy Density for Physics and Technology" (Aug 4, 2025; one of two posts that day opening "Emergent Thermodynamics" — the address is given at step 3).
- Wilson et al. 2011, Nature 479, 376 (step 2): https://arxiv.org/abs/1105.4714
- Dodonov 2010, "Current status of the dynamical Casimir effect", §3.1 (step 2): https://arxiv.org/abs/1004.3301
- Lambrecht, Jaekel & Reynaud 1996, PRL 77, 615, motion-induced radiation (step 2): https://arxiv.org/abs/quant-ph/9606029
- Dilution-refrigerator cooling power: the Bluefors LD and XLD specifications (bluefors.com) and the Quantum Design Oxford ProteoxLX specification (step 2).
- Scharnhorst 1990, Phys. Lett. B 236, 354, and Barton 1990, Phys. Lett. B 237, 559 (step 3); the formulas as reviewed in Scharnhorst 1998, "The velocities of light in modified QED vacua": https://arxiv.org/abs/hep-th/9810221
- Sound speeds: the CRC Handbook (gold) and the Ioffe database (silicon, AlN). Piston-mode thin-film resonators: Pensala & Ylilammi, IEEE TUFFC 56, 1731 (2009) (step 3).
- Separating ε from μ by reflection and transmission (step 1): Nicolson & Ross, IEEE Trans. Instrum. Meas. IM-19, 377 (1970); Weir, Proc. IEEE 62, 33 (1974); Smith et al., PRB 65, 195104 (2002), https://arxiv.org/abs/physics/0111203
- Mirmoosa, Setälä & Norrman 2023, §3.8 (https://arxiv.org/abs/2312.15178), and Galiffi et al. 2022, §3 (https://arxiv.org/abs/2111.08640): in a uniform medium, a change of ε and μ together creates no new photons and sends nothing back (steps 2 and 4).
- Davis, Teofilo, Haisch, Puthoff, Nickisch, Rueda & Cole, "Review of experimental concepts for studying the quantum vacuum field", STAIF 2006, AIP Conf. Proc. 813: https://www.bu.edu/simulation/publications/dcole/PDF/Davis%20et%20al_STAIF06_Log063.pdf. It proposes a laser between Casimir plates compared with beams outside them, and RF-driven piezoelectric cavity mounts; worth reading as prior art for your design (steps 3 and 6).
- Puthoff, Little & Ibison, JBIS 55, 137 (2002), §2.1 (step 5); arXiv copy: https://arxiv.org/abs/1012.5264
- Casimir pressure with real metals (step 5): Decca et al. 2007, Eur. Phys. J. C 51, 963 (https://arxiv.org/abs/0706.3283), measured for gold, 510 mPa at 200 nm; Lambrecht & Reynaud 2000, Eur. Phys. J. D 8, 309 (https://arxiv.org/abs/quant-ph/9907105); Palasantzas et al. 2008, Appl. Phys. Lett. 93, 121912 (https://arxiv.org/abs/0809.1779), the change to the van der Waals regime below about 18 nm.
- Strength and fatigue at GHz (step 2): Eberl et al. 2006, J. Appl. Phys. 99, 113501 (aluminium films fatigued at GHz); Gudlavalleti, Kumar & Anand 2001, MRS Proc. 695 (yield of thin gold films); the 0.25-nm GHz transducer is in US patent 12,042,794 B2.
- Getting light into and out of metal slots (step 6): Kohli et al. 2025, Light Sci. Appl. 14, 399 (https://www.nature.com/articles/s41377-025-02116-1); Chen et al. 2015, J. Lightwave Technol. 33, 535 (https://arxiv.org/abs/1503.03805).
- PVLAS (Ejlli et al. 2020), https://arxiv.org/abs/2005.12913. Context only: the most sensitive vacuum-optics measurement, birefringence in a 2.5 T field, reached about 10⁻²²; it measures a polarization difference, not this kind of index change.
- Optical constants: Johnson & Christy 1972 (gold, silver) and Rakić 1995 (aluminium), via refractiveindex.info (Table B and the script).
