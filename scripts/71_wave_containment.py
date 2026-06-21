#!/usr/bin/env python3
"""
RUNG-37 — WAVE CONTAINMENT: the coupled R35×R36 safety model. Does the tumour-clearing wave stay contained in a
sea of normal tissue, or does it run away once a normal cell falsely fires?

THE HOLE THIS CLOSES (what R35 and R36 each missed)
---------------------------------------------------
R36 proved the bystander wave must be SUPER-critical (b>p_c) to clear the tumour from a partial seed -- but it ran
the wave on a PURE-TUMOUR lattice. R35 proved normal cells falsely fire at a small rate (the leak). Put them
together and a danger neither rung asked appears: the wave is RESISTANCE-AGNOSTIC (it kills a neighbour regardless
of whether it is cancer). The tumour is a small island in a SEA of normal cells. So:
  * tumour clearance NEEDS the wave super-critical (R36),
  * but a super-critical wave in NORMAL tissue, ignited by even one R35 false-fire (or by spillover across the
    tumour boundary), would PERCOLATE through normal tissue -> catastrophic.
=> The wave must be BISTABLE ACROSS TISSUE TYPE: super-critical inside the tumour, SUB-critical in normal tissue.
What creates that difference is a RECOGNITION GATE on the bystander signal: a cell only completes death from the
wave if it ALSO recognises itself as cancer (its own mutation-sensor primed) -- so the death signal that spreads
freely among tumour cells is DAMPED (low effective bond prob) when it reaches a normal cell. (RUNG-12P/B's
"per-hop-gated relay" + RUNG-13/14's "contained", made into the explicit normal-tissue safety test.)

THE MODEL (two-tissue lattice CA; extends scripts/70)
-----------------------------------------------------
LxL lattice. A central tumour DISK (radius R) = TUMOUR cells; everything else = NORMAL cells.
  - TUMOUR ignition: delivered w.p. f, fires w.p. p_kill_t (R35 mutant kill ~0.99) -> DYING.
  - NORMAL ignition: each normal cell falsely fires w.p. p_ig_n  = the R35 normal-cell false-death (the leak).
  - WAVE (type-gated): a DYING cell recruits each ALIVE neighbour w.p. b that depends on the NEIGHBOUR's type:
        b_t  if the neighbour is TUMOUR   (free spread = super-critical, clears the tumour)
        b_n  if the neighbour is NORMAL   (recognition-GATED = damped; the safety knob)
    A cell with k dying neighbours is recruited w.p. 1-(1-b_cell)^k. (b_n << b_t is the recognition gate.)
Measure: tumour_cleared (= dead tumour / tumour) AND normal_killed (= dead normal / normal), the collateral.

THE QUESTION: for what b_n does the wave clear the tumour (b_t super-critical) while keeping normal-tissue death
contained? And does the wave AMPLIFY the R35 leak (a tiny p_ig_n exploding) unless b_n is sub-critical?

HONEST CEILING
--------------
2D lattice CA, sharp tumour/normal boundary (real margins are graded/infiltrative); f, p_kill, b_t, b_n, p_ig_n
are EFFECTIVE probabilities, not measured. The recognition gate is modelled as a lower bond prob into normal
cells, not a molecular mechanism. Robust claim = the percolation REQUIREMENT (b_n<p_c) and the leak-AMPLIFICATION,
not absolute fractions.

USAGE
  python scripts/71_wave_containment.py selftest
  python scripts/71_wave_containment.py run        # -> runs/rung37_wave_containment/ (CPU, ~1-2 min)
"""
import os, sys, json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "runs/rung37_wave_containment")
PC_2D = 0.5  # square-lattice bond-percolation threshold


def make_tumour_mask(L, frac=0.126):
    """Central disk of area `frac` of the lattice = tumour; rest = normal."""
    R = np.sqrt(frac * L * L / np.pi)
    yy, xx = np.mgrid[0:L, 0:L]
    c = (L - 1) / 2.0
    return (xx - c) ** 2 + (yy - c) ** 2 <= R * R


def simulate(L, is_tumour, f, p_kill_t, b_t, b_n, p_ig_n, rng):
    """One realisation. Returns (tumour_cleared, normal_killed)."""
    ALIVE, DYING, DEAD = 0, 1, 2
    state = np.zeros((L, L), dtype=np.int8)
    is_normal = ~is_tumour
    # type-dependent bond prob of being RECRUITED (set by the receiving cell's type)
    b_map = np.where(is_tumour, b_t, b_n).astype(float)
    # ignition
    ig_t = is_tumour & (rng.random((L, L)) < f) & (rng.random((L, L)) < p_kill_t)
    ig_n = is_normal & (rng.random((L, L)) < p_ig_n)
    state[ig_t | ig_n] = DYING
    while np.any(state == DYING):
        dying = (state == DYING)
        k = np.zeros((L, L), dtype=np.int16)
        for ax in (0, 1):
            for s in (1, -1):
                sh = np.roll(dying, s, axis=ax)
                idx = [slice(None), slice(None)]
                idx[ax] = (0 if s == 1 else -1)
                sh[tuple(idx)] = False  # hard edges, no torus
                k += sh.astype(np.int16)
        alive = (state == ALIVE)
        recruit_p = 1.0 - (1.0 - b_map) ** k
        recruited = alive & (rng.random((L, L)) < recruit_p)
        state[dying] = DEAD
        state[recruited] = DYING
    dead = (state == DEAD)
    tum = is_tumour.sum()
    nrm = is_normal.sum()
    return float(dead[is_tumour].sum() / tum), float(dead[is_normal].sum() / nrm)


def mean_run(L, is_tumour, f, p_kill_t, b_t, b_n, p_ig_n, reps=6, seed0=0):
    tc, nk = [], []
    for i in range(reps):
        a, b = simulate(L, is_tumour, f, p_kill_t, b_t, b_n, p_ig_n, np.random.default_rng(seed0 + i))
        tc.append(a); nk.append(b)
    return float(np.mean(tc)), float(np.mean(nk))


def selftest():
    L = 120
    tum = make_tumour_mask(L)
    # b_n = 0: normal cells die only by direct false-ignition -> normal_killed ~ p_ig_n
    _, nk = mean_run(L, tum, 1.0, 0.99, 0.8, 0.0, 0.01, reps=8)
    assert abs(nk - 0.01) < 0.01, nk
    # ungated (b_n = b_t super-critical) + any normal ignition -> runaway normal death
    _, nk_bad = mean_run(L, tum, 1.0, 0.99, 0.8, 0.8, 0.001, reps=6)
    assert nk_bad > 0.5, nk_bad
    # tumour cleared by a super-critical b_t even at low delivery
    tc, _ = mean_run(L, tum, 0.1, 0.99, 0.8, 0.0, 0.0, reps=6)
    assert tc > 0.9, tc
    # sub-critical b_n contains; super-critical b_n explodes (the percolation requirement)
    _, nk_sub = mean_run(L, tum, 1.0, 0.99, 0.8, 0.35, 0.001, reps=6)
    _, nk_sup = mean_run(L, tum, 1.0, 0.99, 0.8, 0.65, 0.001, reps=6)
    assert nk_sup > nk_sub + 0.2, (nk_sub, nk_sup)
    print(f"[selftest] containment limits OK (b_n0->nk={nk:.3f}, ungated->{nk_bad:.3f}, "
          f"tum_clear={tc:.3f}, sub={nk_sub:.3f}<sup={nk_sup:.3f})")


def main():
    os.makedirs(OUT, exist_ok=True)
    selftest()
    res = {"tag": "rung37_wave_containment", "pc_2d": PC_2D,
           "model": "two-tissue lattice CA; wave bond prob set by receiving cell type (recognition gate)"}
    L = 200
    tum = make_tumour_mask(L, frac=0.126)
    print(f"\n=== RUNG-37: wave containment (L={L}, tumour disk = 12.6% centre, b_t=0.80 super-critical) ===")

    # --- 1. the danger: an UNGATED super-critical wave is a runaway in normal tissue ----------------------------
    print("\n-- 1. UNGATED wave (b_n = b_t = 0.80): does one normal false-fire blow up normal tissue? --")
    print("   p_ig_n(R35 leak)   tumour_cleared   normal_KILLED")
    ungated = []
    for p_ig in (0.0, 1e-5, 1e-3, 1e-2):
        tc, nk = mean_run(L, tum, 0.1, 0.99, 0.8, 0.8, p_ig, reps=6)
        ungated.append({"p_ig_n": p_ig, "tumour_cleared": round(tc, 3), "normal_killed": round(nk, 3)})
        print(f"      {p_ig:.0e}            {tc:.3f}           {nk:.3f}")
    res["ungated_wave"] = ungated
    print("   -> a super-critical wave with NO recognition gate kills ~all normal tissue once ANY cell ignites.")

    # --- 2. the fix: a RECOGNITION-GATED wave (sweep b_n) keeps the kill in the tumour --------------------------
    print(f"\n-- 2. RECOGNITION-GATED wave: sweep b_n (normal-tissue bond prob), b_t=0.80, p_ig_n=1e-3 --")
    print("   b_n     tumour_cleared   normal_killed   regime (normal tissue)")
    gated = []
    for b_n in (0.0, 0.1, 0.2, 0.35, 0.45, 0.5, 0.6, 0.8):
        tc, nk = mean_run(L, tum, 0.1, 0.99, 0.8, b_n, 1e-3, reps=6)
        reg = "CONTAINED (sub-critical)" if b_n < PC_2D else ("AT p_c" if b_n == PC_2D else "RUNAWAY (super-critical)")
        gated.append({"b_n": b_n, "tumour_cleared": round(tc, 3), "normal_killed": round(nk, 4), "regime": reg})
        print(f"  {b_n:.2f}      {tc:.3f}          {nk:.4f}        {reg}")
    res["gated_wave_sweep"] = gated
    print("   -> tumour clearance is ~flat (b_t carries it); normal_killed stays tiny ONLY while b_n < p_c=0.5,")
    print("      then explodes. SAFETY REQUIREMENT: the recognition gate must keep b_n SUB-critical.")

    # --- 3. the wave AMPLIFIES the R35 leak unless contained (couple the two rungs numerically) -----------------
    print(f"\n-- 3. does the wave amplify the R35 leak? normal_killed vs the R35 false-death p_ig_n --")
    print("   (b_t=0.80; compare a CONTAINED gate b_n=0.35 vs an UNGATED wave b_n=0.80)")
    print("   p_ig_n (R35)    normal_killed[b_n=0.35]   amplification    normal_killed[b_n=0.80]   amplification")
    amp = []
    for p_ig in (1e-5, 1e-4, 1e-3, 1e-2):
        _, nk_c = mean_run(L, tum, 0.0, 0.0, 0.8, 0.35, p_ig, reps=8)   # no tumour seed: isolate normal-tissue fate
        _, nk_u = mean_run(L, tum, 0.0, 0.0, 0.8, 0.80, p_ig, reps=8)
        amp_c = nk_c / p_ig if p_ig > 0 else float("nan")
        amp_u = nk_u / p_ig if p_ig > 0 else float("nan")
        amp.append({"p_ig_n": p_ig, "nk_contained_bn0.35": nk_c, "amp_contained": round(amp_c, 1),
                    "nk_ungated_bn0.80": nk_u, "amp_ungated": round(amp_u, 1)})
        print(f"     {p_ig:.0e}          {nk_c:.2e}             {amp_c:>6.1f}x          {nk_u:.2e}            {amp_u:>7.1f}x")
    res["leak_amplification"] = amp
    print("   -> the wave amplifies the leak by the MEAN normal-tissue death-cluster size: CONTAINED (b_n=0.35,")
    print("      sub-critical) bounds it to ~tens-fold (->1x only as b_n->0); UNGATED (b_n=0.80) DIVERGES -> a")
    print("      1e-5 leak becomes total death (>1e4 x). So R35's leak is only safe behind a gate WELL below p_c,")
    print("      and the effective normal false-death is ~ (R35 leak) x (normal-tissue cluster size).")

    # --- 4. the unavoidable collateral RIM (even a contained wave kills a margin of normal tissue) --------------
    print(f"\n-- 4. collateral RIM: normal cells killed at the tumour boundary even when CONTAINED (p_ig_n=0) --")
    print("   b_n     tumour_cleared   normal_killed (boundary rim only)")
    rim = []
    for b_n in (0.0, 0.1, 0.2, 0.35, 0.45):
        tc, nk = mean_run(L, tum, 0.5, 0.99, 0.8, b_n, 0.0, reps=8)
        rim.append({"b_n": b_n, "tumour_cleared": round(tc, 3), "normal_rim_killed": round(nk, 4)})
        print(f"  {b_n:.2f}      {tc:.3f}          {nk:.4f}")
    res["collateral_rim"] = rim
    print("   -> even sub-critical, the tumour wave kills a thin RIM of normal cells at the boundary (grows with")
    print("      b_n) = an intrinsic 'surgical margin'. Trade-off: higher b_n clears the tumour edge better but")
    print("      widens the normal rim. The gate sets where on that trade-off you sit.")

    verdict = (
        "The bystander wave must be BISTABLE ACROSS TISSUE TYPE -- super-critical inside the tumour (to clear it, "
        "R36) AND sub-critical in normal tissue (or it runs away). (1) An UNGATED super-critical wave is a "
        "catastrophe by BOUNDARY SPILLOVER alone: the tumour-clearing wave crosses the tumour edge into normal "
        "tissue (b_n super-critical) and percolates through ~all of it -- even with ZERO normal leak (p_ig_n=0 -> "
        "0.998 normal killed); a normal false-fire merely adds ignition sources. The wave does not know the tumour "
        "boundary. (2) A RECOGNITION-GATED wave (the bystander death completes only in cells whose own mutation-"
        "sensor is primed -> a much lower effective bond prob b_n into normal cells) contains the kill: tumour "
        "clearance is unchanged (b_t carries it) while normal-tissue death stays small -- but ONLY while b_n < "
        "p_c=0.5, and really b_n must sit WELL below p_c (normal_killed jumps 0.025->0.16->0.61 across b_n="
        "0.35->0.45->0.50). (3) The wave AMPLIFIES the R35 leak by the MEAN normal-tissue death-cluster size: a "
        "sub-critical gate bounds the amplification (~tens-fold at b_n=0.35, ->1x only as b_n->0) but it DIVERGES "
        "as b_n->p_c, and an ungated wave turns a 1e-5 leak into total death (>1e4 x). So the effective normal "
        "false-death is ~ (R35 leak) x (normal cluster size) -- R35's leak target must be TIGHTENED by that factor, "
        "and is only safe behind a gate well below p_c. (4) Even contained, the wave kills a thin RIM of normal "
        "tissue at the tumour boundary (an intrinsic 'surgical margin' that grows with b_n: ~1.2% at b_n=0.35, "
        "~6% at 0.45). NET: the design must engineer a wave that is super-critical among recognised (tumour) cells "
        "and sub-critical among unrecognised (normal) cells -- a RECOGNITION-GATED bystander signal with b_n well "
        "below p_c. This is a new, named safety requirement that R36's pure-tumour percolation hid; it sits beside "
        "R35's leak-correlation as the second load-bearing residual of the kill-coupling.")
    print("\nVERDICT:\n" + verdict)
    res["verdict"] = verdict
    res["residuals"] = (
        "2D lattice CA with a SHARP tumour/normal boundary (real tumours are infiltrative/graded -> the rim and "
        "containment are softer); b_t, b_n, p_ig_n, f are effective probabilities; the recognition gate is modelled "
        "as a reduced bond prob into normal cells, not a molecular circuit; no 3D / vasculature / immune clearance. "
        "Robust claim = the b_n<p_c containment REQUIREMENT + the leak-amplification, not absolute fractions. Wet "
        "residual: the actual recognition-gating of the bystander signal (does the death factor require the "
        "receiver's mutation-sensor?) -- the molecular basis of b_n << b_t.")
    json.dump(res, open(os.path.join(OUT, "wave_containment.json"), "w"), indent=2)
    print(f"\n[saved] {OUT}/wave_containment.json")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "selftest":
        selftest()
    else:
        main()
