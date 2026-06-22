#!/usr/bin/env python3
"""
RUNG-40 — INDEPENDENCE IS THE CURRENCY: the microenvironment sniffer (Shriya's H5) as the AND-gate's INDEPENDENT
input. Answers R38's open need ("N>=3 INDEPENDENT recognizers") and R35's binding constraint (leak correlation).

THE GAP
-------
R35 showed the AND-gate's normal-cell safety lives or dies on the leak CORRELATION rho between recognizers; R38
showed the safe-and-curative region needs N>=3 INDEPENDENT recognizers. But where do independent recognizers come
FROM? The mutation-sensors (allele-specific CRISPR R27 / RNA toehold R25) all read the SAME signal class -- somatic
nucleic-acid changes -- so their leaks are plausibly CORRELATED (shared delivery vehicle, shared transcriptional
bursting, a shared stressed cell-state that elevates all of them at once). You therefore CANNOT simply stack three
mutation sensors and assume rho=0. The fix is ORTHOGONAL signal MODALITIES that fail independently:
   mutation (DNA/RNA, R27/25)  +  oncometabolite (2-HG, R33)  +  microenvironment (hypoxia/acidosis/ROS, H5).
The microenvironment LEAKS as a standalone gate (hypoxia/low-pH/ROS also occur in normal exercising, inflamed, or
ischaemic tissue -- like every non-mutation window, R15/R23). H5's real value is NOT as a standalone recognizer but
as the AND-gate's INDEPENDENCE-PROVIDER: a different physical signal read by a different mechanism -> its leak is
uncorrelated with the mutation-sensor leaks -> it multiplies its full leak factor into the joint, which correlated
same-modality sensors do not.

THE TEST (block-correlated Gaussian copula -> multivariate-normal orthant probability)
Compare AND-gate normal false-death across designs that differ in correlation STRUCTURE, and show:
  a LEAKY-but-INDEPENDENT microenvironment input beats a CLEAN-but-CORRELATED extra mutation input.

HONEST CEILING: per-sensor leaks and the within/cross correlations are representative parameters, NOT measured
(the microenvironment standalone leak ~20-40% is literature-typical for hypoxia markers, swept here; the
mutation-sensor cross-correlation is the R35 wet residual). Robust claim = the RANKING (independence dominates
individual cleanliness in an AND-gate) and the design conclusion (build the gate from orthogonal modalities).

USAGE
  python scripts/74_independence_andgate.py selftest
  python scripts/74_independence_andgate.py run   # -> runs/rung40_independence/ (CPU, seconds)
"""
import os, sys, json
import numpy as np
from scipy.stats import norm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "runs/rung40_independence")


def block_corr(sizes, rho_within, rho_cross):
    """Correlation matrix: sensors grouped in blocks; rho_within inside a block, rho_cross between blocks."""
    n = sum(sizes)
    blk = np.concatenate([[b] * s for b, s in enumerate(sizes)])
    S = np.full((n, n), rho_cross, float)
    for i in range(n):
        for j in range(n):
            if i == j:
                S[i, j] = 1.0
            elif blk[i] == blk[j]:
                S[i, j] = rho_within
    return S


def is_psd(S, tol=1e-8):
    return np.linalg.eigvalsh(S).min() > -tol


def joint_fire(leaks, S, n=4_000_000, seed=0, chunk=1_000_000):
    """P(all sensors fire in a normal cell): each sensor i fires if latent X_i > t_i (t_i=Phi^-1(1-L_i)),
    X ~ N(0, S). Chunked Monte-Carlo (robust to any correlation structure)."""
    leaks = np.asarray(leaks, float)
    t = norm.ppf(1.0 - leaks)
    L = np.linalg.cholesky(S)
    rng = np.random.default_rng(seed)
    hits, done = 0, 0
    while done < n:
        m = min(chunk, n - done)
        Z = rng.standard_normal((m, len(leaks)))
        X = Z @ L.T
        hits += int(np.all(X > t, axis=1).sum())
        done += m
    return hits / n


def selftest():
    # independent (identity) -> joint == product of leaks
    S = np.eye(3)
    j = joint_fire([0.1, 0.1, 0.1], S, n=4_000_000, seed=1)
    assert abs(j - 0.001) < 1.5e-4, j                      # 0.1^3 = 1e-3
    # adding an INDEPENDENT input multiplies by its leak
    S2 = np.eye(2); S3 = np.eye(3)
    j2 = joint_fire([0.1, 0.1], S2, n=4_000_000, seed=2)
    j3 = joint_fire([0.1, 0.1, 0.3], S3, n=4_000_000, seed=3)
    assert abs(j3 - j2 * 0.3) < 1.5e-4, (j2, j3)
    # correlation INCREASES the joint leak (worse) vs independent, same leaks
    Sc = block_corr([2], rho_within=0.6, rho_cross=0.0)
    jc = joint_fire([0.1, 0.1], Sc, n=4_000_000, seed=4)
    assert jc > j2 + 5e-4, (j2, jc)
    assert is_psd(block_corr([2, 1], 0.5, 0.0)) and is_psd(block_corr([3], 0.5, 0.0))
    print(f"[selftest] copula orthant OK (indep 0.1^3={j:.2e}; +indep x0.3 -> {j3:.2e}; "
          f"corr {jc:.2e} > indep {j2:.2e})")


def main():
    os.makedirs(OUT, exist_ok=True)
    selftest()
    res = {"tag": "rung40_independence", "hypothesis": "H5 microenvironment as independent AND-input"}
    Lmut, L2hg, Ltme = 0.05, 0.02, 0.30   # representative standalone leaks (TME leakiest; the wet residual)
    RW = 0.5                              # within-modality (mutation-block) leak correlation (R35 residual)
    N = 8_000_000

    # 2 mutation sensors, correlated within the mutation modality (the realistic starting AND-gate)
    base2 = joint_fire([Lmut, Lmut], block_corr([2], RW, 0.0), n=N, seed=10)
    print(f"\n=== RUNG-40 (H5): independence is the currency of an AND-gate ===")
    print(f"start: 2 mutation-sensors (L={Lmut}, within-correlation rho={RW}) -> joint normal false-death {base2:.2e}")

    print(f"\n-- add a 3rd input: which kind helps most? (all vs the SAME 2-mutation base) --")
    options = [
        ("3rd MUTATION sensor, CORRELATED (rho=0.5, L=0.05)", block_corr([3], RW, 0.0), [Lmut, Lmut, Lmut]),
        ("3rd MUTATION sensor, hypothetically INDEPENDENT (L=0.05)", block_corr([2, 1], RW, 0.0), [Lmut, Lmut, Lmut]),
        ("MICROENVIRONMENT sensor, INDEPENDENT (L=0.30, leaky!)", block_corr([2, 1], RW, 0.0), [Lmut, Lmut, Ltme]),
        ("2-HG metabolite sensor, INDEPENDENT (L=0.02)", block_corr([2, 1], RW, 0.0), [Lmut, Lmut, L2hg]),
    ]
    rows = []
    print("   3rd input                                                  joint false-death   vs base")
    for name, S, leaks in options:
        j = joint_fire(leaks, S, n=N, seed=hash(name) % 1000)
        rows.append({"third_input": name, "joint_false_death": j, "fold_vs_base": j / base2})
        print(f"   {name:<55} {j:.2e}       {j/base2:.2f}x")
    res["base_2mut"] = base2
    res["add_third_input"] = rows

    # the punchline comparison
    corr_mut = next(r for r in rows if "CORRELATED" in r["third_input"])
    tme = next(r for r in rows if "MICROENVIRONMENT" in r["third_input"])
    print(f"\n   PUNCHLINE: a 30%-leaky INDEPENDENT microenvironment input ({tme['joint_false_death']:.2e}) "
          f"BEATS a 5%-clean CORRELATED mutation input ({corr_mut['joint_false_death']:.2e}) "
          f"by {corr_mut['joint_false_death']/tme['joint_false_death']:.1f}x.")
    res["punchline_independent_leaky_beats_correlated_clean"] = \
        corr_mut["joint_false_death"] / tme["joint_false_death"]

    # full orthogonal-modality design vs all-correlated-mutation design
    print(f"\n-- whole 3-input gate: orthogonal MODALITIES vs same-modality stack --")
    designs = [
        ("3x mutation, correlated rho=0.5", block_corr([3], RW, 0.0), [Lmut, Lmut, Lmut]),
        ("3x mutation, correlated rho=0.7", block_corr([3], 0.7, 0.0), [Lmut, Lmut, Lmut]),
        ("mutation + 2-HG + microenvironment (orthogonal, independent)", np.eye(3), [Lmut, L2hg, Ltme]),
    ]
    drows = []
    print("   design                                                        joint false-death")
    for name, S, leaks in designs:
        j = joint_fire(leaks, S, n=N, seed=(hash(name) % 1000) + 500)
        drows.append({"design": name, "joint_false_death": j})
        print(f"   {name:<58} {j:.2e}")
    res["whole_gate_designs"] = drows

    # sweep: as mutation-sensors get MORE correlated, the orthogonal design's advantage grows
    print(f"\n-- as same-modality correlation rises, the orthogonal-modality design wins by more --")
    print("   rho_within(mutation)   3x-mutation joint   orthogonal joint   advantage")
    ortho = joint_fire([Lmut, L2hg, Ltme], np.eye(3), n=N, seed=999)
    sweep = []
    for rw in (0.0, 0.3, 0.5, 0.7, 0.9):
        jm = joint_fire([Lmut, Lmut, Lmut], block_corr([3], rw, 0.0), n=N, seed=int(rw * 100) + 1)
        sweep.append({"rho_within": rw, "mut3_joint": jm, "ortho_joint": ortho, "advantage": jm / ortho})
        print(f"        {rw:.1f}              {jm:.2e}          {ortho:.2e}        {jm/ortho:6.1f}x")
    res["correlation_sweep"] = sweep

    verdict = (
        "INDEPENDENCE, not individual cleanliness, is the currency of an AND-gate -- and that re-frames Shriya's "
        "microenvironment sniffer (H5). (1) The mutation-sensors that carry our recognition (allele-CRISPR R27, RNA "
        "toehold R25) all read the SAME signal class, so their leaks are plausibly CORRELATED; stacking three of "
        "them does NOT give the independent N=3 that R38 requires -- a correlated 3rd mutation sensor barely lowers "
        "the joint false-death. (2) A microenvironment sensor LEAKS badly as a standalone gate (~30%, since normal "
        "hypoxic/inflamed/ischaemic tissue shares the signal) -- yet as an AND INPUT it multiplies its FULL leak "
        "factor in BECAUSE it is independent, and a 30%-leaky INDEPENDENT microenvironment input beats a 5%-clean "
        "CORRELATED mutation input. (3) The winning design is therefore an AND of ORTHOGONAL MODALITIES -- mutation "
        "(DNA/RNA) + oncometabolite (2-HG, R33) + microenvironment (hypoxia/acidosis, H5) -- whose joint false-death "
        "is the near-full product of the individual leaks (~3e-4 here), and whose advantage over a same-modality "
        "stack GROWS as the same-modality sensors become more correlated. This answers R38's open question (where do "
        ">=3 INDEPENDENT recognizers come from -> orthogonal physical modalities, not more of the same), gives H5 a "
        "concrete role (the independence-provider, not a standalone recognizer), and turns R35's #1 residual (leak "
        "correlation) into a DESIGN PRINCIPLE: pick recognizers that fail for unrelated physical reasons.")
    print("\nVERDICT:\n" + verdict)
    res["verdict"] = verdict
    res["residuals"] = (
        "Per-sensor leaks (Lmut 0.05, L2hg 0.02, Ltme 0.30) and the within/cross correlations are representative, "
        "NOT measured: the microenvironment standalone leak is literature-typical for hypoxia markers (swept), and "
        "the mutation-sensor CROSS-correlation is exactly R35's wet residual (must be measured in stressed normal "
        "cells). Gaussian-copula leak model (a modelling choice). Robust claim = the RANKING (independent>correlated "
        "in an AND-gate) and the orthogonal-modality design principle, not absolute fractions. The microenvironment "
        "sensor itself (a conformational pH/hypoxia switch coupled to the circuit) is an un-built design (cf R33).")
    json.dump(res, open(os.path.join(OUT, "independence.json"), "w"), indent=2)
    print(f"\n[saved] {OUT}/independence.json")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "selftest":
        selftest()
    else:
        main()
