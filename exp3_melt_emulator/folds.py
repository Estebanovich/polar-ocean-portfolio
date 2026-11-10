"""Pre-registered cross-validation folds for Exp. 3.

Training pool: the 25 runs with group in {factorial, geometry}. The reference run
(``ocean0_4km``, different thermocline thickness) and the split == "test" runs
(Tbot = 2.0) never train. Only case, group, split, Tbot, z_c, s and n1_pass of a
run are used.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

ZC_LOW = 325.0    # zc_out_250: z_c < 325 m
ZC_HIGH = 475.0   # zc_out_550: z_c > 475 m; zc_out_400: 325 <= z_c <= 475 m
WARM_TRAIN_MAX = 1.0   # warm_out train: Tbot <= 1.0
WARM_TEST_MIN = 1.4    # warm_out test: Tbot >= 1.4


@dataclass(frozen=True)
class Fold:
    name: str
    kind: str   # "loro" | "zc_out" | "warm_out" | "profile_test" | "final"
    train: tuple[str, ...]
    test: tuple[str, ...]


def make_folds(runs) -> list[Fold]:
    pool = [r for r in runs if r.group != "reference" and r.split != "test"]
    all_train = tuple(r.case for r in pool)
    folds = [Fold(f"loro_{r.case}", "loro",
                  tuple(c for c in all_train if c != r.case), (r.case,))
             for r in pool]

    def zc_fold(name, pred):
        held = tuple(r.case for r in pool if pred(r.z_c))
        return Fold(name, "zc_out", tuple(c for c in all_train if c not in held), held)

    folds += [
        zc_fold("zc_out_250", lambda z: z < ZC_LOW),
        zc_fold("zc_out_400", lambda z: ZC_LOW <= z <= ZC_HIGH),
        zc_fold("zc_out_550", lambda z: z > ZC_HIGH),
    ]
    folds.append(Fold(
        "warm_out", "warm_out",
        tuple(r.case for r in pool if r.Tbot <= WARM_TRAIN_MAX),
        tuple(r.case for r in pool if r.Tbot >= WARM_TEST_MIN)))
    ref = tuple(r.case for r in runs if r.group == "reference")
    folds.append(Fold("profile_test", "profile_test", all_train, ref))
    folds.append(Fold("final", "final", all_train,
                      tuple(r.case for r in runs if r.split == "test")))
    return folds


def drop_n1(fold: Fold, runs) -> Fold:
    """Sensitivity variant: remove N1-failing runs from ``train`` (test unchanged)."""
    bad = {r.case for r in runs if not r.n1_pass}
    return replace(fold, train=tuple(c for c in fold.train if c not in bad))
