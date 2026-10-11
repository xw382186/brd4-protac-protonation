#!/usr/bin/env python
"""RDKit-based pKa estimate for the basic amine of BRD4-PROTAC compounds 5-8.

Method: linear free-energy relationship (LFER). Aliphatic amine basicity is
controlled by (a) degree of N-substitution and (b) the inductive pull of nearby
electronegative atoms, transmitted through sigma bonds with a fall-off factor.
Both are computed from RDKit connectivity only; the regression coefficients are
fitted on amines with tabulated experimental pKa, and validated leave-one-out.

  sigma_ind(N) = SUM_i (chi_i - chi_C) * alpha^(d_i - 1)

  d_i = topological bond distance from the basic N to atom i
  chi = Pauling electronegativity, alpha = 0.4 (standard inductive fall-off)

No per-group empirical constants are introduced: alpha is the only transferred
literature constant, everything else is fitted and reported with its error.
"""
import numpy as np, json, csv
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

ALPHA = 0.40                        # inductive fall-off per sigma bond
CHI = {"C":2.55,"N":3.04,"O":3.44,"S":2.58,"F":3.98,"Cl":3.16,"Br":2.96,"I":2.66,"H":2.20}
CHI_C = CHI["C"]

# ---- reference set: conjugate-acid pKa, standard tabulated values -----------
# (widely tabulated textbook/CRC-type values; pin a primary source before publication)
REF = [
 ("ammonia",                  "N",                      9.25),
 ("methylamine",              "CN",                    10.66),
 ("ethylamine",               "CCN",                   10.70),
 ("dimethylamine",            "CNC",                   10.73),
 ("trimethylamine",           "CN(C)C",                 9.80),
 ("triethylamine",            "CCN(CC)CC",             10.75),
 ("pyrrolidine",              "C1CCNC1",               11.27),
 ("piperidine",               "C1CCNCC1",              11.12),
 ("azetidine",                "C1CNC1",                11.29),
 ("N-methylpyrrolidine",      "CN1CCCC1",              10.32),
 ("N-methylpiperidine",       "CN1CCCCC1",             10.08),
 ("morpholine",               "C1COCCN1",               8.36),
 ("N-methylmorpholine",       "CN1CCOCC1",              7.38),
 ("ethanolamine",             "NCCO",                   9.50),
 ("2-methoxyethylamine",      "NCCOC",                  9.40),
 ("diethanolamine",           "OCCNCCO",                8.88),
 ("triethanolamine",          "OCCN(CCO)CCO",           7.76),
 ("ethylenediamine (pKa1)",   "NCCN",                   9.92),
 ("glycinamide",              "NCC(N)=O",               8.20),
 ("3-aminopropanamide",       "NCCC(N)=O",              9.30),
]

BASIC = Chem.MolFromSmarts(
  "[NX3;!$(N-a);!$(N[CX3]=[OX1,SX1]);!$(N[SX4](=O)=O);!$(N=*);!$(N#*);!$(N[NX3,NX2])]")

def basic_nitrogens(m):
    return [i for (i,) in m.GetSubstructMatches(BASIC)]

def ring_hetero(m, idx):
    """Heteroatoms sharing a ring with the basic N. In a ring the inductive pull
    reaches N along BOTH directions, which a shortest-path sigma term counts once."""
    ri = m.GetRingInfo(); D = rdmolops.GetDistanceMatrix(m)
    tot = 0.0
    for r in ri.AtomRings():
        if idx not in r: continue
        for j in r:
            if j == idx: continue
            sym = m.GetAtomWithIdx(j).GetSymbol()
            if sym in CHI and abs(CHI[sym]-CHI_C) > 1e-9:
                tot += (CHI[sym]-CHI_C) * (ALPHA ** (D[idx][j] - 1))
    return tot

def descriptors(smi, pick=None):
    """Return descriptors for the basic N."""
    m = Chem.MolFromSmiles(smi)
    if m is None: raise ValueError(smi)
    cands = basic_nitrogens(m)
    if not cands: raise ValueError("no basic N in "+smi)
    # for the reference set take the most basic-looking (highest substitution, least EWG);
    # for the PROTACs there is exactly one
    idx = cands[0] if pick is None else pick
    if pick is None and len(cands) > 1:
        # ethylenediamine / diamines: take the N with the smaller inductive pull
        idx = min(cands, key=lambda k: _sigma(m, k))
    mh = Chem.AddHs(m)
    AllChem.ComputeGasteigerCharges(mh)
    # map idx (heavy-atom index is preserved by AddHs)
    q = float(mh.GetAtomWithIdx(idx).GetProp("_GasteigerCharge"))
    a = m.GetAtomWithIdx(idx)
    n_subst = a.GetDegree()                       # heavy neighbours = 1/2/3
    ri = m.GetRingInfo()
    rings = [len(r) for r in ri.AtomRings() if idx in r]
    ring_size = min(rings) if rings else 0
    return n_subst, _sigma(m, idx), ring_size, q, idx, len(cands), ring_hetero(m, idx)

def _sigma(m, idx):
    D = rdmolops.GetDistanceMatrix(m)
    s = 0.0
    for at in m.GetAtoms():
        j = at.GetIdx()
        if j == idx: continue
        sym = at.GetSymbol()
        if sym not in CHI: continue
        dchi = CHI[sym] - CHI_C
        if abs(dchi) < 1e-9: continue
        d = D[idx][j]
        if d < 1 or d > 8: continue
        s += dchi * (ALPHA ** (d - 1))
    return s

# ---- build design matrix ----------------------------------------------------
rows=[]
for name, smi, pka in REF:
    n_subst, sig, ring, q, idx, ncand, rhet = descriptors(smi)
    rows.append(dict(name=name, smi=smi, pka=pka, n_subst=n_subst,
                     sigma=sig, ring=ring, q=q, rhet=rhet))

def design(r):
    # physically motivated, deliberately small basis
    return [1.0, r["sigma"], float(r["n_subst"] == 2), float(r["n_subst"] == 3),
            r.get("rhet", 0.0)]

X = np.array([design(r) for r in rows])
y = np.array([r["pka"] for r in rows])
beta, *_ = np.linalg.lstsq(X, y, rcond=None)
pred = X @ beta
resid = y - pred
r2 = 1 - (resid**2).sum() / ((y - y.mean())**2).sum()
rmse = float(np.sqrt((resid**2).mean()))

# leave-one-out
loo = []
for k in range(len(rows)):
    keep = [i for i in range(len(rows)) if i != k]
    b, *_ = np.linalg.lstsq(X[keep], y[keep], rcond=None)
    loo.append(float(y[k] - X[k] @ b))
loo = np.array(loo)
loo_rmse = float(np.sqrt((loo**2).mean()))
loo_mae  = float(np.abs(loo).mean())

print("="*84)
print("A. 标定集拟合（RDKit 描述符 -> 实验 pKa）")
print("="*84)
print(f"  模型: pKa = b0 + b1*sigma_ind + b2*[2deg] + b3*[3deg] + b4*ring_hetero   (alpha={ALPHA})")
names=["b0(intercept)","b1(sigma_ind)","b2(secondary)","b3(tertiary)","b4(ring hetero)"]
for n,b in zip(names,beta): print(f"    {n:18s} = {b:+8.3f}")
print(f"  n={len(rows)}  R^2={r2:.3f}  RMSE={rmse:.3f}")
print(f"  留一法: LOO-RMSE={loo_rmse:.3f}  LOO-MAE={loo_mae:.3f}  (这就是预测的实际不确定度)")
print()
print(f"  {'化合物':26s}{'实验':>7}{'拟合':>8}{'残差':>8}{'LOO残差':>9}")
for r,p,e,l in zip(rows,pred,resid,loo):
    print(f"  {r['name']:26s}{r['pka']:7.2f}{p:8.2f}{e:+8.2f}{l:+9.2f}")

# ---- apply to compounds 5-8 -------------------------------------------------
SM = {
 "5":"ClC1=CC=C(C2=N[C@H](C(N3C4=C2C(C)=C(C)S4)=NN=C3C)CC(NCCN5CC6(C5)CC(C6)C7=CC=C(CN([C@@H]8C(NC(CC8)=O)=O)C9=O)C9=C7)=O)C=C1",
 "6":"ClC1=CC=C(C2=N[C@H](C(N3C4=C2C(C)=C(C)S4)=NN=C3C)CC(NCCCN5CC6(C5)CC(C6)C7=CC=C(CN([C@@H]8C(NC(CC8)=O)=O)C9=O)C9=C7)=O)C=C1",
 "7":"CC1=C(C2=C(N3C(C)=NN=C3[C@H](CC(NCCN4CCC5(CC4)CN(C(C=C67)=CC=C7CN([C@H]8CCC(NC8=O)=O)C6=O)C5)=O)N=C2C9=CC=C(C=C9)Cl)S1)C",
 "8":"ClC1=CC=C(C2=N[C@H](C(N3C4=C2C(C)=C(C)S4)=NN=C3C)CC(NCCCN5CCC6(CN(C6)C7=CC=C(CN([C@@H]8C(NC(CC8)=O)=O)C9=O)C9=C7)CC5)=O)C=C1",
}
OLD_ANCHOR = {"5":8.5,"6":9.3,"7":8.6,"8":9.4}
print()
print("="*84)
print("B. 外推到化合物 5-8")
print("="*84)
print(f"  {'cpd':<5}{'碱性N数':>8}{'取代度':>7}{'环':>4}{'sigma_ind':>11}{'算得pKa':>10}{'+-LOO':>8}{'原锚点':>8}{'差':>7}")
out={}
for cid,smi in SM.items():
    n_subst, sig, ring, q, idx, ncand, rhet = descriptors(smi)
    r=dict(sigma=sig,n_subst=n_subst,ring=ring,rhet=rhet)
    p=float(np.array(design(r)) @ beta)
    f=lambda x: 100/(1+10**(7.4-x))
    out[cid]=dict(n_basic_N=ncand,n_subst=n_subst,ring=ring,sigma_ind=round(sig,4),ring_hetero=round(rhet,4),
                  pka_calc=round(p,2),loo_rmse=round(loo_rmse,2),
                  pka_lo=round(p-loo_rmse,2),pka_hi=round(p+loo_rmse,2),
                  frac_prot_pct=round(f(p),2),
                  frac_prot_lo_pct=round(f(p-loo_rmse),2),
                  old_anchor=OLD_ANCHOR[cid])
    print(f"  {cid:<5}{ncand:>8}{n_subst:>7}{ring:>4}{sig:>11.4f}{p:>10.2f}{loo_rmse:>8.2f}"
          f"{OLD_ANCHOR[cid]:>8.1f}{p-OLD_ANCHOR[cid]:>+7.2f}")
print()
print("="*84)
print("C. 质子化比例与判定余量（pH 7.4）")
print("="*84)
print(f"  {'cpd':<5}{'算得pKa':>10}{'f(+1)%':>9}{'下界pKa':>10}{'下界f%':>9}{'距7.4余量':>11}")
for cid,d in out.items():
    print(f"  {cid:<5}{d['pka_calc']:>10.2f}{d['frac_prot_pct']:>9.2f}{d['pka_lo']:>10.2f}"
          f"{d['frac_prot_lo_pct']:>9.2f}{d['pka_lo']-7.4:>11.2f}")
print()
print("  判定 +1 只需 pKa > 7.4。上表『下界』已扣掉一个 LOO-RMSE。")
json.dump(dict(alpha=ALPHA,beta=list(map(float,beta)),r2=r2,rmse=rmse,
               loo_rmse=loo_rmse,loo_mae=loo_mae,
               reference=[{k:r[k] for k in ("name","smi","pka","n_subst","sigma","ring","rhet")} for r in rows],
               compounds=out), open("pka_rdkit_lfer.json","w"), indent=1)
with open("pka_rdkit_lfer.csv","w",newline="") as fh:
    w=csv.writer(fh)
    w.writerow(["compound","n_basic_N","n_subst_on_N","ring_size","sigma_ind",
                "pKa_calc","LOO_RMSE","pKa_lower","pKa_upper","frac_prot_pct",
                "frac_prot_at_lower_pct","old_hand_anchor"])
    for cid,d in out.items():
        w.writerow([cid,d["n_basic_N"],d["n_subst"],d["ring"],d["sigma_ind"],d["pka_calc"],
                    d["loo_rmse"],d["pka_lo"],d["pka_hi"],d["frac_prot_pct"],
                    d["frac_prot_lo_pct"],d["old_anchor"]])
print("\n[written] pka_rdkit_lfer.json  pka_rdkit_lfer.csv")
