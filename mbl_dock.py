#!/usr/bin/env python3
# ====
#  mbl_dock.py  —  Docking di-nuclear para Metalo-beta-lactamases (MBL)
# ----
#  Standalone. Metodo padrao-ouro AutoDock4-Zn (pseudo-atomos TZ + AD4Zn.dat).
#  Motor de sampling selecionavel: Vina>=1.2 (--scoring ad4)  OU  autodock4.
#  Ligante: RDKit (SMILES->3D, ETKDGv3+MMFF94, preserva quiralidade) -> obabel.
#  Receptor: prepare_receptor4.py (MGLTools/python2.7) -> Gasteiger por residuo
#            (NAO depende de kekulizacao) -> reinsere Zn com carga +2.0.
#  BRIDGE FILTER (v2, CORRIGIDO): rejeita poses que nao coordenam AMBOS os
#            zincos POR HETEROATOMO COORDENANTE (S/N/O). Carbono e hidrogenio
#            nao coordenam Zn2+ -> alcanos sao rejeitados automaticamente.
#            Tipos AutoDock sao traduzidos p/ elemento (OA->O, NA->N, SA->S,
#            A->C, HD->H), evitando falsos negativos.
#  MODO --redock: redocka o ligante co-cristalizado no proprio receptor e
#            mede o RMSD contra a pose depositada (controle positivo do
#            protocolo; criterio RMSD <= 2.0 A na pose de melhor escore).
#  MODO --reanalyze: re-aplica o filtro corrigido sobre out_*.pdbqt ja
#            existentes, SEM re-executar docking (incorpora o antigo
#            reanalyze_bridges.py).
#  Sem falhas silenciosas. Disclaimer: docking nao-covalente; resultados
#  qualitativos/exploratorios (organossulfurados atuam tb via redox/covalente).
# ====

import os, re, sys, csv, math, shutil, logging, argparse, subprocess
from pathlib import Path
from itertools import combinations
from typing import NamedTuple, Optional, List
import numpy as np

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    _HAS_RDKIT = True
except Exception:
    _HAS_RDKIT = False

logging.basicConfig(level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")
log = logging.getLogger("mbl_dock")


# ====
#  0. QUIMICA DE COORDENACAO — constantes do bridge filter
# ====
# Apenas estes elementos coordenam Zn2+ em metaloproteinas.
# (S-Zn tiolato ~2.3 A; N-Zn His ~2.1 A; O-Zn carboxilato/agua ~2.0 A)
COORDINATING_ELEMENTS = frozenset({"S", "N", "O"})
DEFAULT_BRIDGE_CUT = 3.0          # A — limite generoso p/ coordenacao
MAX_DIZN_DIST = 5.0               # A — acima disso o par Zn-Zn nao e' di-nuclear

# Tipo atomico AutoDock (colunas 78-79 do PDBQT) -> elemento quimico.
# CRITICO: 'OA'/'NA'/'SA' sao aceptores de H-bond (O/N/S) e PRECISAM ser
# reconhecidos, senao quase todo oxigenio/nitrogenio do ligante e' ignorado.
AD_TYPE_TO_ELEMENT = {
    "C": "C", "A": "C",
    "N": "N", "NA": "N", "NS": "N",
    "O": "O", "OA": "O", "OS": "O",
    "S": "S", "SA": "S",
    "H": "H", "HD": "H", "HS": "H",
    "F": "F", "CL": "CL", "BR": "BR", "I": "I", "P": "P",
    "MG": "MG", "MN": "MN", "FE": "FE", "CA": "CA", "ZN": "ZN",
    "TZ": "TZ",   # pseudo-atomo tetraedrico do AutoDock4-Zn
}
_TWO_LETTER = ("CL", "BR", "ZN", "MG", "MN", "FE", "TZ")


def element_of(ln):
    """Elemento quimico de uma linha ATOM/HETATM de PDB ou PDBQT.

    CUIDADO com as colunas: no PDB o campo 'element' ocupa as colunas 77-78
    (indices 76:78), enquanto no PDBQT o tipo AutoDock ocupa as colunas 78-79
    (indices 77:79). Ler o offset errado transforma 'CL' em 'L' e faz o cloro
    ser tratado como elemento desconhecido. Tentamos os dois offsets e ficamos
    com o primeiro que corresponda a um elemento/tipo conhecido; se nenhum
    servir, caimos para o nome do atomo (colunas 13-16)."""
    hits = []
    for a, b in ((77, 79), (76, 78)):             # PDBQT primeiro, depois PDB
        tok = ln[a:b].strip().upper() if len(ln) > a else ""
        if not tok:
            continue
        clean = re.sub(r"[+\-0-9]", "", tok)      # ex.: 'N1' -> 'N', 'CL-' -> 'CL'
        for t in (tok, clean):
            if t in AD_TYPE_TO_ELEMENT:
                hits.append(t); break
    if hits:
        # 'ZN' lido com o offset do PDBQT vira 'N ' -> 'N', que tambem e' valido;
        # o token de dois caracteres e' sempre o mais especifico e vence o empate.
        best = max(hits, key=lambda t: (len(t) == 2 and t in _TWO_LETTER, len(t)))
        return AD_TYPE_TO_ELEMENT[best]
    nm = ln[12:16].strip().upper().lstrip("0123456789")
    if not nm:
        return "?"
    if nm[:2] in _TWO_LETTER:
        return nm[:2]
    return nm[0] if nm else "?"


class Pose(NamedTuple):
    """Uma pose de docking: energia (kcal/mol), coordenadas Nx3 e elementos."""
    energy: Optional[float]
    coords: np.ndarray
    elements: List[str]


def parse_coord_elements(spec):
    """'S,N,O' -> frozenset; 'any'/'all'/'' -> None (desliga o filtro elementar)."""
    if spec is None:
        return COORDINATING_ELEMENTS
    s = spec.strip().lower()
    if s in ("any", "all", "todos", "*", ""):
        return None
    els = frozenset(x.strip().upper() for x in spec.split(",") if x.strip())
    return els or COORDINATING_ELEMENTS


# ====
#  1. PREPARO DE LIGANTE — RDKit (SMILES->3D) -> obabel -> PDBQT
# ====
def _rdkit_3d_pdb_block(smiles, seed=0xF00D):
    if not _HAS_RDKIT:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        log.error("  RDKit nao parseou SMILES: %s", smiles); return None
    mol = Chem.AddHs(mol)
    p = AllChem.ETKDGv3(); p.randomSeed = seed
    p.useSmallRingTorsions = True; p.enforceChirality = True
    if AllChem.EmbedMolecule(mol, p) != 0:
        p.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, p) != 0:
            log.error("  RDKit falhou no embedding 3D: %s", smiles); return None
    try:
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
        else:
            log.warning("  MMFF sem params p/ %s — usando UFF.", smiles)
            AllChem.UFFOptimizeMolecule(mol, maxIters=1000)
    except Exception as e:
        log.warning("  Otimizacao FF falhou (%s): %s — geometria bruta.", smiles, e)
    try:
        Chem.AssignStereochemistryFrom3D(mol)
    except Exception as e:
        log.warning("  AssignStereochemistryFrom3D falhou (%s): %s", smiles, e)
    return Chem.MolToPDBBlock(mol)


def prepare_ligand(smiles, out_pdbqt, obabel="obabel", ph=7.4):
    blk = _rdkit_3d_pdb_block(smiles)
    if blk:
        cmd = [obabel, "-ipdb", "-opdbqt", "-O", out_pdbqt,
               "-p", str(ph), "--partialcharge", "gasteiger"]
        try:
            r = subprocess.run(cmd, input=blk, capture_output=True, text=True, timeout=120)
            if r.returncode == 0 and os.path.exists(out_pdbqt) and os.path.getsize(out_pdbqt) > 0:
                return True
            log.error("  obabel(RDKit->PDBQT) falhou (%s): %s", smiles, r.stderr[:200])
        except Exception as e:
            log.error("  obabel(RDKit->PDBQT) excecao: %s", e)
    log.warning("  Fallback obabel --gen3d (%s) — menos confiavel.", smiles)
    cmd = [obabel, f"-:{smiles}", "-O", out_pdbqt, "--gen3d",
           "-p", str(ph), "--partialcharge", "gasteiger"]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if r.returncode == 0 and os.path.exists(out_pdbqt) and os.path.getsize(out_pdbqt) > 0:
            return True
        log.error("  obabel --gen3d falhou (%s): %s", smiles, r.stderr[:200])
    except Exception as e:
        log.error("  obabel --gen3d excecao: %s", e)
    return False


# ====
#  2. PREPARO DE RECEPTOR — isola cadeia + 2 Zn cataliticos
# ====
def clean_receptor(pdb_in, chain, keep_zn, remove_het, out_pdb):
    """Isola a cadeia + Zn cataliticos. Remove confôrmeros alternativos (altloc):
    mantem apenas altloc em branco ou 'A' e zera o indicador de altloc."""
    kept_prot = 0; zn_xyz = []; keep_zn = set(str(z) for z in keep_zn)
    n_altloc = 0
    with open(pdb_in) as f, open(out_pdb, "w") as o:
        for ln in f:
            rec = ln[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            altloc = ln[16]
            if altloc not in (" ", "A"):
                n_altloc += 1
                continue                       # descarta confôrmeros B/C/...
            if altloc == "A":
                ln = ln[:16] + " " + ln[17:]   # zera o indicador de altloc
            res = ln[17:20].strip(); ch = ln[21].strip(); rs = ln[22:26].strip()
            if rec == "ATOM" and ch == chain:
                o.write(ln); kept_prot += 1
            elif res == "ZN" and ch == chain and rs in keep_zn:
                o.write(ln)
                zn_xyz.append(np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]))
    if n_altloc:
        log.info("  Confôrmeros alternativos descartados (altloc != A): %d", n_altloc)
    if kept_prot == 0:
        raise RuntimeError(f"clean_receptor: 0 atomos proteina na cadeia '{chain}' de {pdb_in}!")
    if len(zn_xyz) != len(keep_zn):
        raise RuntimeError(f"clean_receptor: esperava {len(keep_zn)} Zn, guardei {len(zn_xyz)}.")
    if len(zn_xyz) == 2:
        dist = float(np.linalg.norm(zn_xyz[0] - zn_xyz[1]))
        if not (2.8 <= dist <= MAX_DIZN_DIST):
            log.warning("  Zn1-Zn2 = %.2f A fora da faixa di-nuclear (2.8-%.1f).",
                        dist, MAX_DIZN_DIST)
        else:
            log.info("  Zn1-Zn2 = %.2f A (OK, di-nuclear).", dist)
    log.info("  Receptor limpo: %d atomos proteina + %d Zn -> %s", kept_prot, len(zn_xyz), out_pdb)
    return zn_xyz


def select_catalytic_zn_pair(zn_xyz):
    """Dos Zn mantidos, devolve o par catalitico = par de menor distancia
    dentro da faixa di-nuclear. Util quando a estrutura tem Zn superficial
    extra (ex.: VIM-1/2/7). Com exatamente 2 Zn, devolve-os inalterados."""
    if len(zn_xyz) < 2:
        raise RuntimeError("select_catalytic_zn_pair: menos de 2 Zn — "
                           "o filtro di-nuclear exige 2 zincos.")
    if len(zn_xyz) == 2:
        return list(zn_xyz)
    best, best_d = None, float("inf")
    for a, b in combinations(range(len(zn_xyz)), 2):
        d = float(np.linalg.norm(zn_xyz[a] - zn_xyz[b]))
        if d < best_d and d <= MAX_DIZN_DIST:
            best_d, best = d, [zn_xyz[a], zn_xyz[b]]
    if best is None:
        raise RuntimeError(f"select_catalytic_zn_pair: {len(zn_xyz)} Zn e nenhum par "
                           f"com distancia <= {MAX_DIZN_DIST} A (nao ha sitio di-nuclear).")
    log.info("  %d Zn presentes; par catalitico selecionado (%.2f A).", len(zn_xyz), best_d)
    return best


def _fix_zn_charge(pdbqt_path, charge=2.0):
    """prepare_receptor4 MANTEM os Zn com carga 0.000 (sem param Gasteiger).
    Sobrescreve APENAS a carga [66:76] e o tipo [77:79], preservando o
    alinhamento de colunas do PDBQT (ocupacao/B-factor ate a coluna 66).
    NAO duplica atomos."""
    out, n = [], 0
    for l in open(pdbqt_path):
        if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() == "ZN":
            base = l[:66].rstrip("\n").ljust(66)         # garante 66 chars
            l = base + ("%+.3f" % charge).rjust(10) + " Zn\n"  # [66:76] carga, [77:79] tipo
            n += 1
        out.append(l)
    if n == 0:
        raise RuntimeError("_fix_zn_charge: nenhum Zn encontrado em %s!" % pdbqt_path)
    open(pdbqt_path, "w").writelines(out)
    log.info("  Carga dos Zn corrigida (+%.1f, colunas alinhadas): %d atomos", charge, n)


def prepare_receptor(pdb_in, out_pdbqt, python2, prepare_receptor4, mgl_pckgs):
    """PDB -> PDBQT via prepare_receptor4.py rodado com o python2.7 REAL do MGLTools.
    Gasteiger atribuido por residuo (nao depende de kekulizacao, que trava o obabel
    em proteinas). Reinsere os 2 Zn com +2.0. Aborta se as cargas vierem zeradas."""
    if prepare_receptor4 is None:
        raise RuntimeError("prepare_receptor4 nao informado (use --prepare_receptor4).")
    env = dict(os.environ)
    if mgl_pckgs and os.path.isdir(mgl_pckgs):
        env["PYTHONPATH"] = mgl_pckgs + os.pathsep + env.get("PYTHONPATH", "")
    cmd = [python2, prepare_receptor4, "-r", pdb_in, "-o", out_pdbqt, "-A", "checkhydrogens"]
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if not os.path.exists(out_pdbqt) or os.path.getsize(out_pdbqt) == 0:
        raise RuntimeError(
            f"prepare_receptor4 falhou (arquivo nao gerado).\n"
            f"stdout: {r.stdout[:300]}\nstderr: {r.stderr[:300]}")
    # "no silent failures": cargas nao podem ser todas zero
    charged = False
    for l in open(out_pdbqt):
        if l.startswith("ATOM"):
            try:
                if float(l[66:76]) != 0.0:
                    charged = True; break
            except ValueError:
                continue
    if not charged:
        raise RuntimeError("prepare_receptor4 gerou PDBQT com cargas zeradas (todas 0.0)!")
    _fix_zn_charge(out_pdbqt)
    log.info("  Receptor PDBQT (prepare_receptor4, cargas OK): %s", out_pdbqt)


# ====
#  3. AUTODOCK4-Zn — pseudo-atomos TZ + mapas AD4Zn
# ====
def add_zinc_pseudo(rec_pdbqt, out_pdbqt, zinc_pseudo="zinc_pseudo.py", pythonsh="pythonsh"):
    cmd = [pythonsh, zinc_pseudo, "-r", rec_pdbqt, "-o", out_pdbqt]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if not os.path.exists(out_pdbqt) or os.path.getsize(out_pdbqt) == 0:
        raise RuntimeError("zinc_pseudo.py falhou (instale MGLTools).\n"
                           f"stdout: {r.stdout[:300]}\nstderr: {r.stderr[:300]}")
    ntz = sum(1 for l in open(out_pdbqt) if l[12:16].strip() == "TZ")
    if ntz == 0:
        raise RuntimeError("zinc_pseudo.py rodou mas nao inseriu nenhum TZ!")
    # sanity: cargas ainda presentes apos a insercao dos pseudo-atomos
    charged = any(l.startswith("ATOM") and _safe_charge(l) != 0.0 for l in open(out_pdbqt))
    if not charged:
        raise RuntimeError("receptor_zn.pdbqt sem cargas apos zinc_pseudo.py!")
    log.info("  Pseudo-atomos TZ inseridos: %d (~4 por Zn).", ntz)


def _safe_charge(l):
    try:
        return float(l[66:76])
    except (ValueError, IndexError):
        return 0.0


def _atom_types(pdbqt):
    types = []
    for l in open(pdbqt):
        if l.startswith(("ATOM", "HETATM")):
            t = l[77:79].strip()
            if t and t not in types: types.append(t)
    return types


def write_gpf(rec_pdbqt, lig_pdbqt, center, gpf, npts=40, spacing=0.375, ad4zn_dat="AD4Zn.dat"):
    lig_types = _atom_types(lig_pdbqt); rec_types = _atom_types(rec_pdbqt)
    cx, cy, cz = center
    with open(gpf, "w") as f:
        f.write(f"parameter_file {ad4zn_dat}\n")
        f.write(f"npts {npts} {npts} {npts}\n")
        f.write("gridfld receptor.maps.fld\n")
        f.write(f"spacing {spacing}\n")
        f.write(f"receptor_types {' '.join(rec_types)}\n")
        f.write(f"ligand_types {' '.join(lig_types)}\n")
        f.write(f"receptor {rec_pdbqt}\n")
        f.write(f"gridcenter {cx:.3f} {cy:.3f} {cz:.3f}\n")
        for t in lig_types: f.write(f"map receptor.{t}.map\n")
        f.write("elecmap receptor.e.map\ndsolvmap receptor.d.map\ndielectric -0.1465\n")
    log.info("  GPF escrito: %s (center %.3f %.3f %.3f)", gpf, cx, cy, cz)


def run_autogrid(gpf, autogrid="autogrid4"):
    r = subprocess.run([autogrid, "-p", gpf, "-l", "autogrid.log"], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"autogrid4 falhou: {r.stderr[:300]}")
    log.info("  AutoGrid OK.")


# ---- backend A: Vina >= 1.2 usando force field AD4Zn ----
def run_vina_ad4(lig_pdbqt, center, out_pdbqt, out_log, vina="vina",
                 npts=40, spacing=0.375, exhaustiveness=16):
    size = npts * spacing; cx, cy, cz = center
    cmd = [vina, "--scoring", "ad4", "--maps", "receptor",
           "--ligand", lig_pdbqt,
           "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}", "--center_z", f"{cz:.3f}",
           "--size_x", f"{size:.1f}", "--size_y", f"{size:.1f}", "--size_z", f"{size:.1f}",
           "--exhaustiveness", str(exhaustiveness), "--out", out_pdbqt]
    r = subprocess.run(cmd, capture_output=True, text=True)
    open(out_log, "w").write(r.stdout + "\n" + r.stderr)
    if not os.path.exists(out_pdbqt) or os.path.getsize(out_pdbqt) == 0:
        raise RuntimeError(f"vina(ad4) falhou: {r.stderr[:300]}")
    log.info("  Vina-AD4 OK -> %s", out_pdbqt)


def parse_vina_pdbqt_models(pdbqt):
    """PDBQT multi-modelo do Vina -> [Pose(energy, coords, elements), ...].
    Aceita arquivo com pose unica (sem MODEL/ENDMDL)."""
    models, cur, els, energy = [], [], [], None
    for ln in open(pdbqt):
        if ln.startswith("REMARK VINA RESULT"):
            try: energy = float(ln.split()[3])
            except Exception: energy = None
        elif ln.startswith("REMARK minimizedAffinity") or "Estimated Free Energy" in ln:
            nums = re.findall(r"-?\d+\.\d+", ln)
            if nums and energy is None:
                energy = float(nums[0])
        elif ln.startswith(("ATOM", "HETATM")):
            cur.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            els.append(element_of(ln))
        elif ln.startswith("ENDMDL"):
            if cur: models.append(Pose(energy, np.array(cur), els))
            cur, els, energy = [], [], None
    if cur and not models:                       # pose unica sem MODEL/ENDMDL
        models.append(Pose(energy, np.array(cur), els))
    return models


# ---- backend B: AutoDock4 classico (GA) ----
def write_dpf(lig_pdbqt, rec_pdbqt, dpf, ga_runs=20):
    lig_types = _atom_types(lig_pdbqt)
    with open(dpf, "w") as f:
        f.write("autodock_parameter_version 4.2\n")
        f.write(f"ligand_types {' '.join(lig_types)}\nfld receptor.maps.fld\n")
        for t in lig_types: f.write(f"map receptor.{t}.map\n")
        f.write("elecmap receptor.e.map\ndesolvmap receptor.d.map\n")
        f.write(f"move {lig_pdbqt}\nabout 0.0 0.0 0.0\n")
        f.write("tran0 random\nquaternion0 random\ndihe0 random\n")
        f.write("ga_pop_size 150\nga_num_evals 250000\nga_num_generations 27000\n")
        f.write(f"ga_run {ga_runs}\nrmstol 2.0\nset_ga\nset_sw1\nanalysis\n")
    log.info("  DPF escrito: %s (%d GA runs)", dpf, ga_runs)


def run_autodock(dpf, out_dlg, autodock="autodock4"):
    r = subprocess.run([autodock, "-p", dpf, "-l", out_dlg], capture_output=True, text=True)
    if not os.path.exists(out_dlg):
        raise RuntimeError(f"autodock4 falhou: {r.stderr[:300]}")
    log.info("  AutoDock OK -> %s", out_dlg)


def parse_dlg_models(dlg):
    """DLG do AutoDock4 -> [Pose(energy, coords, elements), ...]."""
    models, cur, els, energy = [], [], [], None
    for ln in open(dlg):
        if "Estimated Free Energy of Binding" in ln:
            try: energy = float(ln.split("=")[1].split("kcal")[0])
            except Exception: energy = None
        elif ln.startswith("DOCKED: ATOM") or ln.startswith("DOCKED: HETATM"):
            b = ln[8:]
            cur.append([float(b[30:38]), float(b[38:46]), float(b[46:54])])
            els.append(element_of(b))
        elif ln.startswith("DOCKED: ENDMDL"):
            if cur: models.append(Pose(energy, np.array(cur), els))
            cur, els, energy = [], [], None
    return models


# ====
#  5. BRIDGE FILTER (v2 — apenas heteroatomos coordenantes S/N/O)
# ====
def bridges_both_zn(pose, zn_xyz, cut=DEFAULT_BRIDGE_CUT,
                    coord_elements=COORDINATING_ELEMENTS):
    """A pose faz ponte entre os DOIS Zn?

    So contam atomos capazes de coordenar Zn2+ (S, N, O por padrao).
    C/H nao coordenam: uma pose puramente hidrocarboneto (alcano) e'
    rejeitada por construcao, mesmo que a distancia geometrica seja curta.

    Aceita um Pose ou, por retrocompatibilidade, um array Nx3 de coordenadas
    (nesse caso nao ha informacao de elemento e o filtro elementar e' ignorado).

    Retorna (is_bridge, [dmin_zn1, dmin_zn2], [elem_zn1, elem_zn2]).
    """
    if isinstance(pose, Pose):
        coords, elements = pose.coords, pose.elements
    else:                                   # array cru (modo legado)
        coords, elements = np.asarray(pose), None

    if coord_elements and elements:
        idx = [i for i, e in enumerate(elements) if e in coord_elements]
        if not idx:                         # nenhum S/N/O -> nunca faz ponte
            return False, [float("inf")] * len(zn_xyz), []
        coords = coords[idx]
        elements = [elements[i] for i in idx]

    mins, hits = [], []
    for z in zn_xyz:
        d = np.linalg.norm(coords - z, axis=1)
        k = int(d.argmin())
        mins.append(round(float(d[k]), 2))
        hits.append(elements[k] if elements else "?")
    return all(m <= cut for m in mins), mins, hits


# ====
#  5b. EXPORTA POSE P/ DISCOVERY STUDIO (Biovia) — PDB legivel
# ====
def _write_model_pdbqt(src_pdbqt, model_idx, dst_pdbqt):
    """Extrai o MODEL model_idx (1-based) de um PDBQT multi-modelo do Vina."""
    grab, idx, cur, out = False, 0, [], []
    for ln in open(src_pdbqt):
        if ln.startswith("MODEL"):
            idx += 1; grab = (idx == model_idx); cur = [ln] if grab else []
        elif ln.startswith("ENDMDL"):
            if grab:
                cur.append(ln); out = cur; break
        elif grab:
            cur.append(ln)
    if not out:                       # fallback: PDBQT sem MODEL (pose unica)
        out = list(open(src_pdbqt))
    open(dst_pdbqt, "w").writelines(out)


def export_pose_for_viewer(out_pdbqt, model_idx, lname, rec_clean, obabel="obabel"):
    """Gera pose_<lig>_best.pdb (melhor modo EM PONTE) e complex_<lig>.pdb
    (receptor limpo + Zn + pose), prontos para arrastar no Discovery Studio."""
    safe = lname.replace(" ", "_")
    pose_pdbqt  = f"pose_{safe}_best.pdbqt"
    pose_pdb    = f"pose_{safe}_best.pdb"
    complex_pdb = f"complex_{safe}.pdb"
    _write_model_pdbqt(out_pdbqt, model_idx, pose_pdbqt)
    r = subprocess.run([obabel, pose_pdbqt, "-O", pose_pdb],
                       capture_output=True, text=True)
    if not os.path.exists(pose_pdb) or os.path.getsize(pose_pdb) == 0:
        raise RuntimeError(f"obabel pose->pdb falhou: {r.stderr[:200]}")
    # complexo = receptor (proteina + Zn) + ligante como HETATM (resname LIG, cadeia Z)
    with open(complex_pdb, "w") as o:
        o.write("REMARK  Complexo receptor + melhor pose EM PONTE (bridged)\n")
        o.write(f"REMARK  ligante: {lname}  |  modelo Vina #{model_idx}\n")
        for ln in open(rec_clean):
            if ln.startswith(("ATOM", "HETATM")):
                o.write(ln)
        o.write("TER\n")
        for ln in open(pose_pdb):
            if ln.startswith(("ATOM", "HETATM")):
                # reescreve como HETATM, resName=LIG, chain=Z, resSeq=900
                o.write("HETATM" + ln[6:17] + "LIG" + " " + "Z" + " 900" + ln[26:])
        o.write("END\n")
    return pose_pdb, complex_pdb


# ====
#  6. PIPELINE POR ALVO
# ====
def dock_target(target, tconf, ligands, tools, workdir, bridge_cut=DEFAULT_BRIDGE_CUT,
                npts=40, ga_runs=20, engine="vina", exhaustiveness=16, dry_run=False,
                coord_elements=COORDINATING_ELEMENTS):
    os.makedirs(workdir, exist_ok=True); cwd0 = os.getcwd(); os.chdir(workdir)
    rows = []
    try:
        log.info("=== %s (%s) ===", target, tconf["pdb_file"])
        log.info("  Bridge filter: cutoff %.1f A | elementos coordenantes: %s",
                 bridge_cut, ",".join(sorted(coord_elements)) if coord_elements else "TODOS")
        rec_clean = "receptor_clean.pdb"
        zn_all = clean_receptor(tconf["pdb_path"], tconf["chain"],
                    tconf["keep_zn"], tconf["remove_het"], rec_clean)
        zn_xyz = select_catalytic_zn_pair(zn_all)
        center = np.mean(zn_xyz, axis=0)
        rec_q = "receptor.pdbqt"
        prepare_receptor(rec_clean, rec_q,
                          python2=tools["pythonsh"],
                          prepare_receptor4=tools["prepare_receptor4"],
                          mgl_pckgs=tools["mgl_pckgs"])
        rec_zn = "receptor_zn.pdbqt"
        add_zinc_pseudo(rec_q, rec_zn, tools["zinc_pseudo"], tools["pythonsh"])

        for lname, smi in ligands:
            lig_q = f"lig_{lname}.pdbqt".replace(" ", "_")
            if not prepare_ligand(smi, lig_q, tools["obabel"]):
                rows.append(dict(target=target, ligand=lname, best_dG=None,
                    bridged_best_dG=None, note="prep_ligante_falhou")); continue
            gpf = "grid.gpf"
            try:
                write_gpf(rec_zn, lig_q, center, gpf, npts=npts, ad4zn_dat=tools["ad4zn_dat"])
                if dry_run:
                    if engine == "vina":
                        size = npts * 0.375
                        cmd = (f"vina --scoring ad4 --maps receptor --ligand {lig_q} "
                               f"--center_x {center[0]:.3f} --center_y {center[1]:.3f} "
                               f"--center_z {center[2]:.3f} --size_x {size:.1f} "
                               f"--size_y {size:.1f} --size_z {size:.1f} "
                               f"--exhaustiveness {exhaustiveness} --out out_{lname}.pdbqt")
                    else:
                        cmd = f"autodock4 -p dock_{lname}.dpf -l dock_{lname}.dlg"
                    log.info("  [DRY-RUN] %s/%s OK. Comando docking seria:\n    %s",
                             target, lname, cmd)
                    rows.append(dict(target=target, ligand=lname, best_dG=None,
                        bridged_best_dG=None, note="dry_run")); continue
                run_autogrid(gpf, tools["autogrid"])
                if engine == "vina":
                    out_q = f"out_{lname}.pdbqt".replace(" ", "_")
                    out_log = f"vina_{lname}.log".replace(" ", "_")
                    run_vina_ad4(lig_q, center, out_q, out_log, vina=tools["vina"],
                                 npts=npts, exhaustiveness=exhaustiveness)
                    models = parse_vina_pdbqt_models(out_q)
                else:
                    dpf = f"dock_{lname}.dpf".replace(" ", "_")
                    dlg = f"dock_{lname}.dlg".replace(" ", "_")
                    write_dpf(lig_q, rec_zn, dpf, ga_runs=ga_runs)
                    run_autodock(dpf, dlg, tools["autodock"])
                    models = parse_dlg_models(dlg)
            except Exception as e:
                log.error("  [%s/%s] docking falhou: %s", target, lname, e)
                rows.append(dict(target=target, ligand=lname, best_dG=None,
                    bridged_best_dG=None, note=f"docking_erro:{e}")); continue
            if not models:
                rows.append(dict(target=target, ligand=lname, best_dG=None,
                    bridged_best_dG=None, note="sem_poses")); continue
            energies = [m.energy for m in models if m.energy is not None]
            best = min(energies) if energies else None
            # ligante sem nenhum heteroatomo coordenante (alcano puro): ponte impossivel
            has_coord = (not coord_elements) or any(
                e in coord_elements for m in models for e in m.elements)
            # avalia ponte por modelo, guardando o indice (1-based = numero do MODEL)
            valid = []  # (energy, dmins, hits, model_idx)
            for i, m in enumerate(models, start=1):
                if m.energy is None:
                    continue
                ok, dmin, hits = bridges_both_zn(m, zn_xyz, bridge_cut, coord_elements)
                if ok:
                    valid.append((m.energy, dmin, hits, i))
            bb = min(valid, key=lambda x: x[0]) if valid else (None, None, None, None)
            note = "ok" if valid else ("sem_heteroatomo_coordenante" if not has_coord
                                       else "zero_pontes")
            rows.append(dict(target=target, ligand=lname, best_dG=best,
                    bridged_best_dG=bb[0], zn_min_dists=bb[1],
                    bridge_atoms=None if bb[2] is None else "/".join(bb[2]),
                    n_bridged=len(valid), n_poses=len(models), note=note))
            log.info("  %-24s best=%s  bridged_best=%s  (%d/%d em ponte)%s",
                    lname, "n/a" if best is None else f"{best:.2f}",
                    "n/a" if bb[0] is None else f"{bb[0]:.2f}", len(valid), len(models),
                    "  [sem S/N/O — rejeitado]" if not has_coord else "")
            # ---- exporta melhor pose EM PONTE p/ Discovery Studio (so engine vina) ----
            if engine == "vina" and bb[3] is not None:
                try:
                    p_pdb, c_pdb = export_pose_for_viewer(
                        out_q, bb[3], lname, rec_clean, obabel=tools["obabel"])
                    log.info("  Biovia: %s + %s (MODEL #%d, dG=%.2f, ponte via %s)",
                            p_pdb, c_pdb, bb[3], bb[0], "/".join(bb[2]))
                except Exception as e:
                    log.warning("  Falha ao exportar pose p/ viewer (%s): %s", lname, e)
            elif engine == "vina":
                log.info("  Sem pose em ponte p/ %s — nada exportado ao viewer.", lname)
    finally:
        os.chdir(cwd0)
    return rows


# ====
#  SAIDAS
# ====
RESULT_COLS = ["target","ligand","best_dG","bridged_best_dG","zn_min_dists",
               "bridge_atoms","n_bridged","n_poses","note"]


def write_results(all_rows, out_csv="mbl_results.csv", cols=None):
    cols = cols or RESULT_COLS
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in all_rows: w.writerow(r)
    log.info("Resultados -> %s", out_csv)


def write_comparison(all_rows, control, out_csv="comparison_vs_control.csv"):
    by_t = {}
    for r in all_rows:
        by_t.setdefault(r["target"], {})[r["ligand"]] = r.get("bridged_best_dG")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["target","ligand","bridged_dG","control","control_dG",
                    "delta_vs_control","melhor_que_controle"])
        for t, ligs in by_t.items():
            ctrl = ligs.get(control)
            for lname, dG in ligs.items():
                if lname == control: continue
                delta = (dG - ctrl) if (dG is not None and ctrl is not None) else None
                w.writerow([t, lname, dG, control, ctrl,
                    None if delta is None else round(delta, 2),
                    delta is not None and delta < 0])
    log.info("Comparacao vs controle -> %s", out_csv)


# ====
#  7. REANALISE — re-aplica o bridge filter sobre out_*.pdbqt existentes
#     (incorpora o antigo reanalyze_bridges.py; NAO re-executa docking)
# ====
def read_zn_from_receptor(path):
    """Le os Zn de receptor_zn.pdbqt (ou receptor_clean.pdb).
    Ignora os pseudo-atomos TZ. Retorna [{'resnum':int,'xyz':array}, ...]."""
    zns = []
    for ln in open(path):
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        el = element_of(ln)
        if el == "TZ":
            continue
        if el != "ZN" and ln[17:20].strip().upper() != "ZN":
            continue
        try:
            xyz = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
        except ValueError:
            continue
        rs = ln[22:26].strip()
        zns.append({"resnum": int(rs) if rs.lstrip("-").isdigit() else 0, "xyz": xyz})
    return zns


def get_zn_coords(target_dir, keep_zn=None):
    """Par de Zn cataliticos do alvo, lido do receptor ja preparado.
    keep_zn: lista/str de numeros de residuo a preservar (opcional)."""
    target_dir = Path(target_dir)
    src = None
    for cand in ("receptor_zn.pdbqt", "receptor.pdbqt", "receptor_clean.pdb"):
        if (target_dir / cand).exists():
            src = target_dir / cand; break
    if src is None:
        return None
    zns = read_zn_from_receptor(src)
    if len(zns) < 2:
        return None
    if keep_zn:
        if isinstance(keep_zn, str):
            keep_zn = [x for x in re.split(r"[,\s]+", keep_zn) if x]
        want = {int(x) for x in keep_zn if str(x).lstrip("-").isdigit()}
        sel = [a for a in zns if a["resnum"] in want]
        if len(sel) == 2:
            return [sel[0]["xyz"], sel[1]["xyz"]]
    return select_catalytic_zn_pair([a["xyz"] for a in zns])


def find_out_pdbqt(target_dir, ligand_name):
    """Localiza out_<LIGANTE>.pdbqt tolerando variacoes de nome."""
    target_dir = Path(target_dir)
    variants = [ligand_name,
                ligand_name.replace(" ", "_"),
                ligand_name.replace(" ", ""),
                ligand_name.replace(",", "_"),
                ligand_name.replace(" ", "_").replace(",", "_")]
    for v in variants:
        c = target_dir / f"out_{v}.pdbqt"
        if c.exists():
            return c
    stem = re.split(r"[\s,]+", ligand_name.strip())[0].lower()
    if stem:
        for f in sorted(target_dir.glob("out_*.pdbqt")):
            if stem in f.stem.lower():
                return f
    return None


def load_keep_zn_map(path):
    """{target: keep_zn} a partir do CSV de alvos (colunas keep_zn ou zn_keep)."""
    if not path or not os.path.exists(path):
        return {}
    out = {}
    for r in csv.DictReader(open(path)):
        name = (r.get("target") or r.get("name") or "").strip()
        keep = (r.get("keep_zn") or r.get("zn_keep") or "").strip()
        if name:
            out[name] = keep
    return out


def reanalyze(outdir, original_csv=None, bridge_cut=DEFAULT_BRIDGE_CUT,
              coord_elements=COORDINATING_ELEMENTS, keep_zn_map=None,
              export_poses=False, obabel="obabel"):
    """Refaz o bridge filter sobre os PDBQT ja produzidos. Retorna rows."""
    outdir = Path(outdir)
    keep_zn_map = keep_zn_map or {}
    orig = {}
    if original_csv and os.path.exists(original_csv):
        for r in csv.DictReader(open(original_csv)):
            if not r.get("target") or r["target"].strip() == "target":
                continue
            orig.setdefault(r["target"].strip(), {})[r["ligand"].strip()] = r
        log.info("Resultados originais lidos: %s", original_csv)
    else:
        log.info("Sem CSV original — ligantes serao descobertos por out_*.pdbqt.")

    log.info("Reanalise: cutoff %.1f A | elementos coordenantes: %s",
             bridge_cut, ",".join(sorted(coord_elements)) if coord_elements else "TODOS")

    rows = []
    tdirs = sorted([d for d in outdir.iterdir() if d.is_dir()]) if outdir.is_dir() else []
    if not tdirs:
        raise RuntimeError(f"Nenhuma subpasta de alvo em {outdir}")
    for tdir in tdirs:
        target = tdir.name
        try:
            zn_xyz = get_zn_coords(tdir, keep_zn_map.get(target))
        except Exception as e:
            log.error("[%s] Zn nao resolvido: %s", target, e); zn_xyz = None
        if not zn_xyz:
            log.warning("[%s] receptor com Zn nao encontrado — alvo pulado.", target)
            for lname in sorted(orig.get(target, {})):
                rows.append(dict(target=target, ligand=lname, note="zn_nao_encontrado"))
            continue
        d = float(np.linalg.norm(zn_xyz[0] - zn_xyz[1]))
        log.info("=== %s === Zn1-Zn2 = %.2f A", target, d)

        ligs = sorted(orig.get(target, {})) or \
               sorted(f.stem[4:] for f in tdir.glob("out_*.pdbqt"))
        if not ligs:
            log.warning("[%s] nenhum out_*.pdbqt encontrado.", target); continue

        for lname in ligs:
            o = orig.get(target, {}).get(lname, {})
            path = find_out_pdbqt(tdir, lname)
            if path is None:
                log.warning("  %-30s out_*.pdbqt nao localizado", lname)
                rows.append(dict(target=target, ligand=lname,
                                 best_dG=o.get("best_dG"), note="pdbqt_nao_localizado",
                                 bridged_best_dG_ORIGINAL=o.get("bridged_best_dG"),
                                 n_bridged_ORIGINAL=o.get("n_bridged")))
                continue
            models = parse_vina_pdbqt_models(str(path))
            if not models:
                log.warning("  %-30s sem poses em %s", lname, path.name)
                rows.append(dict(target=target, ligand=lname, note="sem_poses",
                                 bridged_best_dG_ORIGINAL=o.get("bridged_best_dG"),
                                 n_bridged_ORIGINAL=o.get("n_bridged")))
                continue
            energies = [m.energy for m in models if m.energy is not None]
            best = min(energies) if energies else None
            if best is None and o.get("best_dG"):
                try: best = float(o["best_dG"])
                except ValueError: pass
            has_coord = (not coord_elements) or any(
                e in coord_elements for m in models for e in m.elements)
            valid = []
            for i, m in enumerate(models, start=1):
                ok, dmin, hits = bridges_both_zn(m, zn_xyz, bridge_cut, coord_elements)
                if ok and m.energy is not None:
                    valid.append((m.energy, dmin, hits, i))
            bb = min(valid, key=lambda x: x[0]) if valid else (None, None, None, None)
            note = "ok" if valid else ("sem_heteroatomo_coordenante" if not has_coord
                                       else "zero_pontes")
            rows.append(dict(target=target, ligand=lname, best_dG=best,
                bridged_best_dG=bb[0], zn_min_dists=bb[1],
                bridge_atoms=None if bb[2] is None else "/".join(bb[2]),
                n_bridged=len(valid), n_poses=len(models), note=note,
                bridged_best_dG_ORIGINAL=o.get("bridged_best_dG"),
                n_bridged_ORIGINAL=o.get("n_bridged")))
            log.info("  %s %-30s %d/%d em ponte  bridged_best=%s%s",
                     "OK " if valid else "REJ", lname, len(valid), len(models),
                     "n/a" if bb[0] is None else f"{bb[0]:.2f}",
                     "  [sem S/N/O]" if not has_coord else "")
            if export_poses and bb[3] is not None:
                rec_clean = tdir / "receptor_clean.pdb"
                if not rec_clean.exists():
                    log.warning("  receptor_clean.pdb ausente — pose nao exportada.")
                    continue
                cwd0 = os.getcwd(); os.chdir(tdir)
                try:
                    p_pdb, c_pdb = export_pose_for_viewer(
                        path.name, bb[3], lname, "receptor_clean.pdb", obabel=obabel)
                    log.info("  Biovia: %s + %s (MODEL #%d)", p_pdb, c_pdb, bb[3])
                except Exception as e:
                    log.warning("  Falha ao exportar pose (%s): %s", lname, e)
                finally:
                    os.chdir(cwd0)
    return rows


# ====
#  8. REDOCKING DO LIGANTE CRISTALOGRAFICO (controle positivo do protocolo)
#     Extrai o HETATM co-cristalizado, redocka no proprio receptor e mede o
#     RMSD contra a pose depositada. Criterio aceito: RMSD <= 2.0 A.
# ====
REDOCK_RMSD_PASS = 2.0        # A — limiar convencional de sucesso


def extract_hetatm_ligand(pdb_in, resname, chain=None, resnum=None, out_pdb=None):
    """Escreve as linhas HETATM do ligante co-cristalizado num PDB isolado.
    Mantem apenas altloc ' '/'A'. Retorna (caminho, n_atomos_pesados)."""
    resname = resname.strip().upper()
    lines, heavy, seen_alt = [], 0, set()
    for ln in open(pdb_in):
        if not ln.startswith("HETATM"):
            continue
        if ln[17:20].strip().upper() != resname:
            continue
        if chain and ln[21] != chain:
            continue
        if resnum is not None and ln[22:26].strip() != str(resnum):
            continue
        alt = ln[16]
        if alt not in (" ", "A"):
            seen_alt.add(alt); continue
        lines.append(ln)
        if element_of(ln) != "H":
            heavy += 1
    if not lines:
        raise RuntimeError(f"Ligante '{resname}' nao encontrado em {pdb_in}"
                           + (f" (cadeia {chain})" if chain else ""))
    if seen_alt:
        log.info("    altloc %s descartado(s); mantido A/' '.", ",".join(sorted(seen_alt)))
    out_pdb = out_pdb or f"xtal_{resname}.pdb"
    with open(out_pdb, "w") as fh:
        fh.writelines(lines); fh.write("END\n")
    return out_pdb, heavy


def load_pdb_heavy_coords(path):
    """Coordenadas Nx3 dos atomos pesados + lista de elementos (PDB ou PDBQT)."""
    xyz, els = [], []
    for ln in open(path):
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        e = element_of(ln)
        if e in ("H", "?"):
            continue
        xyz.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]); els.append(e)
    return np.array(xyz), els


def symmetry_corrected_rmsd(ref_pdb, pose_pdb, obrms="obrms"):
    """RMSD em atomos pesados, sem sobreposicao (in place), corrigido para
    simetria topologica via `obrms` do Open Babel. Se obrms nao existir, usa
    o fallback ingenuo por ordem de atomo, que SUPERESTIMA o RMSD em ligantes
    com grupos equivalentes (fenilas, carboxilatos) — sinalizado no retorno.

    Retorna (rmsd, metodo)."""
    if shutil.which(obrms):
        try:
            r = subprocess.run([obrms, ref_pdb, pose_pdb],
                               capture_output=True, text=True, timeout=120)
            m = re.search(r"([0-9]*\.?[0-9]+)\s*$", r.stdout.strip().splitlines()[-1]) \
                if r.stdout.strip() else None
            if r.returncode == 0 and m:
                return float(m.group(1)), "obrms"
            log.warning("    obrms nao retornou valor (%s); usando fallback.",
                        (r.stderr or r.stdout)[:120])
        except Exception as e:
            log.warning("    obrms falhou (%s); usando fallback.", e)
    a, ea = load_pdb_heavy_coords(ref_pdb)
    b, eb = load_pdb_heavy_coords(pose_pdb)
    if len(a) != len(b):
        raise RuntimeError(f"RMSD: contagem de atomos pesados difere "
                           f"(referencia {len(a)}, pose {len(b)}). "
                           f"Instale o Open Babel (obrms) ou verifique a protonacao.")
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean())), "naive_order"


def redock_target(target, tconf, tools, workdir, bridge_cut=DEFAULT_BRIDGE_CUT,
                  npts=40, exhaustiveness=16, coord_elements=COORDINATING_ELEMENTS,
                  resname=None, dry_run=False):
    """Redocka o ligante co-cristalizado do alvo no seu proprio receptor.

    O ligante e' extraido do PDB original DEPOIS de o receptor ter sido limpo
    (o mesmo remove_het que o retira do receptor identifica-o aqui), convertido
    a PDBQT com obabel preservando as coordenadas cristalograficas como
    referencia, e redockado a partir de conformacao re-gerada.
    """
    os.makedirs(workdir, exist_ok=True); cwd0 = os.getcwd(); os.chdir(workdir)
    try:
        # ligante alvo: explicito, ou o primeiro remove_het que nao seja solvente
        # HETATM que NAO sao ligantes: solventes, tampoes, crioprotetores, ions e
        # placeholders. UNX/UNK/UNL sao "atomo/ligante desconhecido" — densidade
        # nao interpretada, sem topologia: redockar isso nao significa nada.
        SOLV = {"HOH", "WAT", "GOL", "EDO", "PEG", "PG4", "1PE", "SO4", "PO4", "CL",
                "NA", "K", "MG", "CA", "ZN", "NI", "CD", "ACT", "DMS", "TRS", "MPD",
                "BCN", "IMD", "FMT", "NO3", "ACY", "EPE", "MES", "CIT", "TAR", "BME",
                "IOD", "BR", "FLC", "UNX", "UNK", "UNL", "PEO", "OXY"}
        MIN_HEAVY = 6          # abaixo disso o RMSD e' pouco informativo
        cands = [h.strip().upper() for h in tconf.get("remove_het", []) if h.strip()]
        lig_res = (resname or "").strip().upper() or next(
            (h for h in cands if h not in SOLV), None)
        if not lig_res:
            log.info("[%s] sem ligante co-cristalizado em remove_het — redocking pulado.", target)
            return None
        log.info("=== REDOCK %s (%s, ligante %s) ===", target, tconf["pdb_file"], lig_res)

        rec_clean = "receptor_clean.pdb"
        zn_all = clean_receptor(tconf["pdb_path"], tconf["chain"],
                                tconf["keep_zn"], tconf["remove_het"], rec_clean)
        zn_xyz = select_catalytic_zn_pair(zn_all)
        center = np.mean(zn_xyz, axis=0)

        xtal_pdb, n_heavy = extract_hetatm_ligand(
            tconf["pdb_path"], lig_res, chain=tconf["chain"], out_pdb=f"xtal_{lig_res}.pdb")
        log.info("  Ligante cristalografico: %d atomos pesados -> %s", n_heavy, xtal_pdb)
        if n_heavy < MIN_HEAVY:
            log.warning("  [%s] '%s' tem apenas %d atomos pesados (<%d): fragmento ou "
                        "ion, RMSD pouco informativo — redocking pulado. Use "
                        "--redock_resname para forcar outro HETATM.",
                        target, lig_res, n_heavy, MIN_HEAVY)
            return dict(target=target, pdb=tconf["pdb_file"], ligand=lig_res,
                        heavy_atoms=n_heavy, note="ligante_pequeno_demais")

        # distancia minima heteroatomo->Zn na pose CRISTALOGRAFICA
        xr, xe = load_pdb_heavy_coords(xtal_pdb)
        idx = [i for i, e in enumerate(xe) if (not coord_elements) or e in coord_elements]
        xtal_bridge, xtal_d, xtal_hits = (False, [None, None], [])
        if idx:
            xtal_bridge, xtal_d, xtal_hits = bridges_both_zn(
                Pose(None, xr[idx], [xe[i] for i in idx]), zn_xyz, bridge_cut, coord_elements)
        log.info("  Pose cristalografica: ponte=%s dist=%s atomos=%s",
                 xtal_bridge, xtal_d, "/".join(xtal_hits) if xtal_hits else "-")

        # PDBQT do ligante a partir das coordenadas cristalograficas
        lig_q = f"lig_redock_{lig_res}.pdbqt"
        r = subprocess.run([tools["obabel"], xtal_pdb, "-O", lig_q, "-h",
                            "--partialcharge", "gasteiger"],
                           capture_output=True, text=True, timeout=180)
        if r.returncode != 0 or not os.path.exists(lig_q) or os.path.getsize(lig_q) == 0:
            raise RuntimeError(f"obabel falhou ao converter {xtal_pdb}: {r.stderr[:200]}")

        if dry_run:
            log.info("  [DRY-RUN] extracao e preparo OK; docking nao executado.")
            return dict(target=target, pdb=tconf["pdb_file"], ligand=lig_res,
                        heavy_atoms=n_heavy, note="dry_run")

        rec_q = "receptor.pdbqt"
        prepare_receptor(rec_clean, rec_q, python2=tools["pythonsh"],
                         prepare_receptor4=tools["prepare_receptor4"],
                         mgl_pckgs=tools["mgl_pckgs"])
        rec_zn = "receptor_zn.pdbqt"
        add_zinc_pseudo(rec_q, rec_zn, tools["zinc_pseudo"], tools["pythonsh"])
        write_gpf(rec_zn, lig_q, center, "grid.gpf", npts=npts, ad4zn_dat=tools["ad4zn_dat"])
        run_autogrid("grid.gpf", tools["autogrid"])

        out_q = f"redock_{lig_res}.pdbqt"
        run_vina_ad4(lig_q, center, out_q, f"redock_{lig_res}.log", vina=tools["vina"],
                     npts=npts, exhaustiveness=exhaustiveness)
        models = parse_vina_pdbqt_models(out_q)
        if not models:
            raise RuntimeError(f"redocking sem poses para {lig_res}")

        # RMSD de cada pose contra a referencia cristalografica
        rmsds, method = [], None
        for i in range(1, len(models) + 1):
            _write_model_pdbqt(out_q, i, f"_m{i}.pdbqt")
            subprocess.run([tools["obabel"], f"_m{i}.pdbqt", "-O", f"_m{i}.pdb"],
                           capture_output=True, text=True, timeout=120)
            try:
                v, method = symmetry_corrected_rmsd(xtal_pdb, f"_m{i}.pdb",
                                                    obrms=tools.get("obrms", "obrms"))
            except Exception as e:
                log.warning("    RMSD do modelo %d falhou: %s", i, e); v = None
            rmsds.append(v)
        valid = [(v, i) for i, v in enumerate(rmsds, 1) if v is not None]
        if not valid:
            raise RuntimeError("nenhum RMSD pode ser calculado")
        top_rmsd = rmsds[0]                       # pose de melhor escore
        best_rmsd, best_i = min(valid)            # pose mais proxima da referencia
        e_top = models[0].energy
        ok, dmin, hits = bridges_both_zn(models[0], zn_xyz, bridge_cut, coord_elements)
        log.info("  RMSD top1 = %s A | melhor RMSD = %.2f A (pose #%d) | metodo=%s",
                 "n/a" if top_rmsd is None else f"{top_rmsd:.2f}", best_rmsd, best_i, method)
        log.info("  %s (criterio: RMSD top1 <= %.1f A)",
                 "APROVADO" if (top_rmsd is not None and top_rmsd <= REDOCK_RMSD_PASS)
                 else "REPROVADO", REDOCK_RMSD_PASS)
        if method == "naive_order":
            log.warning("  RMSD sem correcao de simetria (obrms ausente) — valor e' um "
                        "LIMITE SUPERIOR; instale o Open Babel para o valor correto.")
        return dict(target=target, pdb=tconf["pdb_file"], ligand=lig_res,
                    heavy_atoms=n_heavy, top_pose_dG=e_top,
                    rmsd_top_pose=None if top_rmsd is None else round(top_rmsd, 2),
                    rmsd_best_pose=round(best_rmsd, 2), best_pose_index=best_i,
                    rmsd_method=method, n_poses=len(models),
                    xtal_bridges=xtal_bridge, xtal_zn_dists=xtal_d,
                    xtal_bridge_atoms="/".join(xtal_hits) if xtal_hits else None,
                    redock_bridges=ok, redock_zn_dists=dmin,
                    passes=bool(top_rmsd is not None and top_rmsd <= REDOCK_RMSD_PASS))
    finally:
        os.chdir(cwd0)


REDOCK_COLS = ["target", "pdb", "ligand", "heavy_atoms", "top_pose_dG", "rmsd_top_pose",
               "rmsd_best_pose", "best_pose_index", "rmsd_method", "n_poses",
               "xtal_bridges", "xtal_zn_dists", "xtal_bridge_atoms",
               "redock_bridges", "redock_zn_dists", "passes", "note"]


# ====
#  CLI
# ====
def load_targets_csv(path, pdb_dir):
    targets = {}
    pdb_dir = os.path.abspath(pdb_dir)        # o pipeline faz chdir p/ o workdir:
                                              # caminhos relativos quebrariam depois
    for r in csv.DictReader(open(path)):
        if not r.get("center_x") or r.get("status","").startswith("PENDENTE"): continue
        code = r["pdb"]; pdb_path = None
        for fn in os.listdir(pdb_dir):
            if fn.lower().endswith(".pdb") and code.lower() in fn.lower():
                pdb_path = os.path.join(pdb_dir, fn); break
        if pdb_path is None:
            log.warning("PDB '%s' nao encontrado em %s — pulando %s.", code, pdb_dir, r["target"]); continue
        targets[r["target"]] = dict(pdb_file=code, pdb_path=pdb_path, chain=r["chain"],
            keep_zn=[x for x in r["keep_zn"].split(",") if x],
            remove_het=[x for x in r.get("remove_het","").split(",") if x])
    return targets


def load_ligands_csv(path):
    return [(r["name"].strip(), r["smiles"].strip()) for r in csv.DictReader(open(path))]


def _auto_mgl_pckgs(prepare_receptor4):
    """Localiza o diretorio MGLToolsPckgs a partir do caminho do prepare_receptor4.py.

    Dois layouts sao comuns e ambos precisam funcionar:
      (a) .../MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py
      (b) .../envs/<env>/bin/prepare_receptor4.py   (conda; wrapper em bin/)
    No caso (b) o MGLToolsPckgs e' irmao de bin/, entao subimos o caminho
    procurando um diretorio que contenha MGLToolsPckgs."""
    if not prepare_receptor4:
        return None
    p = os.path.abspath(prepare_receptor4)
    # (a) algum ancestral chamado MGLToolsPckgs
    q = p
    while True:
        if os.path.basename(q) == "MGLToolsPckgs" and os.path.isdir(q):
            return q
        parent = os.path.dirname(q)
        if parent == q:
            break
        q = parent
    # (b) algum ancestral que CONTENHA um MGLToolsPckgs (ex.: prefixo do conda)
    q = os.path.dirname(p)
    while True:
        cand = os.path.join(q, "MGLToolsPckgs")
        if os.path.isdir(cand):
            return cand
        parent = os.path.dirname(q)
        if parent == q:
            return None
        q = parent


def _auto_zinc_pseudo(zinc_pseudo, mgl_pckgs, prepare_receptor4):
    """Resolve o caminho de zinc_pseudo.py. Procura, nesta ordem: o valor dado
    (se existir como arquivo ou no PATH), a pasta Utilities24 do MGLToolsPckgs,
    e a pasta do proprio prepare_receptor4.py."""
    if zinc_pseudo and (os.path.exists(zinc_pseudo) or shutil.which(zinc_pseudo)):
        return zinc_pseudo
    name = os.path.basename(zinc_pseudo or "zinc_pseudo.py")
    cands = [os.path.abspath(name),                                  # diretorio atual
             os.path.join(os.path.dirname(os.path.abspath(__file__)), name)]  # junto ao script
    if mgl_pckgs:
        cands += [os.path.join(mgl_pckgs, "AutoDockTools", "Utilities24", name),
                  os.path.join(mgl_pckgs, "AutoDockTools", name)]
    if prepare_receptor4:
        d = os.path.dirname(os.path.abspath(prepare_receptor4))
        cands += [os.path.join(d, name),
                  os.path.join(os.path.dirname(d), "bin", name)]
    for c in cands:
        if os.path.exists(c):
            return c
    return zinc_pseudo
    return p if os.path.isdir(p) else None


def main():
    ap = argparse.ArgumentParser(
        description="Docking di-nuclear MBL (AD4Zn + bridge filter S/N/O).",
        epilog="Modo reanalise: mbl_dock.py --reanalyze --outdir mbl_out "
               "--targets mbl_targets.csv")
    ap.add_argument("--targets", default=None,
                    help="CSV de alvos (obrigatorio no docking; opcional em --reanalyze).")
    ap.add_argument("--ligands", default=None,
                    help="CSV de ligantes (obrigatorio no docking).")
    ap.add_argument("--pdb_dir", default=".")
    ap.add_argument("--outdir", default="mbl_out")
    ap.add_argument("--control", default="Chlorhexidine")
    ap.add_argument("--bridge_cut", type=float, default=DEFAULT_BRIDGE_CUT,
                    help=f"Cutoff heteroatomo->Zn em A (default: {DEFAULT_BRIDGE_CUT}).")
    ap.add_argument("--coord_elements", default="S,N,O",
                    help="Elementos aceitos como coordenantes do Zn no bridge filter "
                         "(default: S,N,O). Use 'any' para o comportamento antigo, "
                         "geometrico, que aceitava carbono.")
    ap.add_argument("--reanalyze", action="store_true",
                    help="NAO roda docking: re-aplica o bridge filter nos out_*.pdbqt "
                         "ja existentes em --outdir e grava --reanalyze_out.")
    ap.add_argument("--reanalyze_out", default=None,
                    help="CSV de saida da reanalise (default: <outdir>/mbl_results_fixed.csv).")
    ap.add_argument("--original_csv", default=None,
                    help="CSV de resultados original p/ comparar (default: "
                         "<outdir>/mbl_results.csv).")
    ap.add_argument("--export_poses", action="store_true",
                    help="Na reanalise, reexporta a melhor pose em ponte p/ Discovery Studio.")
    ap.add_argument("--redock", action="store_true",
                    help="CONTROLE POSITIVO: redocka o ligante co-cristalizado de cada alvo "
                         "no proprio receptor e mede o RMSD contra a pose depositada. "
                         "Nao docka os ligantes de --ligands.")
    ap.add_argument("--redock_resname", default=None,
                    help="Codigo do HETATM a redockar (ex.: QZH). Se omitido, usa o primeiro "
                         "item de remove_het que nao seja solvente/tampao.")
    ap.add_argument("--redock_out", default=None,
                    help="CSV de saida do redocking (default: <outdir>/redock_validation.csv).")
    ap.add_argument("--npts", type=int, default=40)
    ap.add_argument("--ga_runs", type=int, default=20)
    ap.add_argument("--engine", choices=["vina","autodock4"], default="vina")
    ap.add_argument("--exhaustiveness", type=int, default=16)
    ap.add_argument("--dry_run", action="store_true",
                    help="Valida binarios, limpa receptor, prepara ligantes e gera GPF; "
                    "mostra o comando de docking — SEM rodar autogrid/docking.")
    ap.add_argument("--obabel", default="obabel")
    ap.add_argument("--obrms", default="obrms",
                    help="Executavel obrms (Open Babel) p/ RMSD com correcao de simetria.")
    ap.add_argument("--autogrid", default="autogrid4")
    ap.add_argument("--autodock", default="autodock4")
    ap.add_argument("--vina", default="vina")
    ap.add_argument("--pythonsh", default="python2.7",
                    help="Interpretador Python 2.7 REAL do MGLTools "
                         "(ex.: .../envs/mbl/bin/python2.7). NAO use o wrapper 'pythonsh'.")
    ap.add_argument("--zinc_pseudo", default="zinc_pseudo.py")
    ap.add_argument("--ad4zn_dat", default="AD4Zn.dat")
    ap.add_argument("--prepare_receptor4", default=None,
                help="Caminho para prepare_receptor4.py (MGLTools). "
                     "Obrigatorio no docking; dispensavel em --reanalyze.")
    ap.add_argument("--mgl_pckgs", default=None,
                help="Diretorio MGLToolsPckgs (para PYTHONPATH). "
                     "Se omitido, deriva do caminho do prepare_receptor4.")
    args = ap.parse_args()
    coord_elements = parse_coord_elements(args.coord_elements)

    # ---------- MODO REANALISE (sem docking) ----------
    if args.reanalyze:
        original = args.original_csv or os.path.join(args.outdir, "mbl_results.csv")
        out_csv = args.reanalyze_out or os.path.join(args.outdir, "mbl_results_fixed.csv")
        keep_map = load_keep_zn_map(args.targets) if args.targets else {}
        rows = reanalyze(args.outdir, original_csv=original,
                         bridge_cut=args.bridge_cut, coord_elements=coord_elements,
                         keep_zn_map=keep_map, export_poses=args.export_poses,
                         obabel=args.obabel)
        cols = RESULT_COLS + ["bridged_best_dG_ORIGINAL", "n_bridged_ORIGINAL"]
        write_results(rows, out_csv, cols=cols)
        n_ok = sum(1 for r in rows if r.get("n_bridged"))
        n_sem = sum(1 for r in rows if r.get("note") == "sem_heteroatomo_coordenante")
        log.info("Reanalise: %d/%d ligantes com pelo menos uma pose em ponte; "
                 "%d rejeitados por nao terem S/N/O.", n_ok, len(rows), n_sem)
        if any(r.get("ligand") == args.control for r in rows):
            write_comparison(rows, args.control,
                             os.path.join(args.outdir, "comparison_vs_control_fixed.csv"))
        log.info("Concluido (reanalise).")
        return

    # ---------- MODO DOCKING / REDOCKING ----------
    req = [("--targets", args.targets), ("--prepare_receptor4", args.prepare_receptor4)]
    if not args.redock:
        req.append(("--ligands", args.ligands))
    faltando = [n for n, v in req if not v]
    if faltando:
        log.error("Argumentos obrigatorios ausentes: %s", ", ".join(faltando))
        sys.exit(1)

    need = [args.obabel]
    if not args.dry_run:
        need += [args.autogrid, args.vina if args.engine == "vina" else args.autodock]
    for exe in need:
        if shutil.which(exe) is None and not os.path.exists(exe):
            log.error("Executavel nao encontrado no PATH: %s", exe); sys.exit(1)
    if not args.dry_run:
        if not os.path.exists(args.pythonsh) and shutil.which(args.pythonsh) is None:
            log.error("python2.7 nao encontrado: %s", args.pythonsh); sys.exit(1)
        if not os.path.exists(args.prepare_receptor4):
            log.error("prepare_receptor4.py nao encontrado: %s", args.prepare_receptor4); sys.exit(1)
    if not _HAS_RDKIT:
        log.warning("RDKit ausente — usando fallback obabel --gen3d.")

    mgl_pckgs = args.mgl_pckgs or _auto_mgl_pckgs(args.prepare_receptor4)
    if mgl_pckgs:
        log.info("MGLToolsPckgs: %s", mgl_pckgs)
    else:
        log.warning("MGLToolsPckgs nao localizado — prepare_receptor4 pode falhar por PYTHONPATH.")

    zinc_pseudo = _auto_zinc_pseudo(args.zinc_pseudo, mgl_pckgs, args.prepare_receptor4)
    if not args.dry_run and not (os.path.exists(zinc_pseudo) or shutil.which(zinc_pseudo)):
        log.error("zinc_pseudo.py nao encontrado (procurado no diretorio atual, ao lado "
                  "deste script, em MGLToolsPckgs/AutoDockTools/Utilities24 e junto ao "
                  "prepare_receptor4.py). Ele faz parte do AutoDock4Zn: baixe-o em "
                  "https://autodock.scripps.edu/resources/autodock4zn/ e passe o caminho "
                  "com --zinc_pseudo. Localize com: find $HOME -name zinc_pseudo.py")
        sys.exit(1)

    # o pipeline faz chdir p/ o workdir: estes caminhos precisam ser absolutos
    ad4zn_dat = args.ad4zn_dat
    if os.path.exists(ad4zn_dat):
        ad4zn_dat = os.path.abspath(ad4zn_dat)
    elif not args.dry_run:
        log.error("AD4Zn.dat nao encontrado em '%s'. Ele acompanha o pacote AutoDock4Zn; "
                  "passe o caminho com --ad4zn_dat.", args.ad4zn_dat)
        sys.exit(1)
    log.info("AD4Zn.dat: %s", ad4zn_dat)
    if os.path.exists(zinc_pseudo):
        zinc_pseudo = os.path.abspath(zinc_pseudo)
    log.info("zinc_pseudo.py: %s", zinc_pseudo)

    tools = dict(obabel=args.obabel, autogrid=args.autogrid, autodock=args.autodock,
                 vina=args.vina, pythonsh=args.pythonsh,
                 zinc_pseudo=zinc_pseudo, ad4zn_dat=ad4zn_dat,
                 prepare_receptor4=args.prepare_receptor4, mgl_pckgs=mgl_pckgs,
                 obrms=args.obrms)
    targets = load_targets_csv(args.targets, args.pdb_dir)

    # ---------- MODO REDOCKING (controle positivo do protocolo) ----------
    if args.redock:
        if shutil.which(args.obrms) is None:
            log.warning("'%s' nao encontrado — o RMSD sera calculado sem correcao de "
                        "simetria e superestimado. Instale o Open Babel.", args.obrms)
        out_csv = args.redock_out or os.path.join(args.outdir, "redock_validation.csv")
        os.makedirs(args.outdir, exist_ok=True)
        rows = []
        for tname, tconf in targets.items():
            try:
                r = redock_target(tname, tconf, tools, os.path.join(args.outdir, tname),
                                  bridge_cut=args.bridge_cut, npts=args.npts,
                                  exhaustiveness=args.exhaustiveness,
                                  coord_elements=coord_elements,
                                  resname=args.redock_resname, dry_run=args.dry_run)
                if r: rows.append(r)
            except Exception as e:
                log.error("[%s] redocking falhou: %s", tname, e)
                rows.append(dict(target=tname, pdb=tconf["pdb_file"], note=f"erro: {e}"))
        if rows:
            write_results(rows, out_csv, cols=REDOCK_COLS)
            done = [r for r in rows if r.get("rmsd_top_pose") is not None]
            if done:
                ok = sum(1 for r in done if r.get("passes"))
                vals = [r["rmsd_top_pose"] for r in done]
                log.info("Redocking: %d/%d alvos com RMSD(top1) <= %.1f A "
                         "(mediana %.2f A, faixa %.2f-%.2f).",
                         ok, len(done), REDOCK_RMSD_PASS,
                         float(np.median(vals)), min(vals), max(vals))
        log.info("Concluido (redocking).")
        return

    ligands = load_ligands_csv(args.ligands)
    log.info("Motor: %s | Alvos: %s", args.engine, list(targets))
    log.info("Ligantes: %s", [l for l, _ in ligands])

    os.makedirs(args.outdir, exist_ok=True); all_rows = []
    for tname, tconf in targets.items():
        wd = os.path.join(args.outdir, tname)
        try:
            all_rows += dock_target(tname, tconf, ligands, tools, wd,
                bridge_cut=args.bridge_cut, npts=args.npts, ga_runs=args.ga_runs,
                engine=args.engine, exhaustiveness=args.exhaustiveness,
                dry_run=args.dry_run, coord_elements=coord_elements)
        except Exception as e:
            log.error("Alvo %s abortado: %s", tname, e)
            all_rows.append(dict(target=tname, ligand="-", best_dG=None,
                    bridged_best_dG=None, note=f"alvo_erro:{e}"))

    write_results(all_rows, os.path.join(args.outdir, "mbl_results.csv"))
    if not args.dry_run and any(r["ligand"] == args.control for r in all_rows):
        write_comparison(all_rows, args.control, os.path.join(args.outdir, "comparison_vs_control.csv"))
    log.info("Concluido.")


if __name__ == "__main__":
    main()
