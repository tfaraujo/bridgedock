# Geração da tabela de alvos — `mbl_targets.csv`

Este documento descreve como a tabela `mbl_targets.csv` usada pelo `mbl_dock.py`
foi produzida para o estudo de ancoragem em metalo-β-lactamases (MBL), o que nela
é medição geométrica e o que é curadoria humana, e quais são as limitações
conhecidas dos scripts aqui incluídos.

O objetivo é reprodutibilidade honesta: quem baixar os PDBs deve conseguir
reobter os números publicados e identificar sem ambiguidade quais campos
dependeram de decisão do autor.

---

## 1. O que a tabela contém

Uma linha por alvo, com os parâmetros que o `mbl_dock.py` consome para preparar o
receptor e centrar a caixa de busca:

| Coluna | Significado |
|---|---|
| `target` | Nome da enzima (NDM-1, VIM-2, SPM-1, …) |
| `pdb` | Código PDB da estrutura escolhida |
| `chain` | Cadeia usada para o docking |
| `keep_zn` | `resSeq` dos dois Zn²⁺ catalíticos preservados |
| `drop_zn` | `resSeq` de Zn²⁺ adicionais (cristalográficos/adventícios) removidos |
| `remove_het` | HETATM removidos do receptor (solventes, tampões, crioprotetores, inibidores co-cristalizados) |
| `center_x/y/z` | Centro da caixa de busca — ponto médio entre os dois Zn²⁺ preservados |
| `znzn_A` | Distância Zn1–Zn2 em Å |
| `status` | Nota de curadoria (motivo da escolha de cadeia, de estrutura, do que foi removido) |

Alvos incluídos: NDM-1 (8B1W), NDM-5 (6MGY), NDM-7 (7AEZ), VIM-1 (5N5G),
VIM-2 (7A5Z), VIM-7 (2Y87), IMP-1 (7YH9), SPM-1 (5NDB). Todos MBL de subclasse B1,
dinucleares.

---

## 2. Como a tabela foi realmente produzida

A tabela foi montada alvo a alvo com scripts de inspeção estrutural escritos ad hoc:
para cada PDB, listavam-se os Zn²⁺ presentes, as distâncias entre eles, os átomos
doadores (N, O, S) a ≤ 2,9 Å de cada íon e os HETATM do arquivo. A partir desse
levantamento, o autor decidia a cadeia, o par catalítico, os Zn e HET a remover, e
registrava o motivo na coluna `status`.

**Aviso de proveniência.** O conjunto de scripts que gerou a tabela final não foi
preservado. O que sobreviveu está aqui:

- **`gerador_5NDB_SPM-1.py`** — fragmento remanescente daquele fluxo, na forma em
  que foi aplicado ao SPM-1 (5NDB), com verificações residuais de VIM-1 e VIM-2.
  Imprime a esfera de coordenação de cada Zn e o centro do par. Não escreve CSV e
  tem nomes de arquivo fixos no código. É documental, não de produção.
- **`make_mbl_targets.py`** — **reconstrução posterior**, escrita por engenharia
  reversa da tabela publicada. Automatiza a parte puramente geométrica do
  procedimento original e a aplica de forma uniforme a todos os alvos.

O `make_mbl_targets.py` **não é** o script que gerou o `mbl_targets.csv`, e não o
reproduz integralmente. As diferenças estão na seção 4 e são conhecidas, não
acidentais.

---

## 3. Proveniência por coluna

| Coluna | Origem | Reproduzível por `make_mbl_targets.py` |
|---|---|---|
| `target`, `pdb` | Seleção de literatura | — (mapeamento fixo no script) |
| `keep_zn` | Geometria: par de Zn na mesma cadeia com 2,8 ≤ d ≤ 5,0 Å | **Sim** |
| `center_x/y/z` | Geometria: ponto médio do par | **Sim** |
| `znzn_A` | Geometria: distância do par | **Sim** |
| `chain` | Heurística + revisão humana | Parcial (ver 4.1) |
| `drop_zn` | Curadoria | **Não** — coluna não emitida |
| `remove_het` | Curadoria sobre o inventário do arquivo | **Não** — ver 4.2 |
| `status` | Literatura e decisão do autor | **Não** — por natureza |

---

## 4. Divergências conhecidas da reconstrução

Executando

```bash
python make_mbl_targets.py --pdb_dir ./pdb --out mbl_targets_auto.csv
```

sobre os oito PDBs, obtêm-se as colunas geométricas idênticas às publicadas em
**sete dos oito alvos**. As divergências são estas.

### 4.1 SPM-1 (5NDB) — cadeia

O `make_mbl_targets.py` escolhe a cadeia pela ordem alfabética entre os pares de Zn
válidos, e devolve a cadeia A. A tabela publicada usa a cadeia **B**:

| | Cadeia | Zn | Centro (x, y, z) |
|---|---|---|---|
| Publicado (curado) | B | 401, 402 | −43,106 / 16,023 / 8,747 |
| `make_mbl_targets.py` | A | 401, 402 | −45,690 / 21,040 / −12,390 |

O 5NDB tem duas cópias na unidade assimétrica, com três Zn²⁺ em cada
(A: 401, 402, 407; B: 401, 402, 406). A escolha da cadeia B foi feita por inspeção
do sítio, não por ordenação alfabética. Ordem alfabética é um critério arbitrário e
está sinalizado como tal — a alternativa correta (qualidade real do sítio: número de
doadores N/O/S a ≤ 2,9 Å de cada Zn e B-factor médio da primeira esfera) está
prevista para a v2.

**Para reproduzir a linha publicada, a cadeia do SPM-1 deve ser fixada manualmente.**

### 4.2 `remove_het` — lista fixa em vez de inventário

O `make_mbl_targets.py` escreve, para todos os alvos, a mesma lista de solventes e
tampões comuns (`HOH, GOL, EDO, SO4, PO4, MES, TRS, ACT, DMS, PEG`). A tabela
publicada lista os HETATM efetivamente presentes em cada estrutura, o que inclui os
inibidores co-cristalizados: `OQU` (8B1W), `R8W` (7AEZ), `BCN` (5N5G), `QZH` (7A5Z),
`IT0` (7YH9), `8TW` (5NDB).

**Consequência prática:** usada sem revisão, a saída automática deixa o ligante do
cristal dentro do receptor, e a caixa de busca cai sobre um sítio ocupado. Isto não
é um detalhe cosmético — invalida o docking do alvo afetado.

### 4.3 `drop_zn` — não emitido

Zn²⁺ excedentes ao par catalítico não são identificados pela reconstrução:
VIM-1 (303), VIM-2 (302) e SPM-1 (406/407) foram anotados manualmente. Sem essa
coluna, um terceiro íon permanece no receptor.

### 4.4 Nome de coluna

A reconstrução escreve `zn_zn_dist`; a tabela publicada usa `znzn_A`. Renomear
antes de alimentar o `mbl_dock.py`.

---

## 5. Como usar hoje

1. Baixar os PDBs listados na seção 1 e colocá-los em um diretório.
2. Rodar o `make_mbl_targets.py` para obter as colunas geométricas.
3. Conferir a saída contra o `mbl_targets.csv` publicado.
4. Completar `drop_zn`, `remove_het` e `status`, e corrigir a cadeia do SPM-1,
   conforme a tabela publicada.

Para o conjunto deste estudo, o caminho recomendado é simplesmente **usar o
`mbl_targets.csv` versionado neste repositório**, que é a tabela efetivamente
empregada nos cálculos. O `make_mbl_targets.py` serve para auditar as colunas
geométricas e para estender o procedimento a novos alvos, com revisão manual dos
campos de curadoria.

---

## 6. Trabalho futuro (`scaffold_targets.py`, v2)

Esta versão documenta o procedimento tal como foi executado. Uma versão
generalizada, para ser rodada antes do `mbl_dock.py` sobre alvos arbitrários, está
planejada para depois da publicação:

- **`--fetch 8B1W,6MGY`** — download direto do RCSB, partindo apenas de códigos PDB.
- **`--mononuclear`** — suporte a MBL de subclasse B2 (CphA, Sfh-I), com um único
  Zn²⁺ catalítico, hoje classificadas como falha pelo critério de par dinuclear.
- **Escolha de cadeia por qualidade do sítio** — contagem de doadores N/O/S a
  ≤ 2,9 Å de cada Zn e B-factor médio da primeira esfera, substituindo a ordem
  alfabética.
- **Inventário e classificação de HETATM** — solvente/tampão vs. ligante, por lista
  de referência combinada a número de átomos pesados e distância ao sítio.
- **`drop_zn` automático** — Zn²⁺ da cadeia escolhida fora do par catalítico.
- **`status` como diagnóstico factual** — número de cadeias, Zn excedentes, HET
  volumosos próximos ao sítio, coordenação incompleta. O campo permanece de
  preenchimento humano: anotações como *"5NDB substitui 2FHX oxidado"* vêm da
  literatura, não do arquivo, e nenhuma heurística deve fabricá-las.

---

## 7. Limitações declaradas

- A escolha de cadeia e a classificação de HETATM são heurísticas e exigem revisão.
  Nenhuma versão futura eliminará essa exigência; a v2 apenas fornecerá o
  diagnóstico que a torna mais rápida e rastreável.
- Estruturas apo, com ocupação parcial de Zn, ou com metal substituído (Cd, Co)
  não são tratadas. O critério 2,8–5,0 Å pressupõe sítio dinuclear íntegro.
- Confôrmeros alternativos: apenas `altLoc` em branco ou `A` são lidos.
- A tabela reflete o estado das estruturas na data de acesso ao PDB; revisões de
  depósito podem alterar numeração de resíduos.
