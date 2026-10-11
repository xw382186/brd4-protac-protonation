# BRD4 PROTAC 质子化计算

本仓库对应论文 *Computational Evaluation of Linker-Varied BRD4-Targeting PROTACs Using Deep-Learning Activity Classification and Molecular Dynamics* 最新补充信息中的质子化评估及 Supplementary Tables S4A–S4C。

[English README](README.md) · [详细计算流程与公式（英文）](docs/METHODS.md) · [可选状态枚举说明（英文）](docs/ENUMERATION.md)

仓库提供原始 RDKit 回归脚本、20 个标定输入、化合物 5–8 的输入 SMILES、修改后的 Dimorphite-DL 规则、全精度参考复算结果及便于运行的复算脚本。原始脚本保留不变，新增脚本用于组织输出与记录计算来源。

## 复算方法

SI 记载的核验环境为 **Python 3.11.16、RDKit 2026.03.5、NumPy 2.4.6**。数值 pKa 回归需要 NumPy 与 RDKit；Dimorphite-DL 用于独立的候选质子化状态枚举。

```bash
python -m pip install -r requirements.txt
python scripts/reproduce.py --output-dir results/reproduced
python -m unittest discover -s tests -v
```

每次复算应指定尚不存在的输出目录。脚本先校验原始计算脚本的 SHA256，再按其原有描述符、标定数据、拟合与预测步骤执行计算，将新结果写入指定目录，并记录实际软件版本。原始脚本最后的固定文件写入部分不执行，参考结果不会被覆盖。

可选的候选状态枚举应在**单独的环境**中运行，因为 Dimorphite-DL 2.0.2 的依赖要求为 `rdkit<2026`。这样可保留上述数值 pKa 复算环境：

```bash
python -m venv .venv-enumeration
source .venv-enumeration/bin/activate
python -m pip install -r requirements-enumeration.txt
python scripts/enumerate_states.py --inputs data/molecular_inputs.smi --output results/enumerated_states
```

已测试的可选环境使用 RDKit 2025.09.6。按照存档的定制规则和 SI 参数，5–8 每个化合物均返回一个净形式电荷为 +1 的唯一状态。另一次 RDKit 2026.03.5 兼容性运行获得了相同的状态记录。枚举不产生这里报告的数值 pKa，也不为各状态分配概率，详见 [ENUMERATION.md](docs/ENUMERATION.md)。

## 计算了什么？

1. 根据 RDKit 分子连接关系选择符合规则的脂肪胺位点：5、6 为 azetidine N，7、8 为 piperidine N。
2. 用原脚本 `REF` 列表中的 20 个数值输入拟合普通最小二乘回归。
3. 以距离衰减的诱导描述符、氮取代度指示变量及含氮环描述符预测所选胺的共轭酸 pKa。
4. 将未四舍五入的 pKa 代入单个位点 Henderson–Hasselbalch 公式，计算 pH 7.4 下所选胺的质子化百分比。
5. 导出全精度预测、拟合结果与逐条留一法结果。

| 化合物 | 所选氮及输入 SMILES 中的零基原子编号 | 估计 pKa | pH 7.4 下所选胺质子化比例 | 工作净电荷 |
| --- | --- | ---: | ---: | ---: |
| 5 | Azetidine N；26 | 9.59 | 99.4% | +1 |
| 6 | Azetidine N；27 | 10.01 | 99.8% | +1 |
| 7 | Piperidine N；17 | 9.51 | 99.2% | +1 |
| 8 | Piperidine N；27 | 9.93 | 99.7% | +1 |

上述数值按 SI Table S4A 的精度显示。其余化合物沿用工作中性电荷赋值，SI 未报告其数值 pKa。原子编号对应 `data/molecular_inputs.smi` 中未经改写的输入字符串。

## 文件说明

| 路径 | 内容 |
| --- | --- |
| `source/pka_rdkit_lfer.py` | 保留不变的原始计算脚本 |
| `scripts/reproduce.py` | 可移植的数值复算脚本，记录输出与计算来源 |
| `scripts/enumerate_states.py` | 可选枚举复现脚本，在临时包副本中应用规则 |
| `data/calibration_inputs.csv` | 原脚本的 20 个标定输入，对应 SI Table S4B |
| `data/molecular_inputs.smi` | 化合物 5–8 的输入 SMILES |
| `data/si_table_s4a.csv` | SI 的工作电荷赋值及四舍五入后的胺位点结果 |
| `rules/` | 原始和修改后的 SMARTS 规则表、二者重建差异及记载的枚举调用参数 |
| `results/reference/` | 已提供的全精度数值复算结果、留一法结果及来源记录 |
| `docs/METHODS.md` | 公式、参数、位点筛选、逐步流程及结果解释 |
| `docs/ENUMERATION.md` | 可选枚举的环境、执行步骤、输出与验证 |
| `tests/` | 可重复性检查 |

## 数值来源及解释

化合物的数值 pKa 来自原回归脚本的预测。用于标定的 20 个 pKa 则直接来自脚本的 `REF` 列表，SI Table S4B 也列出了这些输入。现有材料没有为每条输入指定一次文献出处；仓库忠实保留输入，不替换数值或补造来源。

该回归只用分子连接关系，不使用 MD 坐标。Dimorphite-DL 候选状态枚举与数值回归是独立步骤。Henderson–Hasselbalch 百分比只描述单个位点模型下所选胺的质子化；工作净电荷是实际采用的配体状态赋值，不是按总体分子状态分布加权的电荷。留一法误差用于描述标定集内的预测表现，不代表各 PROTAC 预测的置信区间。

原始脚本中的历史比较字段与输出标签保留原貌。`OLD_ANCHOR` 不参与拟合或预测；Gasteiger 氮电荷和环大小只是诊断输出，不是回归特征。完整步骤与原脚本校验值见 [METHODS.md](docs/METHODS.md)。
