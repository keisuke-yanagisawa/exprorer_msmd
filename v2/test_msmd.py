from unittest import TestCase
import MSMD

class TestMSMDSystemGenerator(TestCase):
  def __init__(self, *args, **kwargs):
    super(TestMSMDSystemGenerator, self).__init__(*args, **kwargs)
  def test_protein_only_system(self):
    "タンパク質しかない系を作るテスト。原子数で判断。"
    pass
  def test_probe_only_system(self):
    "プローブ1分子しかない系を作るテスト。原子数で判断。"
    pass
  def test_periodic_boundary(self):
    "境界条件を持つ系を作るテスト"
    pass
  def test_protein_and_single_probe(self):
    "タンパク質とプローブ1分子しかない系を作るテスト。原子数で判断。"
    pass
  def test_protein_and_multi_probes(self):
    "タンパク質とプローブ10分子からなる系を作るテスト。原子数で判断"
    pass
  def test_system_size(self):
    "タンパク質のx, y, z方面のもっとも長い距離 + 10A を1辺とするボックスを作るテスト"
    pass
  def test_no_molecule_passed_to_packmol_error(self):
    "packmolの入力として1つも分子が指定されていない場合はValueError"  
    pass
  def test_add_zero_molecules_packmol_errory(self):
    "packmolの入力として、合計で0この分子を追加する、と設定された場合はValueError"
    pass
  def test_zero_concentration_probe_warning(self):
    "probeの濃度が0 M と指定された場合は、proteinのみからなる系を作成し、warningを出す"
    pass
  def test_too_high_concentration_probe_warning(self):
    "probeの濃度が異常に高い（1 M より高い）場合は、系を作成し、warningを出す"
    pass
  def test_same_system_with_same_seed(self):
    "シード値が一緒の場合は常に同じ系を作る"
    pass
  def test_different_system_with_different_seed(self):
    "シード値が異なる場合には常に異なる系を作る"
    pass
  def test_negative_concentration_probe_error(self):
    "probe濃度が負を指定された場合はValueError"
    pass
  def test_malformed_protein_file_error(self):
    """
    タンパク質のファイルの中身が不正であった場合はValue?Error
    TODO: どのような不正かを考える必要がある
    """
    pass
  def test_no_protein_file_error(self):
    "タンパク質ファイルが存在しなかったらValueError"
    pass
  def test_malformed_probe_file_error(self):
    """
    プローブのファイルの中身が不正であった場合はValue?Error
    TODO: どのような不正かを考える必要がある
    """
    pass
  def test_no_probe_file_error(self):
    "プローブファイルが存在しなかったらValueError"
    pass
  def test_add_extra_particle_with_single_probe(self):
    "probeが1つだけ存在するシステムで、仮想原子が追加されている"
    pass
  def test_add_extra_particles_with_multi_probes(self):
    "probeが10こ存在するシステムで、仮想原子が10こ追加されている"
    pass
  def test_extra_particle_position(self):
    "仮想原子の座標がprobeの重原子の中心になっている"
    pass
  def test_repulsion_between_extra_particles(self):
    "仮想原子間に反発項が働いている"
    pass
  def test_ions_with_neutral_probe(self):
    "probeが中性の場合のイオンはCl-とNa+が同数"
  def test_ions_with_positive_probe(self):
    "probeが1価の陽イオンの場合はCl-の数はprobe数 + Na+ 数"
    pass
  def test_ions_with_negative_probe(self):
    "probeが1価の陰イオンの場合はNa+の数はprobe数 + Cl- 数"
    pass
  def test_protein_with_no_hydrogen_error(self):
    "タンパク質に水素が含まれていない場合はValueError"
    pass
  def test_probe_with_no_hydrogen_error(self):
    "プローブに水素が含まれていない場合はValueError"
    pass
  def test_sum_of_probe_partial_charge_is_not_zero_error(self):
    "プローブの部分電荷の合計が0でない場合はValueError"
    pass

class TestSimulationSequence(TestCase):
  def __init__(self, *args, **kwargs):
    super(TestSimulationSequence, self).__init__(*args, **kwargs)
  def test_read_sequence_minimization(self):
    pass
  def test_read_seuqnece_heating(self):
    pass
  def test_read_sequence_equilibration(self):
    pass
  def test_read_sequence_production(self):
    pass
