このディレクトリには、サンプルファイルが含まれています。

- `example_protocol.yaml`
  - テスト実行用のプロトコルファイルです。
  - タンパク質として `1F47B.pdb` 、プローブとして `A11` を用います。
- `1F47B.pdb`
  - タンパク質の構造ファイルです。 PDB ID `1F47` の B chain を切り出しています。
- `A11.mol2` および `A11.pdb`
  - isopropanol プローブの構造ファイルです。
  - これらの構造は Gaussian を利用して作成されています。
    - 構造最適化および部分電荷の割り当てに利用した Gaussian の入力ファイルは `A11_opt.gjf` および `A11_elec.gjf` として同梱しています。
    - `g16 < A11_opt.gjf` とすることで、構造最適化を行うことが出来ます。
    - `A11_elec.gjf` は、例として構造最適化前の構造に対して部分電荷を割り当てています。
      - 実際には `A11_opt.gjf` で得られた構造を用いて部分電荷を割り当てるべきであることに注意してください。

-----

- `yanagisawa2022_protocol.yaml`
  - Yanagisawa et al. 2022 [^1] で用いられたプロトコルファイルです。
  - 再現実験を行う際には、このファイルを参照してください。
- `debug_protocol.yaml`
  - 開発時のデバッグ用プロトコルファイルです。
  - 通常の実行には使用しないでください。



  [^1]:**Keisuke Yanagisawa**, Ryunosuke Yoshino, Genki Kudo, Takatsugu Hirokawa. "Inverse Mixed-Solvent Molecular Dynamics for Visualization of the Residue Interaction Profile of Molecular Probes", *International Journal of Molecular Sciences*, **23**: 4749, 2022/04. DOI: [10.3390/ijms23094749](https://doi.org/10.3390/ijms23094749)