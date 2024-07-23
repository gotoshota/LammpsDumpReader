# lammps_io
`lammps_io.py`を実行すればとりあえずデータを読み込めるようになってる。
```
python3 lammps_io.py --input data.lammpstrj
```

基本的には、以下の様に1スナップショットずつ読み込んでいく。
ダイナミクスなどの計算がしたい場合は、別に配列を用意して、その中にデータを格納していく様にする。
```
import lammps_io as lmp
import numpy as np

traj = lmp.LammpsTrajectory('data.lammpstrj')
```
このように`traj`を経由してデータにアクセスしていく。
例えば、
```
while True:
    snapshot = traj.read_snapshot()
    print(snapshot['box_bounds'])
    print(snapshot['atom_data'])
```
のようにすれば、スナップショットの最後まで読み込むことができる。
`read_snapshot`は辞書型の配列を返す。それぞれ、
- `timestep`: タイムステップ
- `num_atoms`: 原子の数
- `box_bounds`: ボックスの境界条件
- `atoms`: ATOM行の情報全て。空白文字でsplitしている。
注意すべきは、`atoms`行である。

## ["atoms"]の中身とアクセスの仕方
`snapshot["atoms"][idx_atom]`で、`idx_atom`番目の原子の情報にアクセスできる。
例えば、lammpsのdump形式が
`id x y z`である場合、
```
$ print(snapshot["atoms"][idx_atom])
[id, x, y, z]
```
と返ってくる。
つまり、このとき座標にアクセスしたければ、
```
$ print(snapshot["atoms"][idx_atom][1:4])
```
とすればよい。

dump形式によって、アクセスするインデックスが異なるので、注意。
将来的には、`snapshot["x"]`みたいな感じで自動的にパースするようにしたい。
やる気が出たら。
