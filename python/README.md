# lammpsIO
`lammpsIO.py`を実行すればとりあえずデータを読み込めるようになってる。
```
python3 lammpsIO.py --input data.lammpstrj
```

基本的には、以下の様に1スナップショットずつ読み込んでいく。
ダイナミクスなどの計算がしたい場合は、別に配列を用意して、その中にデータを格納していく様にする。
```
import lammps_io as lmp
import numpy as np

lmp = lmp.lammpstrjReader('dump.lammpstrj')
```
このように`lmp`を経由してデータにアクセスしていく。
例えば、
```
while True:
    lmp.read_snapshot()
    print(lmp.timestep)
    print(lmp.num_atoms)
```
のようにすれば、スナップショットの最後まで読み込むことができる。
lmpmpstrjReaderの中身は以下の通り。
- `timestep`: タイムステップ
- `num_atoms`: 原子の数
- `box_bounds`: ボックスの境界条件
- `mol`: 分子の情報
- `types`: 原子の種類
- `coords`: 原子の座標
- `image_flags`: 周期境界条件のフラグ

## 速度
`ITEM: ATOMS id mol x y z ix iy iz` のファイルを読み込むのに、4.46フレーム/秒くらいの速度である。
ちょっと遅いかもしれないが、まあ、許容範囲かな。
型判定が多くなるやり方なので、もう少し高速化する方法があるかもしれない。

