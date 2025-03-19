# LAMMPSトラジェクトリ読み書きライブラリ / LAMMPS Trajectory Reader/Writer Library

## 概要 / Overview

日本語 | [English](#english)

このライブラリは、LAMMPSの分子動力学シミュレーションの出力ファイル（dump形式）を読み込み、処理、書き込みするための機能を提供します。高分子シミュレーションなどのトラジェクトリデータを効率的に扱うために設計されています。FortranとPythonの両方で実装されています。

### 主な機能

- LAMMPSトラジェクトリファイルの読み込み
- LAMMPSトラジェクトリファイルの書き込み
- 読み込み専用と書き込み専用のインターフェース
- 新規ファイル作成および既存ファイルへの追記
- mol_idによるフィルタリング

## 使用方法

### Fortran版

#### コンパイル方法

Fortranコンパイラ（gfortran, ifortなど）を使用してコンパイルできます：

```bash
# モジュールのコンパイル
gfortran -c fortran/lammpsio.f90

# テストプログラムのコンパイル
gfortran -o test/lampstrjReader fortran/lammpsio.f90 test/main.f90
```

または、テストディレクトリのMakefileを使用する：

```bash
cd test
make
```

#### 実行方法

```bash
# テストプログラムの実行（デフォルトデータファイル使用）
./test/lampstrjReader

# 特定のファイルを指定して実行
./test/lampstrjReader /path/to/dump.lammpstrj
```

#### 基本的な使用例

```fortran
! 読み込み例
type(lammpstrjReader) :: reader
call reader%open("trajectory.lammpstrj")

do
  call reader%read()
  if (reader%end_of_file) exit
  
  ! ここでデータを処理
  print *, "Timestep:", reader%timestep
  print *, "Particles:", reader%nparticles
end do

call reader%close()

! 書き込み例
type(lammpstrjWriter) :: writer
call writer%create("output.lammpstrj")

writer%timestep = 1000
writer%nparticles = 10
! その他のデータを設定
call writer%write()
call writer%close()
```

### Python版

#### 実行方法

Python版は追加の依存パッケージのインストールなしで実行できます：

```bash
# 読み込みのみ
python python/lammpsIO.py --input /path/to/dump.lammpstrj

# 読み込みと書き込み
python python/lammpsIO.py --input /path/to/dump.lammpstrj --output /path/to/output.lammpstrj

# mol_id=1の粒子だけを書き出し
python python/lammpsIO.py --input /path/to/dump.lammpstrj --output /path/to/output.lammpstrj --mol_id 1
```

#### 基本的な使用例

```python
# 読み込み例
from lammpsIO import lammpstrjReader

reader = lammpstrjReader("trajectory.lammpstrj")

while True:
    reader.read()
    if reader.end_of_file:
        break
    
    # ここでデータを処理
    print(f"Timestep: {reader.timestep}")
    print(f"Particles: {reader.nparticles}")
    
reader.close()

# 書き込み例
from lammpsIO import lammpstrjWriter

writer = lammpstrjWriter()
writer.create("output.lammpstrj")

writer.timestep = 1000
writer.nparticles = 10
# 座標データを設定（3行×粒子数列の配列）
writer.coords = coords  # numpy array with shape (3, nparticles)
writer.write()
writer.close()

# mol_idによるフィルタリング例
for i in range(reader.nparticles):
    if reader.mol[i] == 1:  # mol_id = 1の粒子だけを処理
        # 処理内容
```

---

<a name="english"></a>

## English

This library provides functionality for reading, processing, and writing LAMMPS molecular dynamics simulation output files (dump format). It is designed for efficiently handling trajectory data from polymer simulations and other molecular systems. Both Fortran and Python implementations are available.

### Main Features

- Reading LAMMPS trajectory files
- Writing LAMMPS trajectory files
- Separate reader and writer interfaces
- Support for creating new files and appending to existing files
- Filtering by mol_id

## Usage

### Fortran Version

#### Compilation

You can compile using Fortran compilers (gfortran, ifort, etc.):

```bash
# Compile the module
gfortran -c fortran/lammpsio.f90

# Compile the test program
gfortran -o test/lampstrjReader fortran/lammpsio.f90 test/main.f90
```

Or using the Makefile in the test directory:

```bash
cd test
make
```

#### Running

```bash
# Run test program (using default data file)
./test/lampstrjReader

# Run with specific file
./test/lampstrjReader /path/to/dump.lammpstrj
```

#### Basic Usage Examples

```fortran
! Reading example
type(lammpstrjReader) :: reader
call reader%open("trajectory.lammpstrj")

do
  call reader%read()
  if (reader%end_of_file) exit
  
  ! Process data here
  print *, "Timestep:", reader%timestep
  print *, "Particles:", reader%nparticles
end do

call reader%close()

! Writing example
type(lammpstrjWriter) :: writer
call writer%create("output.lammpstrj")

writer%timestep = 1000
writer%nparticles = 10
! Set other data
call writer%write()
call writer%close()
```

### Python Version

#### Running

The Python version can be run without additional package installation:

```bash
# Reading only
python python/lammpsIO.py --input /path/to/dump.lammpstrj

# Reading and writing
python python/lammpsIO.py --input /path/to/dump.lammpstrj --output /path/to/output.lammpstrj

# Write only particles with mol_id=1
python python/lammpsIO.py --input /path/to/dump.lammpstrj --output /path/to/output.lammpstrj --mol_id 1
```

#### Basic Usage Examples

```python
# Reading example
from lammpsIO import lammpstrjReader

reader = lammpstrjReader("trajectory.lammpstrj")

while True:
    reader.read()
    if reader.end_of_file:
        break
    
    # Process data here
    print(f"Timestep: {reader.timestep}")
    print(f"Particles: {reader.nparticles}")
    
reader.close()

# Writing example
from lammpsIO import lammpstrjWriter

writer = lammpstrjWriter()
writer.create("output.lammpstrj")

writer.timestep = 1000
writer.nparticles = 10
# Set coordinate data (3 rows × nparticles columns array)
writer.coords = coords  # numpy array with shape (3, nparticles)
writer.write()
writer.close()

# Filtering by mol_id example
for i in range(reader.nparticles):
    if reader.mol[i] == 1:  # Process only particles with mol_id = 1
        # Processing code
```

## Documentation Generation

You can generate documentation using Doxygen:

```bash
doxygen Doxyfile
```

The generated documentation can be viewed at `doc/html/index.html`.
