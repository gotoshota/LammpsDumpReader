# LAMMPSトラジェクトリ読み書きライブラリ / LAMMPS Trajectory Reader/Writer Library

## 概要 / Overview

日本語 | [English](#english)

このライブラリは、LAMMPSの分子動力学シミュレーションの出力ファイル（dump形式）を読み込み、処理、書き込みするための機能を提供します。高分子シミュレーションなどのトラジェクトリデータを効率的に扱うために設計されています。

### 主な機能

- LAMMPSトラジェクトリファイルの読み込み
- LAMMPSトラジェクトリファイルの書き込み
- 読み込み専用と書き込み専用のインターフェース
- 新規ファイル作成および既存ファイルへの追記

## 使用方法

### コンパイル方法

Fortranコンパイラ（gfortran, ifortなど）を使用してコンパイルできます：

```bash
# モジュールのコンパイル
gfortran -c fortran/lammpsio.f90

# テストプログラムのコンパイル
gfortran -c test/main.f90

# 実行ファイルの生成
gfortran lammpsio.o main.o -o test_lammps
```

### 基本的な使用例

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

## ドキュメント生成

Doxygenを使用してドキュメントを生成できます：

```bash
doxygen Doxyfile
```

生成されたドキュメントは `doc/html/index.html` で閲覧できます。

---

<a name="english"></a>

## English

This library provides functionality for reading, processing, and writing LAMMPS molecular dynamics simulation output files (dump format). It is designed for efficiently handling trajectory data from polymer simulations and other molecular systems.

### Main Features

- Reading LAMMPS trajectory files
- Writing LAMMPS trajectory files
- Separate reader and writer interfaces
- Support for creating new files and appending to existing files

## Usage

### Compilation

You can compile using Fortran compilers (gfortran, ifort, etc.):

```bash
# Compile the module
gfortran -c fortran/lammpsio.f90

# Compile the test program
gfortran -c test/main.f90

# Generate executable
gfortran lammpsio.o main.o -o test_lammps
```

### Basic Usage Examples

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

## Documentation Generation

You can generate documentation using Doxygen:

```bash
doxygen Doxyfile
```

The generated documentation can be viewed at `doc/html/index.html`.
