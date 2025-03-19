#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@file lammpsIO.py
@brief LAMMPSトラジェクトリファイルの読み込みと書き込みを行うモジュール
@details このモジュールはLAMMPSの分子動力学シミュレーションの出力ファイル（dump形式）を
読み込み、処理、書き込みするための機能を提供します。lammpsio.f90のPython実装版です。
@author 開発者名
@date 作成日
"""

import numpy as np
import math


class AtomIndex:
    """
    @brief LAMMPSダンプファイル内の各カラム（id, x, y, z など）の位置情報を保持するクラス
    @details LAMMPSダンプファイル内の各カラム（原子ID、座標、イメージフラグなど）の
    位置情報を保持します。カラムのインデックスは0から始まります。
    """
    def __init__(self):
        """
        @brief AtomIndexオブジェクトを初期化します
        @details 全てのインデックスは0で初期化されます
        """
        self.id = 0    #!< 原子IDのカラムインデックス
        self.mol = 0   #!< 分子IDのカラムインデックス
        self.type = 0  #!< 原子タイプのカラムインデックス
        self.xu = 0    #!< アンラップしたx座標のカラムインデックス
        self.yu = 0    #!< アンラップしたy座標のカラムインデックス
        self.zu = 0    #!< アンラップしたz座標のカラムインデックス
        self.x = 0     #!< x座標のカラムインデックス
        self.y = 0     #!< y座標のカラムインデックス
        self.z = 0     #!< z座標のカラムインデックス
        self.ix = 0    #!< xイメージフラグのカラムインデックス
        self.iy = 0    #!< yイメージフラグのカラムインデックス
        self.iz = 0    #!< zイメージフラグのカラムインデックス
        self.xs = 0    #!< スケーリングされたx座標のカラムインデックス
        self.ys = 0    #!< スケーリングされたy座標のカラムインデックス
        self.zs = 0    #!< スケーリングされたz座標のカラムインデックス


class lammpstrj:
    """
    @brief LAMMPSトラジェクトリデータを管理する基本クラス
    @details トラジェクトリファイルの読み込みと書き込み機能を提供する親クラスです。
    このクラスは直接使用せず、lammpstrjReaderかlammpstrjWriterを使用してください。
    """
    def __init__(self):
        """
        @brief lammpstrjオブジェクトを初期化します
        @details 全てのフィールドが初期化され、後で値が設定されます
        """
        self.filename = None      #!< ファイル名 
        self.file = None          #!< ファイルオブジェクト
        self.end_of_file = False  #!< EOFフラグ
        self.has_atom_idx = False #!< 原子インデックスが初期化されているかのフラグ
        self.is_writing = False   #!< 書き込みモードフラグ
        self.timestep = None      #!< 現在のタイムステップ
        self.nparticles = None    #!< 粒子数（Fortranでのnparticles）
        self.box_bounds = None    #!< シミュレーションボックスの境界 (3,3)のnumpy配列
        self.coords = None        #!< 座標データ (3, nparticles)のnumpy配列
        self.image_flags = None   #!< イメージフラグ (3, nparticles)のnumpy配列
        self.id = None            #!< 原子ID (nparticles)のnumpy配列
        self.mol = None           #!< 分子ID (nparticles)のnumpy配列
        self.type = None          #!< 原子タイプ (nparticles)のnumpy配列
        self.atom_idx = AtomIndex() #!< カラムインデックス情報
    
    def open(self, filename, mode='read'):
        """
        @brief トラジェクトリファイルを開きます
        @param filename 開くファイル名
        @param mode オプションのファイルモード ('read', 'write', 'append')
        @exception ValueError 不明なモードが指定された場合
        """
        self.filename = filename
        
        if mode == 'read':
            self.is_writing = False
            self.file = open(filename, 'r')
        elif mode == 'write':
            self.is_writing = True
            self.file = open(filename, 'w')
        elif mode == 'append':
            self.is_writing = True
            self.file = open(filename, 'a')
        else:
            raise ValueError(f"エラー: 不明なモードです: {mode}")
    
    def close(self):
        """
        @brief トラジェクトリファイルを閉じます
        @details ファイルが開かれている場合、それを閉じてNoneに設定します
        """
        if self.file:
            self.file.close()
            self.file = None
    
    def __del__(self):
        """
        @brief オブジェクトが削除されるときにファイルを閉じます
        @details デストラクタはガベージコレクション時に呼ばれます
        """
        self.close()
    
    def read(self):
        """
        @brief トラジェクトリファイルから次のフレームを読み込みます
        @details ファイルから次のタイムステップのデータを読み込み、クラスのフィールドを更新します。
        これには粒子の座標、イメージフラグ、各種IDなどが含まれます。
        @exception Exception ファイル読み込み中にエラーが発生した場合
        """
        if self.is_writing:
            print("エラー: このトラジェクトリは書き込みモードのため読み込めません")
            return
        
        if self.end_of_file:
            return
        
        try:
            line = self.file.readline()
            if not line:
                self.end_of_file = True
                return
            
            # タイムステップを読み込み
            self.timestep = int(self.file.readline().strip())
            
            # 粒子数を読み込み
            self.file.readline()  # ITEM: NUMBER OF ATOMS
            self.nparticles = int(self.file.readline().strip())
            
            # ボックス境界を読み込み
            self.file.readline()  # ITEM: BOX BOUNDS
            self.box_bounds = np.array([
                list(map(float, self.file.readline().strip().split())) for _ in range(3)
            ])
            
            # ATOMS行を読み込み、カラム情報を解析
            line = self.file.readline()  # ITEM: ATOMS
            
            if not self.has_atom_idx:
                # カラム情報を解析
                atom_header_parts = line.strip().split()
                # "ITEM:" と "ATOMS" をスキップ
                ncols = len(atom_header_parts) - 2
                
                for i in range(ncols):
                    col_name = atom_header_parts[i + 2].strip()
                    if col_name == "id":
                        self.atom_idx.id = i
                    elif col_name == "mol":
                        self.atom_idx.mol = i
                    elif col_name == "type":
                        self.atom_idx.type = i
                    elif col_name == "xu":
                        self.atom_idx.xu = i
                    elif col_name == "yu":
                        self.atom_idx.yu = i
                    elif col_name == "zu":
                        self.atom_idx.zu = i
                    elif col_name == "x":
                        self.atom_idx.x = i
                    elif col_name == "y":
                        self.atom_idx.y = i
                    elif col_name == "z":
                        self.atom_idx.z = i
                    elif col_name == "ix":
                        self.atom_idx.ix = i
                    elif col_name == "iy":
                        self.atom_idx.iy = i
                    elif col_name == "iz":
                        self.atom_idx.iz = i
                    elif col_name == "xs":
                        self.atom_idx.xs = i
                    elif col_name == "ys":
                        self.atom_idx.ys = i
                    elif col_name == "zs":
                        self.atom_idx.zs = i
                
                self.has_atom_idx = True
                
                # 必要な配列を確保
                self.coords = np.zeros((3, self.nparticles))
                if self.atom_idx.id > 0:
                    self.id = np.zeros(self.nparticles, dtype=int)
                if self.atom_idx.mol > 0:
                    self.mol = np.zeros(self.nparticles, dtype=int)
                if self.atom_idx.type > 0:
                    self.type = np.zeros(self.nparticles, dtype=int)
                if self.atom_idx.ix > 0:
                    self.image_flags = np.zeros((3, self.nparticles), dtype=int)
            
            # 粒子データを読み込み
            for i in range(self.nparticles):
                line = self.file.readline()
                parts = line.strip().split()
                
                # 座標の読み込み
                if self.atom_idx.x > 0:
                    self.coords[0, i] = float(parts[self.atom_idx.x])
                elif self.atom_idx.xu > 0:
                    self.coords[0, i] = float(parts[self.atom_idx.xu])
                elif self.atom_idx.xs > 0:
                    self.coords[0, i] = float(parts[self.atom_idx.xs])
                
                if self.atom_idx.y > 0:
                    self.coords[1, i] = float(parts[self.atom_idx.y])
                elif self.atom_idx.yu > 0:
                    self.coords[1, i] = float(parts[self.atom_idx.yu])
                elif self.atom_idx.ys > 0:
                    self.coords[1, i] = float(parts[self.atom_idx.ys])
                
                if self.atom_idx.z > 0:
                    self.coords[2, i] = float(parts[self.atom_idx.z])
                elif self.atom_idx.zu > 0:
                    self.coords[2, i] = float(parts[self.atom_idx.zu])
                elif self.atom_idx.zs > 0:
                    self.coords[2, i] = float(parts[self.atom_idx.zs])
                
                # ID, mol, type の読み込み
                if self.atom_idx.id > 0 and self.id is not None:
                    self.id[i] = int(parts[self.atom_idx.id])
                if self.atom_idx.mol > 0 and self.mol is not None:
                    self.mol[i] = int(parts[self.atom_idx.mol])
                if self.atom_idx.type > 0 and self.type is not None:
                    self.type[i] = int(parts[self.atom_idx.type])
                
                # イメージフラグの読み込み
                if self.atom_idx.ix > 0 and self.image_flags is not None:
                    self.image_flags[0, i] = int(parts[self.atom_idx.ix])
                if self.atom_idx.iy > 0 and self.image_flags is not None:
                    self.image_flags[1, i] = int(parts[self.atom_idx.iy])
                if self.atom_idx.iz > 0 and self.image_flags is not None:
                    self.image_flags[2, i] = int(parts[self.atom_idx.iz])
        
        except Exception as e:
            print(f"Error reading file: {e}")
            self.end_of_file = True
    
    def write(self):
        """
        @brief トラジェクトリデータをファイルに書き込みます
        @details 現在のフレームデータをLAMMPS形式でファイルに書き込みます。
        タイムステップ、粒子数、ボックス境界、各粒子のデータが含まれます。
        """
        if not self.is_writing:
            print("エラー: このトラジェクトリは読み込みモードのため書き込めません")
            return
        
        # タイムステップ情報を書き込み
        self.file.write("ITEM: TIMESTEP\n")
        self.file.write(f"{self.timestep}\n")
        
        # 粒子数を書き込み
        self.file.write("ITEM: NUMBER OF ATOMS\n")
        self.file.write(f"{self.nparticles}\n")
        
        # ボックス境界情報を書き込み
        self.file.write("ITEM: BOX BOUNDS pp pp pp\n")
        for i in range(3):
            bounds = " ".join(map(str, self.box_bounds[i, :]))
            self.file.write(f"{bounds}\n")
        
        # 原子データのヘッダーを書き込み
        header = "ITEM: ATOMS id"
        if self.type is not None:
            header += " type"
        if self.mol is not None:
            header += " mol"
        header += " x y z"
        if self.image_flags is not None:
            header += " ix iy iz"
        self.file.write(f"{header}\n")
        
        # 各原子のデータを書き込み
        for i in range(self.nparticles):
            line = ""
            
            # ID を書き込み (必須)
            if self.id is not None:
                line += f"{self.id[i]} "
            else:
                line += f"{i+1} "  # IDがない場合は位置を使用
            
            # タイプを書き込み (オプション)
            if self.type is not None:
                line += f"{self.type[i]} "
            
            # 分子 ID を書き込み (オプション)
            if self.mol is not None:
                line += f"{self.mol[i]} "
            
            # 座標を書き込み (必須)
            line += f"{self.coords[0, i]} {self.coords[1, i]} {self.coords[2, i]}"
            
            # イメージフラグを書き込み (オプション)
            if self.image_flags is not None:
                line += f" {self.image_flags[0, i]} {self.image_flags[1, i]} {self.image_flags[2, i]}"
            
            self.file.write(f"{line}\n")


class lammpstrjReader(lammpstrj):
    """
    @brief 読み込み専用のLAMMPSトラジェクトリ管理クラス
    @details lammpstrjを継承し、読み込み機能に特化した子クラスです。
    LAMMPSトラジェクトリファイルを読み込むために使用します。
    """
    def __init__(self, filename=None):
        """
        @brief lammpstrjReaderオブジェクトを初期化します
        @param filename オプションのファイル名。指定された場合、そのファイルが開かれます
        """
        super().__init__()
        if filename:
            self.open(filename)
    
    def open(self, filename, mode='read'):
        """
        @brief 読み込み専用でトラジェクトリファイルを開きます
        @param filename 開くファイル名
        @param mode 無視されるパラメータ（互換性のために残されています）
        @details modeパラメータは無視され、常に'read'モードでファイルを開きます
        """
        super().open(filename, 'read')  # 常に読み込みモードで開く
    
    def read_next_frame(self):
        """
        @brief 次のフレームを読み込みます（旧バージョンとの互換性のためのメソッド）
        @return self.read()の戻り値
        @details 旧バージョンとの互換性のために提供されているメソッドです。
        内部的にはread()メソッドを呼び出します。
        """
        return self.read()


class lammpstrjWriter(lammpstrj):
    """
    @brief 書き込み専用のLAMMPSトラジェクトリ管理クラス
    @details lammpstrjを継承し、書き込み機能に特化した子クラスです。
    LAMMPSトラジェクトリファイルを作成または追記するために使用します。
    """
    def __init__(self):
        """
        @brief lammpstrjWriterオブジェクトを初期化します
        """
        super().__init__()
    
    def open(self, filename, mode='write'):
        """
        @brief 書き込み専用でトラジェクトリファイルを開きます
        @param filename 開くファイル名
        @param mode オプションのファイルモード ('write'または'append')
        @details 'read'モードが指定された場合、警告を出力して'write'モードに変更します
        """
        if mode == 'read':
            print("警告: Writerは読み込みモードではオープンできません。書き込みモードに変更します。")
            mode = 'write'
        
        super().open(filename, mode)
    
    def create(self, filename):
        """
        @brief 新規書き込み用のトラジェクトリファイルを作成します
        @param filename 作成するファイル名
        @details 既存のファイルがある場合は上書きされます
        """
        self.open(filename, 'write')
    
    def append(self, filename):
        """
        @brief 既存のトラジェクトリファイルに追記するために開きます
        @param filename 追記するファイル名
        @details 指定したファイルが存在しない場合は新規作成されます
        """
        self.open(filename, 'append')


def write_lammpstrj(filename, coords, box_bounds, image_flags=None):
    """
    @brief シンプルなヘルパー関数：LAMMPSトラジェクトリファイルを書き込みます
    @param filename 出力ファイル名
    @param coords 座標データ（N, 3）のnumpy配列
    @param box_bounds ボックス境界（3, 2）または（3, 3）のnumpy配列
    @param image_flags オプションのイメージフラグ（N, 3）のnumpy配列
    @details この関数はシンプルなヘルパー関数で、lammpstrjWriterクラスを使わずに
    直接ファイルに書き込みます。タイムステップは0に設定されます。
    """
    with open(filename, "w") as f:
        f.write("ITEM: TIMESTEP\n")
        f.write("0\n")
        f.write("ITEM: NUMBER OF ATOMS\n")
        f.write(f"{len(coords)}\n")
        f.write("ITEM: BOX BOUNDS\n")
        # box_bounds は (3, 2) or (3, 3) の numpy 配列
        f.write(" ".join(map(str, box_bounds[0, :])) + "\n")
        f.write(" ".join(map(str, box_bounds[1, :])) + "\n")
        f.write(" ".join(map(str, box_bounds[2, :])) + "\n")
        f.write("ITEM: ATOMS id x y z")
        if image_flags is not None:
            f.write(" ix iy iz")
        f.write("\n")
        for i, coord in enumerate(coords):
            line = f"{i} {coord[0]} {coord[1]} {coord[2]}"
            if image_flags is not None:
                line += f" {image_flags[i, 0]} {image_flags[i, 1]} {image_flags[i, 2]}"
            f.write(line + "\n")


class LammpsData:
    """
    LAMMPSのデータファイルを読み込み，データを格納するクラス

    Attributes:
        filename (str): LAMMPSデータファイルのパス
        box (Box): ボックス情報を保持するオブジェクト
        atoms (Atom): 原子情報を保持するオブジェクト
        bonds (Bond): ボンド情報を保持するオブジェクト
        angles (Angle): 角度情報を保持するオブジェクト
        現在，2面角 (Dihedrals) や不整合角 (Impropers) は未対応． また，将来の拡張も考慮していない．
    """

    def __init__(self, filename=None):
        """
        LammpsDataオブジェクトを初期化し，ファイルが指定されていればデータを読み込む

        Args:
            filename (str): LAMMPSデータファイルのパス
        """
        self.filename = filename
        self.atoms = self.Atom()
        self.bonds = self.Bond()
        self.angles = self.Angle()
        self.box = self.Box()
        self.masses = self.Mass()

        if filename is not None:
            self.read(filename)

    class Mass:
        """
        原子の質量情報を保持するクラス

        Attributes:
            id (list of int): 原子タイプのリスト
            mass (list of float): 原子の質量のリスト
        """

        def __init__(self):
            self.id = []
            self.mass = []

    class Box:
        """
        ボックス情報を保持するクラス

        Attributes:
            x (tuple of float): x軸の範囲 (xlo, xhi)
            y (tuple of float): y軸の範囲 (ylo, yhi)
            z (tuple of float): z軸の範囲 (zlo, zhi)
            lx (float): x軸の長さ
            ly (float): y軸の長さ
            lz (float): z軸の長さ
        """

        def __init__(self):
            self.x = (0.0, 0.0)
            self.y = (0.0, 0.0)
            self.z = (0.0, 0.0)
            self.lx = 0.0
            self.ly = 0.0
            self.lz = 0.0

    class Atom:
        """
        原子情報を保持するクラス

        Attributes:
            id (list of int): 原子IDのリスト
            mol_id (list of int): 分子IDのリスト
            type (list of int): 原子タイプのリスト
            coords (list of list of float): 座標 (x, y, z) のリスト
            num_atoms (int): 原子の総数
            num_mols (int): 分子の総数
            num_types (int): 原子タイプの数
            image_flag (list of tuple of float): 原子のイメージフラグのリスト
        """

        def __init__(self):
            self.id = []
            self.mol_id = []
            self.type = []
            self.coords = []
            self.num_atoms = 0
            self.num_mols = 0
            self.num_types = 0
            self.image_flag = []  # list of tuple of float

    class Bond:
        """
        ボンド情報を保持するクラス

        Attributes:
            id (list of int): ボンドIDのリスト
            type (list of int): ボンドタイプのリスト
            atoms (list of tuple): ボンドを構成する原子のリスト
            num_bonds (int): ボンドの総数
            num_types (int): ボンドタイプの数
        """

        def __init__(self):
            self.id = []
            self.type = []
            self.atoms = []
            self.num_bonds = 0
            self.num_types = 0

    class Angle:
        """
        角度情報を保持するクラス

        Attributes:
            id (list of int): 角度IDのリスト
            type (list of int): 角度タイプのリスト
            atoms (list of tuple): 角度を構成する原子のリスト
            num_angles (int): 角度の総数
            num_types (int): 角度タイプの数
        """

        def __init__(self):
            self.id = []
            self.type = []
            self.atoms = []
            self.num_angles = 0
            self.num_types = 0

    def read(self, filename):
        """
        LAMMPSデータファイルからデータを読み込み，各属性に格納する

        Args:
            filename (str): 読み込むファイル名
        """
        if filename is not None:
            self.filename = filename

        with open(self.filename, "r") as f:
            lines = f.readlines()

        # --- 1. ヘッダー部からボックス情報や型数を取得 ---
        for line in lines:
            line_strip = line.strip()
            if "xlo" in line_strip and "xhi" in line_strip:
                parts = line_strip.split()
                self.box.x = (float(parts[0]), float(parts[1]))
                self.box.lx = float(parts[1]) - float(parts[0])
            elif "ylo" in line_strip and "yhi" in line_strip:
                parts = line_strip.split()
                self.box.y = (float(parts[0]), float(parts[1]))
                self.box.ly = float(parts[1]) - float(parts[0])
            elif "zlo" in line_strip and "zhi" in line_strip:
                parts = line_strip.split()
                self.box.z = (float(parts[0]), float(parts[1]))
                self.box.lz = float(parts[1]) - float(parts[0])
            elif "atoms" in line_strip:
                parts = line_strip.split()
                self.atoms.num_atoms = int(parts[0])
            elif "bonds" in line_strip:
                parts = line_strip.split()
                self.bonds.num_bonds = int(parts[0])
            elif "angles" in line_strip:
                parts = line_strip.split()
                self.angles.num_angles = int(parts[0])
            elif "atom types" in line_strip:
                parts = line_strip.split()
                self.atoms.num_types = int(parts[0])
            elif "bond types" in line_strip:
                parts = line_strip.split()
                self.bonds.num_types = int(parts[0])
            elif "angle types" in line_strip:
                parts = line_strip.split()
                self.angles.num_types = int(parts[0])

        # --- 2. セクション毎にデータをパース ---
        section_names = [
            "Masses",
            "Atoms",
            "Bonds",
            "Angles",
            "Dihedrals",
            "Impropers",
            "Pair Coeffs",
            "Bond Coeffs",
            "Angle Coeffs",
            "Dihedral Coeffs",
            "Improper Coeffs",
        ]
        current_section = None
        i = 0
        while i < len(lines):
            line = lines[i].strip()

            # セクションヘッダーの検出
            if any(line.startswith(sec) for sec in section_names):
                # セクション名を現在のセクションとして記憶
                current_section = line.split()[0]
                # セクション名の行とその次の空行（またはコメント行）をスキップ
                i += 2
                continue

            # セクション外なら次の行へ
            if current_section is None:
                i += 1
                continue

            # 空行やコメント行はスキップ
            if not line or line.startswith("#"):
                i += 1
                continue

            # 新たなセクションが始まった場合は current_section をリセット
            if any(line.startswith(sec) for sec in section_names):
                current_section = None
                continue

            # データ行のパース
            parts = line.split()
            if current_section == "Atoms":
                # 例: "1 1 1 0.0 0.0 0.0 ..." (原子ID, 分子ID, タイプ, x, y, z, ...)
                if len(parts) >= 6:
                    self.atoms.id.append(int(parts[0]))
                    self.atoms.mol_id.append(int(parts[1]))
                    self.atoms.type.append(int(parts[2]))
                    self.atoms.coords.append(
                        [float(parts[3]), float(parts[4]), float(parts[5])]
                    )
                    self.atoms.image_flag.append(
                        (float(parts[6]), float(parts[7]), float(parts[8]))
            elif current_section == "Masses":
                if len(parts) >= 2:
                    self.masses.id.append(int(parts[0]))
                    self.masses.mass.append(float(parts[1]))
            elif current_section == "Bonds":
                if len(parts) >= 4:
                    self.bonds.id.append(int(parts[0]))
                    self.bonds.type.append(int(parts[1]))
                    self.bonds.atoms.append((int(parts[2]), int(parts[3])))
            # TODO: 他のセクション (Masses, Angle, など) のパース処理も同様に追加

            i += 1
        # 被ってない，mol_id の数を数える
        self.atoms.num_mols = len(set(self.atoms.mol_id))
        # Atom, Bond セクションのデータを sort
        (
            self.atoms.id,
            self.atoms.mol_id,
            self.atoms.type,
            self.atoms.coords,
            self.atoms.image_flag,
        ) = map(
            list,
            zip(
                *sorted(
                    zip(
                        self.atoms.id,
                        self.atoms.mol_id,
                        self.atoms.type,
                        self.atoms.coords,
                        self.atoms.image_flag,
                    ),
                    key=lambda x: x[0],
                )
            ),
        )
        self.bonds.id, self.bonds.type, self.bonds.atoms = map(
            list,
            zip(
                *sorted(
                    zip(self.bonds.id, self.bonds.type, self.bonds.atoms),
                    key=lambda x: x[0],
                )
            ),
        )

    def write(self, filename):
        """
        LammpsDataオブジェクトのデータをLAMMPSデータファイルとして出力

        Args:
            filename (str): 出力ファイル名
        """
        with open(filename, "w") as f:
            # ヘッダー部
            f.write(f"LAMMPS data file generated by LammpsData class\n\n")
            f.write(f"{self.atoms.num_atoms} atoms\n")
            f.write(f"{self.bonds.num_bonds} bonds\n")
            f.write(f"{self.atoms.num_types} atom types\n")
            f.write(f"{self.bonds.num_types} bond types\n")
            f.write("\n")
            f.write(f"{self.box.x[0]} {self.box.x[1]} xlo xhi\n")
            f.write(f"{self.box.y[0]} {self.box.y[1]} ylo yhi\n")
            f.write(f"{self.box.z[0]} {self.box.z[1]} zlo zhi\n\n")

            # Masses セクション
            f.write("Masses\n\n")
            for i in range(len(self.masses.id)):
                f.write(f"{self.masses.id[i]} {self.masses.mass[i]}\n")
            f.write("\n")

            # Atoms セクション
            f.write(f"Atoms\n\n")
            for i in range(self.atoms.num_atoms):
                f.write(
                    f"{self.atoms.id[i]} {self.atoms.mol_id[i]} {
                        self.atoms.type[i]} "
                    f"{self.atoms.coords[i][0]} {self.atoms.coords[i][1]} "
                    f"{self.atoms.coords[i][2]} "
                    f"{self.atoms.image_flag[i][0]} {self.atoms.image_flag[i][1]} {
                        self.atoms.image_flag[i][2]}\n"
                )
            f.write("\n")
            f.write("Bonds\n\n")
            for i in range(self.bonds.num_bonds):
                f.write(
                    f"{self.bonds.id[i]} {self.bonds.type[i]} "
                    f"{self.bonds.atoms[i][0]} {self.bonds.atoms[i][1]}\n"
                )

    def __str__(self):
        return f"LammpsData({self.filename})"

    def __repr__(self):
        return self.__str__()

    def polyWrap(self):
        """
        ポリマーの原子座標を周期境界条件でラップする
        このとき，結合が切れないように原子座標を調整し，
        分子の重心（COM）がシミュレーションセル内に収まるように全体を平行移動する
        """
        # 分子IDは 1 から始まる
        for i in range(1, self.atoms.num_mols + 1):
            # 分子IDが i の原子のインデックスを取得
            idx = [j for j, mol_id in enumerate(self.atoms.mol_id) if mol_id == i]
            if not idx:
                continue

            # シミュレーションセル各方向の境界とセルサイズを取得
            # ※ここではセル情報を self.box.x, self.box.y, self.box.z としてアクセスしています
            box_bounds = {"x": self.box.x, "y": self.box.y, "z": self.box.z}
            box_lengths = {
                axis: box_bounds[axis][1] - box_bounds[axis][0] for axis in box_bounds
            }

            # 1. 分子内で「アンラップ」する（隣接原子との連続性を保持）
            # 最初の原子はそのまま採用
            unwrapped_coords = [list(self.atoms.coords[idx[0]])]
            unwrapped_image_flag = [list(self.atoms.image_flag[idx[0]])]
            for k in range(1, len(idx)):
                prev_coord = unwrapped_coords[k - 1]
                current_coord = list(self.atoms.coords[idx[k]])
                new_coord = []
                new_image = []
                for d, axis in enumerate(["x", "y", "z"]):
                    lo, hi = box_bounds[axis]
                    L = box_lengths[axis]
                    diff = current_coord[d] - prev_coord[d]
                    # 補正：差がセルサイズの半分を超えているかどうかで判断
                    if diff > 0.5 * L:
                        shift_val = 1
                    elif diff < -0.5 * L:
                        shift_val = -1
                    else:
                        shift_val = 0
                    new_coord.append(current_coord[d] - shift_val * L)
                    new_image.append(float(shift_val))
                unwrapped_coords.append(new_coord)
                # image_flag はタプルの場合が多いので、リスト内各成分同士を足し合わせたタプルにする
                unwrapped_image_flag.append(
                    tuple(
                        float(a) + float(b)
                        for a, b in zip(self.atoms.image_flag[idx[k]], new_image)
                    )
                )

            # 2. 分子の重心（COM）を計算（各原子の質量が等しいと仮定）
            n_atoms = len(unwrapped_coords)
            com = [
                sum(coord[d] for coord in unwrapped_coords) / n_atoms for d in range(3)
            ]

            # 3. COM がセル内に入るように全体を平行移動するシフト量を算出
            shift = []  # 座標のシフト（float 値）
            shift_image = []  # 画像フラグのシフト（整数値）
            for d, axis in enumerate(["x", "y", "z"]):
                lo, hi = box_bounds[axis]
                L = box_lengths[axis]
                # modulo 演算で COM を [lo, hi) に写像
                new_com = ((com[d] - lo) % L) + lo
                shift_val = new_com - com[d]
                shift.append(shift_val)
                # 負の場合も正しく扱うため math.floor を利用
                shift_image.append(math.floor((com[d] - lo) / L))

            # 4. 各原子座標と image_flag にシフトを適用
            wrapped_coords = []
            new_image_flag = []
            for coord in unwrapped_coords:
                new_coord = [coord[d] + shift[d] for d in range(3)]
                wrapped_coords.append(new_coord)
            for image_flag in unwrapped_image_flag:
                new_img = [
                    float(image_flag[d]) + float(shift_image[d]) for d in range(3)
                ]
                new_image_flag.append(new_img)

            # 5. 元の座標リストと image_flag を更新
            for index, new_coord in zip(idx, wrapped_coords):
                self.atoms.coords[index] = new_coord
            for index, new_img in zip(idx, new_image_flag):
                self.atoms.image_flag[index] = tuple(new_img)


if __name__ == "__main__":
    import argparse
    import time

    # 使用例
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, required=True)
    parser.add_argument("--output", type=str, help="出力ファイル名（指定すると書き込みテストも実行）")
    parser.add_argument("--mol_id", type=int, default=None, help="指定したmol_idの粒子だけを書き出し")

    args = parser.parse_args()
    
    # トラジェクトリの読み込みテスト
    print(f"入力ファイル: {args.input} の読み込みテスト")
    reader = lammpstrjReader(args.input)

    ts = time.time()
    nframes = 0
    
    # 書き込みテスト用の設定
    writer = None
    if args.output:
        print(f"出力ファイル: {args.output} への書き込みテスト")
        writer = lammpstrjWriter()
        writer.create(args.output)
    
    # 全フレームを処理
    while True:
        reader.read()
        if reader.end_of_file:
            break
        
        # フレームごとの処理
        print(f"Timestep: {reader.timestep}, Number of atoms: {reader.nparticles}")
        if reader.coords is not None:
            print(f"coords of atom 1: {reader.coords[:, 0]}")
            if reader.image_flags is not None:
                print(f"image_flags of atom 1: {reader.image_flags[:, 0]}")
        
        # 書き込みテスト
        if writer:
            if args.mol_id is not None and reader.mol is not None:
                # 指定されたmol_idの粒子だけを抽出して書き出し
                mol_mask = reader.mol == args.mol_id
                count_mol = np.sum(mol_mask)
                
                if count_mol > 0:
                    print(f"mol_id={args.mol_id}の粒子数: {count_mol}")
                    
                    # 新しいデータを作成
                    writer.timestep = reader.timestep
                    writer.nparticles = count_mol
                    writer.box_bounds = reader.box_bounds
                    
                    # 必要な配列を割り当て
                    if writer.coords is None or writer.coords.shape[1] != count_mol:
                        writer.coords = np.zeros((3, count_mol))
                    
                    if reader.id is not None:
                        if writer.id is None or len(writer.id) != count_mol:
                            writer.id = np.zeros(count_mol, dtype=int)
                    
                    if reader.type is not None:
                        if writer.type is None or len(writer.type) != count_mol:
                            writer.type = np.zeros(count_mol, dtype=int)
                    
                    if reader.mol is not None:
                        if writer.mol is None or len(writer.mol) != count_mol:
                            writer.mol = np.zeros(count_mol, dtype=int)
                    
                    if reader.image_flags is not None:
                        if writer.image_flags is None or writer.image_flags.shape[1] != count_mol:
                            writer.image_flags = np.zeros((3, count_mol), dtype=int)
                    
                    # データをコピー
                    new_idx = 0
                    for i in range(reader.nparticles):
                        if reader.mol[i] == args.mol_id:
                            writer.coords[:, new_idx] = reader.coords[:, i]
                            
                            if reader.id is not None and writer.id is not None:
                                writer.id[new_idx] = reader.id[i]
                            
                            if reader.type is not None and writer.type is not None:
                                writer.type[new_idx] = reader.type[i]
                            
                            if reader.mol is not None and writer.mol is not None:
                                writer.mol[new_idx] = reader.mol[i]
                            
                            if reader.image_flags is not None and writer.image_flags is not None:
                                writer.image_flags[:, new_idx] = reader.image_flags[:, i]
                            
                            new_idx += 1
                    
                    # 書き込み
                    writer.write()
                else:
                    print(f"mol_id={args.mol_id}の粒子がありません")
            else:
                # 全粒子を書き出し
                writer.timestep = reader.timestep
                writer.nparticles = reader.nparticles
                writer.box_bounds = reader.box_bounds
                writer.coords = reader.coords
                writer.id = reader.id
                writer.type = reader.type
                writer.mol = reader.mol
                writer.image_flags = reader.image_flags
                
                writer.write()
        
        nframes += 1
    
    # ファイルを閉じる
    reader.close()
    if writer:
        writer.close()
        
    te = time.time()
    print(f"処理時間: {te - ts:.2f}秒")
    print(f"フレーム数: {nframes}")
    print(f"1秒あたりの処理フレーム数: {nframes / (te - ts):.2f}")
