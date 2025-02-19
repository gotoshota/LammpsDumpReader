import numpy as np


class lammpstrjReader:
    def __init__(self, filename):
        self.filename = filename
        self.file = open(filename, "r")
        self.end_of_file = False
        self.num_columns = None
        self.timestep = None
        self.num_atoms = None
        self.box_bounds = None
        self.atom_data_labels = None
        self.coords = None
        self.image_flags = None
        self.mol = None
        self.type = None

    def __del__(self):
        if self.file:
            self.file.close()

    def read_next_frame(self):
        if self.end_of_file:
            return None

        try:
            line = self.file.readline()
            if not line:
                self.end_of_file = True
                return None
            self.timestep = int(self.file.readline().strip())
            self.file.readline()  # ITEM: NUMBER OF ATOMS
            self.num_atoms = int(self.file.readline().strip())
            self.file.readline()  # ITEM: BOX BOUNDS
            # self.box_bounds = np.array([self.file.readline().strip().split() for _ in range(3)], dtype=float)
            self.box_bounds = np.array(
                [
                    list(map(float, self.file.readline().strip().split()))
                    for _ in range(3)
                ]
            )
            line = self.file.readline()  # ITEM: ATOMS

            # ITEM: ATOMS 行からカラム名を読み込み
            labels_line = line.strip()
            self.atom_data_labels = labels_line.split()
            # ITEM:, ATOMS は除く
            self.atom_data_labels = self.atom_data_labels[2:]
            self.num_columns = len(self.atom_data_labels)

            # データの初期化
            self.coords = np.zeros((self.num_atoms, 3))
            self.image_flags = np.zeros((self.num_atoms, 3), dtype=int)
            # データの読み込み
            for i in range(self.num_atoms):
                atom_data = self.file.readline().strip().split()

                if len(atom_data) != self.num_columns:
                    raise ValueError(
                        f"Expected {self.num_columns} columns, but got {len(atom_data)}"
                    )

                # 各カラムに基づいて適切な型に変換
                parsed_data = []
                for label, value in zip(self.atom_data_labels, atom_data):
                    if label in ["id", "mol"]:  # ここでは'id'や'mol'は整数として処理
                        parsed_data.append(int(value))
                    else:  # その他のデータは浮動小数点数として処理
                        parsed_data.append(float(value))

                parsed_data = np.array(parsed_data)  # NumPy配列に変換

                # カラム名に基づいてデータを格納
                if "x" in self.atom_data_labels:
                    x_index = self.atom_data_labels.index("x")
                    self.coords[i, 0] = atom_data[x_index]
                if "y" in self.atom_data_labels:
                    y_index = self.atom_data_labels.index("y")
                    self.coords[i, 1] = atom_data[y_index]
                if "z" in self.atom_data_labels:
                    z_index = self.atom_data_labels.index("z")
                    self.coords[i, 2] = atom_data[z_index]
                if "xu" in self.atom_data_labels:
                    xu_index = self.atom_data_labels.index("xu")
                    self.coords[i, 0] = atom_data[xu_index]
                if "yu" in self.atom_data_labels:
                    yu_index = self.atom_data_labels.index("yu")
                    self.coords[i, 1] = atom_data[yu_index]
                if "zu" in self.atom_data_labels:
                    zu_index = self.atom_data_labels.index("zu")
                    self.coords[i, 2] = atom_data[zu_index]
                if "ix" in self.atom_data_labels:
                    ix_index = self.atom_data_labels.index("ix")
                    self.image_flags[i, 0] = int(atom_data[ix_index])
                if "iy" in self.atom_data_labels:
                    iy_index = self.atom_data_labels.index("iy")
                    self.image_flags[i, 1] = int(atom_data[iy_index])
                if "iz" in self.atom_data_labels:
                    iz_index = self.atom_data_labels.index("iz")
                    self.image_flags[i, 2] = int(atom_data[iz_index])
                if "mol" in self.atom_data_labels:
                    mol_index = self.atom_data_labels.index("mol")
                    self.mol = int(atom_data[mol_index])
                if "type" in self.atom_data_labels:
                    type_index = self.atom_data_labels.index("type")
                    self.type = int(atom_data[type_index])

            return {
                "timestep": self.timestep,
                "num_atoms": self.num_atoms,
                "box_bounds": self.box_bounds,
                "coords": self.coords,
                "image_flags": self.image_flags,
            }
        except Exception as e:
            print(f"Error reading file: {e}")
            self.end_of_file = True
            return None


def write_lammpstrj(filename, coords, box_bounds, image_flags=None):
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
                    )
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

    args = parser.parse_args()
    lmp = lammpstrjReader(args.input)

    ts = time.time()
    nframes = 0
    while True:
        lmp.read_next_frame()
        if lmp.end_of_file:
            break
        # フレームごとの処理
        print(f"Timestep: {lmp.timestep}, Number of atoms: {lmp.num_atoms}")
        if lmp.coords is not None:
            print(f"coords: {lmp.coords[1, :]}")
            print(f"image_flags: {lmp.image_flags[1, :]}")
        nframes += 1
    te = time.time()
    print(f"Time: {te - ts}")
    print(f"Number of frames: {nframes}")
    print(f"read frames per second: {nframes / (te - ts)}")
