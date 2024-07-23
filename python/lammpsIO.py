import numpy as np

class lammpstrjReader:
    def __init__(self, filename):
        self.filename = filename
        self.file = open(filename, 'r')
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
            #self.box_bounds = np.array([self.file.readline().strip().split() for _ in range(3)], dtype=float)
            self.box_bounds = np.array([list(map(float, self.file.readline().strip().split())) for _ in range(3)])
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
                    raise ValueError(f"Expected {self.num_columns} columns, but got {len(atom_data)}")
                
                # 各カラムに基づいて適切な型に変換
                parsed_data = []
                for label, value in zip(self.atom_data_labels, atom_data):
                    if label in ['id', 'mol']:  # ここでは'id'や'mol'は整数として処理
                        parsed_data.append(int(value))
                    else:  # その他のデータは浮動小数点数として処理
                        parsed_data.append(float(value))

                parsed_data = np.array(parsed_data)  # NumPy配列に変換

                # カラム名に基づいてデータを格納
                if 'x' in self.atom_data_labels:
                    x_index = self.atom_data_labels.index('x')
                    self.coords[i, 0] = atom_data[x_index]
                if 'y' in self.atom_data_labels:
                    y_index = self.atom_data_labels.index('y')
                    self.coords[i, 1] = atom_data[y_index]
                if 'z' in self.atom_data_labels:
                    z_index = self.atom_data_labels.index('z')
                    self.coords[i, 2] = atom_data[z_index]
                if 'xu' in self.atom_data_labels:
                    xu_index = self.atom_data_labels.index('xu')
                    self.coords[i, 0] = atom_data[xu_index]
                if 'yu' in self.atom_data_labels:
                    yu_index = self.atom_data_labels.index('yu')
                    self.coords[i, 1] = atom_data[yu_index]
                if 'zu' in self.atom_data_labels:
                    zu_index = self.atom_data_labels.index('zu')
                    self.coords[i, 2] = atom_data[zu_index]
                if 'ix' in self.atom_data_labels:
                    ix_index = self.atom_data_labels.index('ix')
                    self.image_flags[i, 0] = int(atom_data[ix_index])
                if 'iy' in self.atom_data_labels:
                    iy_index = self.atom_data_labels.index('iy')
                    self.image_flags[i, 1] = int(atom_data[iy_index])
                if 'iz' in self.atom_data_labels:
                    iz_index = self.atom_data_labels.index('iz')
                    self.image_flags[i, 2] = int(atom_data[iz_index])
                if 'mol' in self.atom_data_labels:
                    mol_index = self.atom_data_labels.index('mol')
                    self.mol = int(atom_data[mol_index])
                if 'type' in self.atom_data_labels:
                    type_index = self.atom_data_labels.index('type')
                    self.type = int(atom_data[type_index])
                

            return {
                'timestep': self.timestep,
                'num_atoms': self.num_atoms,
                'box_bounds': self.box_bounds,
                'coords': self.coords,
                'image_flags': self.image_flags
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

if __name__ == '__main__':
    import argparse
    import time
    # 使用例
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, required=True)
    
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

