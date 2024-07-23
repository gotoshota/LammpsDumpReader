import numpy as np
import argparse

class lammpstrjReader:
    def __init__(self, filename):
        self.filename = filename
        self.file = open(filename, 'r')
        self.end_of_file = False
        self.num_columns = None

    def __del__(self):
        self.file.close()

    def read_next_frame(self):
        if self.end_of_file:
            return None

        line = self.file.readline()
        if not line:
            self.end_of_file = True
            return None

        # 読み飛ばす
        if not line.startswith('ITEM: TIMESTEP'):
            self.end_of_file = True
            return None

        timestep = int(self.file.readline().strip())
        self.file.readline()  # ITEM: NUMBER OF ATOMS
        num_atoms = int(self.file.readline().strip())
        self.file.readline()  # ITEM: BOX BOUNDS
        box_bounds = np.array([self.file.readline().strip().split() for _ in range(3)], dtype=float)
        self.file.readline()  # ITEM: ATOMS

        # 次の行を一時的に読み込み、列数を取得する
        if not self.num_columns:
            atom_data = self.file.readline().strip().split()
            num_columns = len(atom_data)

        # 先程読み込んだ行のデータをリストに格納
        atoms = np.zeros((num_atoms, num_columns))
        atoms[0, :] = list(map(float, atom_data))

        # 残りの行を読み込む
        for i in range(1, num_atoms):
            atom_data = list(map(float, self.file.readline().strip().split()))
            atoms[i, :] = atom_data

        return {
            'timestep': timestep,
            'num_atoms': num_atoms,
            'box_bounds': box_bounds,
            'atoms': atoms
        }

def write_lammpstrj(filename, coords, box_bounds):
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
        f.write("ITEM: ATOMS id xu yu zu\n")
        for i, coord in enumerate(coords):
            f.write(f"{i} {coord[0]} {coord[1]} {coord[2]}\n")

if __name__ == '__main__':
    # 使用例
    args = argparse.ArgumentParser()
    args.add_argument('--input', type=str, required=True)
    
    args = args.parse
    traj = lammpstrjReader(args.input)
    
    while True:
        frame = traj.read_next_frame()
        if frame is None:
            break
        # フレームごとの処理
        print(f"Timestep: {frame['timestep']}, Number of atoms: {frame['num_atoms']}")
        # "atoms" には行の全てのデータが入っているため、インデックスを指定する。
        idx_x = 3 # 例えばx座標のインデックス
        print(f"ATOMS: {frame['atoms'][1][idx_x]}")
        print(f"{len(frame['atoms'])}")


