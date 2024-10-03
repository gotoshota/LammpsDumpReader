import argparse

def parse_args():
    """
    コマンドライン引数を解析する関数。
    必要な引数（trj1, trj2, output）を取得。
    """
    parser = argparse.ArgumentParser(description="2つのLAMMPSトラジェクトリーファイルを結合します。")
    parser.add_argument("--trj1", required=True, help="1つ目のLAMMPSトラジェクトリーファイルのパス。")
    parser.add_argument("--trj2", required=True, help="2つ目のLAMMPSトラジェクトリーファイルのパス。")
    parser.add_argument("--output", required=True, help="出力ファイルのパス。")
    parser.add_argument("--help", "-h", action="store_true", help="このヘルプメッセージを表示します。")
    return parser.parse_args()

def display_usage():
    """
    使用方法を表示する関数。
    """
    print("使い方:")
    print("  python merge_lammpstrj.py --trj1 <trj1のパス> --trj2 <trj2のパス> --output <出力ファイルのパス>")
    print("オプション:")
    print("  --trj1 <trj1のパス>  1つ目のLAMMPSトラジェクトリーファイルのパス。")
    print("  --trj2 <trj2のパス>  2つ目のLAMMPSトラジェクトリーファイルのパス。")
    print("  --output <出力ファイルのパス> 出力ファイルのパス。")
    print("  --help                 このヘルプメッセージを表示します。")

def find_timestep_position(file_path):
    """
    LAMMPSダンプファイルからタイムステップの位置を見つける関数。
    各タイムステップのファイル内の位置を辞書形式で返す。
    """
    timestep_positions = {}
    with open(file_path, 'r') as file:
        position = 0
        while True:
            line = file.readline()
            if not line:
                break
            if line.startswith("ITEM: TIMESTEP"):
                timestep = int(file.readline().strip())
                timestep_positions[timestep] = position
            position = file.tell()
    return timestep_positions

def find_merge_timestep(timesteps1, timesteps2):
    """
    2つのファイルのタイムステップを比較し、一致する最初のタイムステップを見つける関数。
    """
    for timestep in timesteps1:
        if timestep in timesteps2:
            return timestep
    return None

def write_merged_file(file1_path, file2_path, timesteps1, timesteps2, merge_timestep, output_path):
    """
    2つのファイルを結合し、出力ファイルに書き込む関数。
    """
    with open(file1_path, 'r') as file1, open(file2_path, 'r') as file2, open(output_path, 'w') as output_file:
        # file1からの最初の部分を出力ファイルに書き込む
        file1.seek(0)
        while file1.tell() < timesteps1[merge_timestep]:
            output_file.write(file1.readline())
        
        # file2からの後半部分を出力ファイルに書き込む
        file2.seek(timesteps2[merge_timestep])
        while True:
            line = file2.readline()
            if not line:
                break
            output_file.write(line)

def main():
    """
    メイン関数。
    コマンドライン引数を解析し、ファイルのタイムステップ位置を見つけ、ファイルを結合する。
    """
    args = parse_args()
    if args.help:
        display_usage()
        return
    
    # タイムステップの位置を見つける
    timesteps1 = find_timestep_position(args.trj1)
    timesteps2 = find_timestep_position(args.trj2)

    # 結合ポイントを見つける
    merge_timestep = find_merge_timestep(timesteps1, timesteps2)
    if merge_timestep is None:
        raise ValueError("一致するTIMESTEPが見つかりませんでした。")

    # 結合ファイルを書き出す
    write_merged_file(args.trj1, args.trj2, timesteps1, timesteps2, merge_timestep, args.output)

    print(f"結合されたファイルは '{args.output}' に保存されました。")

if __name__ == "__main__":
    main()

