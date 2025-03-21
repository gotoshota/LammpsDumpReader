# コンパイラとフラグの設定
FC = gfortran
FFLAGS = -O3 -Wall

# ディレクトリの設定
SRCDIR = ../fortran
TESTDIR = .
DATADIR = ../data

# ソースファイルとオブジェクトファイル
SRC = $(SRCDIR)/lammpsio.f90

# テストプログラム
TEST_SRC = $(TESTDIR)/test.f90
TEST_EXE = test_fort

# デフォルトターゲット
all: $(TEST_EXE)

# テストプログラムのビルド
$(TEST_EXE): $(TEST_SRC)
	$(FC) $(FFLAGS) -o $(TEST_EXE) $(SRC) $(TEST_SRC)

# クリーンアップ
clean:
	rm -f $(OBJ) $(MOD) $(TEST_EXE)
	rm -f $(SRCDIR)/*.mod $(TESTDIR)/*.mod $(SRCDIR)/*.o

# ファイルの実行
run: $(TEST_EXE)
	./$(TEST_EXE) $(DATADIR)/dump.lammpstrj

.PHONY: all clean run

