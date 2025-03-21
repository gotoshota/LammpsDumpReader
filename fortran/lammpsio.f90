!> @file lammpsio.f90
!> @brief LAMMPSトラジェクトリファイルの読み込みと書き込みを行うモジュール
!> @details このモジュールはLAMMPSの分子動力学シミュレーションの出力ファイル（dump形式）を
!> 読み込み、処理、書き込みするための機能を提供します。
module lammpsio
    implicit none

    !> @brief 座標データを管理する型
    type :: coordinates
        real, allocatable :: coords(:, :)
        logical :: is_wrapped = .false. !< ラップされた座標を保持しているかのフラグ
        logical :: is_scaled = .false. !< スケーリングされた座標を保持しているかのフラグ
        logical :: has_image_flags = .false. !< イメージフラグを保持しているかのフラグ
    contains
        procedure :: wrap => wrap_coordinates
        procedure :: unwrap => unwrap_coordinates
    end type coordinates

    !> @brief 原子データのカラムインデックスを管理する構造体
    !> @details LAMMPSダンプファイル内の各カラム（id, x, y, z など）の位置情報を保持します
    type :: AtomIndex
        INTEGER(KIND=4) :: id = 0 !< 原子IDのカラムインデックス
        INTEGER(KIND=4) :: mol = 0 !< 分子IDのカラムインデックス
        INTEGER(KIND=4) :: type = 0 !< 原子タイプのカラムインデックス
        INTEGER(KIND=4) :: xu = 0 !< アンラップしたx座標のカラムインデックス
        INTEGER(KIND=4) :: yu = 0 !< アンラップしたy座標のカラムインデックス
        INTEGER(KIND=4) :: zu = 0 !< アンラップしたz座標のカラムインデックス
        INTEGER(KIND=4) :: x = 0 !< x座標のカラムインデックス
        INTEGER(KIND=4) :: y = 0 !< y座標のカラムインデックス
        INTEGER(KIND=4) :: z = 0 !< z座標のカラムインデックス
        INTEGER(KIND=4) :: ix = 0 !< xイメージフラグのカラムインデックス
        INTEGER(KIND=4) :: iy = 0 !< yイメージフラグのカラムインデックス
        INTEGER(KIND=4) :: iz = 0 !< zイメージフラグのカラムインデックス
        integer(KIND=4) :: xs = 0 !< スケーリングされたx座標のカラムインデックス
        integer(KIND=4) :: ys = 0 !< スケーリングされたy座標のカラムインデックス
        integer(KIND=4) :: zs = 0 !< スケーリングされたz座標のカラムインデックス
        integer(KIND=4) :: xsu = 0 !< スケーリングされたアンラップしたx座標のカラムインデックス
        integer(KIND=4) :: ysu = 0 !< スケーリングされたアンラップしたy座標のカラムインデックス
        integer(KIND=4) :: zsu = 0 !< スケーリングされたアンラップしたz座標のカラムインデックス
    end type AtomIndex

    !> @brief LAMMPSトラジェクトリデータを管理する基本構造体
    !> @details トラジェクトリファイルの読み込みと書き込み機能を提供する親クラス
    type :: lammpstrj
        character(len=256) :: filename !< ファイル名
        integer :: unit !< ファイルユニット番号
        logical :: end_of_file = .false. !< EOFフラグ
        logical :: has_atom_idx = .false. !< 原子インデックスが初期化されているかのフラグ
        logical :: is_writing = .false. !< 書き込みモードフラグ
        integer :: timestep !< 現在のタイムステップ
        integer :: nparticles !< 粒子数
        double precision :: box_bounds(3, 3) = 0.0d0 !< シミュレーションボックスの境界
        type(coordinates) :: coords
        integer, allocatable :: image_flags(:, :) !< イメージフラグ (3, nparticles)
        integer, allocatable :: id(:) !< 原子ID
        integer, allocatable :: mol(:) !< 分子ID
        integer, allocatable :: type(:) !< 原子タイプ
        type(AtomIndex) :: atom_idx !< カラムインデックス情報
    contains
        procedure :: open => open_trajectory !< ファイルを開く
        procedure :: close => close_trajectory !< ファイルを閉じる
        procedure :: read => read_next_frame !< 次のフレームを読み込む
        procedure :: write => write_lammpstrj !< フレームを書き込む
    end type lammpstrj

    !> @brief 読み込み専用のLAMMPSトラジェクトリ管理クラス
    !> @details lammpstrjを継承し、読み込み機能に特化した子クラス
    type, extends(lammpstrj) :: lammpstrjReader
    contains
        procedure :: open => open_reader !< 読み込み用にファイルを開く
    end type lammpstrjReader

    !> @brief 書き込み専用のLAMMPSトラジェクトリ管理クラス
    !> @details lammpstrjを継承し、書き込み機能に特化した子クラス
    type, extends(lammpstrj) :: lammpstrjWriter
    contains
        procedure :: open => open_writer !< 書き込み用にファイルを開く
        procedure :: create => create_trajectory !< 新規ファイルを作成する
        procedure :: append => append_trajectory !< 既存ファイルに追記する
    end type lammpstrjWriter

contains

    !> @brief トラジェクトリファイルを開く
    !> @param[in,out] this lammpstrjオブジェクト
    !> @param[in] filename 開くファイル名
    !> @param[in] mode オプションのファイルモード ('read', 'write', 'append')
    subroutine open_trajectory(this, filename, mode)
        class(lammpstrj), intent(inout) :: this
        character(len=*), intent(in) :: filename
        character(len=*), intent(in), optional :: mode ! 'read', 'write', 'append'
        integer :: ios
        character(len=10) :: use_mode

        this%filename = filename

        ! デフォルトは読み込みモード
        use_mode = 'read'
        if (present(mode)) use_mode = mode

        select case (trim(use_mode))
        case ('read')
            this%is_writing = .false.
            open (newunit=this%unit, file=filename, status='old', action='read', iostat=ios)
        case ('write')
            this%is_writing = .true.
            open (newunit=this%unit, file=filename, status='replace', action='write', iostat=ios)
        case ('append')
            this%is_writing = .true.
            open (newunit=this%unit, file=filename, status='unknown', position='append', action='write', iostat=ios)
        case default
            print *, "エラー: 不明なモードです: ", use_mode
            stop
        end select

        if (ios /= 0) then
            print *, "エラー: ファイルを開けません: ", filename
            stop
        end if
    end subroutine open_trajectory

    !> @brief 読み込み専用でトラジェクトリファイルを開く
    !> @param[in,out] this lammpstrjReaderオブジェクト
    !> @param[in] filename 開くファイル名
    !> @param[in] mode 無視されるパラメータ（互換性のために残されている）
    subroutine open_reader(this, filename, mode)
        class(lammpstrjReader), intent(inout) :: this
        character(len=*), intent(in) :: filename
        character(len=*), intent(in), optional :: mode

        ! モードを無視して読み込みモードで開く
        this%is_writing = .false.
        call open_trajectory(this, filename, 'read')
    end subroutine open_reader

    !> @brief 書き込み専用でトラジェクトリファイルを開く
    !> @param[in,out] this lammpstrjWriterオブジェクト
    !> @param[in] filename 開くファイル名
    !> @param[in] mode オプションのファイルモード ('write'または'append')
    subroutine open_writer(this, filename, mode)
        class(lammpstrjWriter), intent(inout) :: this
        character(len=*), intent(in) :: filename
        character(len=*), intent(in), optional :: mode
        character(len=10) :: use_mode

        ! デフォルトは新規書き込みモード
        use_mode = 'write'
        if (present(mode)) use_mode = mode

        if (trim(use_mode) == 'read') then
            print *, "警告: Writerは読み込みモードではオープンできません。書き込みモードに変更します。"
            use_mode = 'write'
        end if

        this%is_writing = .true.
        call open_trajectory(this, filename, use_mode)
    end subroutine open_writer

    !> @brief 新規書き込み用のトラジェクトリファイルを作成する
    !> @param[in,out] this lammpstrjWriterオブジェクト
    !> @param[in] filename 作成するファイル名
    subroutine create_trajectory(this, filename)
        class(lammpstrjWriter), intent(inout) :: this
        character(len=*), intent(in) :: filename

        call open_writer(this, filename, 'write')
    end subroutine create_trajectory

    !> @brief 既存のトラジェクトリファイルに追記するために開く
    !> @param[in,out] this lammpstrjWriterオブジェクト
    !> @param[in] filename 追記するファイル名
    subroutine append_trajectory(this, filename)
        class(lammpstrjWriter), intent(inout) :: this
        character(len=*), intent(in) :: filename

        call open_writer(this, filename, 'append')
    end subroutine append_trajectory

    !> @brief トラジェクトリファイルを閉じる
    !> @param[in,out] this lammpstrjオブジェクト
    subroutine close_trajectory(this)
        class(lammpstrj), intent(inout) :: this
        close (this%unit)
    end subroutine close_trajectory

    !> @brief トラジェクトリファイルから次のフレームを読み込む
    !> @param[in,out] this lammpstrjオブジェクト
    !> @details ファイルから次のタイムステップのデータを読み込み、構造体のフィールドを更新します
    subroutine read_next_frame(this)
        class(lammpstrj), intent(inout) :: this
        character(len=256) :: line
        CHARACTER(LEN=:), allocatable :: atom_header_parts(:), atom_parts(:), box_parts(:)
        integer :: ios, i, j
        integer :: nColums

        if (this%is_writing) then
            print *, "エラー: このトラジェクトリは書き込みモードのため読み込めません"
            return
        end if

        if (this%end_of_file) return

        read (this%unit, '(A)', iostat=ios) line ! ITEM: TIMESTEP
        if (ios /= 0) then
            this%end_of_file = .true.
            return
        end if
        read (this%unit, '(I10)', iostat=ios) this%timestep
        read (this%unit, '(A)', iostat=ios) line ! ITEM: NUMBER OF ATOMS
        read (this%unit, '(I10)', iostat=ios) this%nparticles
        read (this%unit, '(A)', iostat=ios) line ! ITEM: BOX BOUNDS
        read (this%unit, '(A)', iostat=ios) line ! xlo xhi (xy)
        box_parts = split_line(line)
        read (box_parts(:), *) this%box_bounds(1:size(box_parts), 1)
        read (this%unit, '(A)', iostat=ios) line ! ylo yhi (yz)
        box_parts = split_line(line)
        read (box_parts(:), *) this%box_bounds(1:size(box_parts), 2)
        read (this%unit, '(A)', iostat=ios) line ! zlo zhi (xz)
        box_parts = split_line(line)
        read (box_parts(:), *) this%box_bounds(1:size(box_parts), 3)
        if (this%has_atom_idx) then
            read (this%unit, '(A)', iostat=ios) line ! Skip Atom header line and use the same atom header
        else if (this%has_atom_idx .eqv. .false.) then
            read (this%unit, '(A)', iostat=ios) line ! Read Atom header line and aplit it
            atom_header_parts = split_line(line)
            nColums = size(atom_header_parts) - 2
            do i = 1, nColums ! 最初の "ITEM:" "ATOMS" はスキップ
                select case (trim(adjustl(atom_header_parts(i + 2))))
                case ("id")
                    this%atom_idx%id = i
                case ("mol")
                    this%atom_idx%mol = i
                case ("type")
                    this%atom_idx%type = i
                case ("xu")
                    this%atom_idx%xu = i
                case ("yu")
                    this%atom_idx%yu = i
                case ("zu")
                    this%atom_idx%zu = i
                case ("x")
                    this%atom_idx%x = i
                case ("y")
                    this%atom_idx%y = i
                case ("z")
                    this%atom_idx%z = i
                case ("ix")
                    this%atom_idx%ix = i
                case ("iy")
                    this%atom_idx%iy = i
                case ("iz")
                    this%atom_idx%iz = i
                case ("xs")
                    this%atom_idx%xs = i
                case ("ys")
                    this%atom_idx%ys = i
                case ("zs")
                    this%atom_idx%zs = i
                case ("xsu")
                    this%atom_idx%xsu = i
                case ("ysu")
                    this%atom_idx%ysu = i
                case ("zsu")
                    this%atom_idx%zsu = i
                case default
                    print *, "Error: Unknown column to read header"
                    print *, "Column: ", i
                end select
            end do
            this%has_atom_idx = .true.
            if (.not. allocated(this%coords%coords)) allocate (this%coords%coords(3, this%nparticles)) ! x, y, z
            if (this%atom_idx%id > 0) allocate (this%id(this%nparticles))
            if (this%atom_idx%mol > 0) allocate (this%mol(this%nparticles))
            if (this%atom_idx%type > 0) allocate (this%type(this%nparticles))
            if (this%atom_idx%ix > 0) then
                allocate (this%image_flags(3, this%nparticles))
                this%coords%has_image_flags = .true.
            end if
        end if

        if (.not. allocated(this%coords%coords)) then
            allocate (this%coords%coords(3, this%nparticles))
        end if

        do i = 1, this%nparticles
            read (this%unit, "(A)") line
            atom_parts = split_line(adjustl(line))
            do j = 1, size(atom_parts)
                if (j == this%atom_idx%x) then
                    read (atom_parts(j), *) this%coords%coords(1, i)
                else if (j == this%atom_idx%y) then
                    read (atom_parts(j), *) this%coords%coords(2, i)
                else if (j == this%atom_idx%z) then
                    read (atom_parts(j), *) this%coords%coords(3, i)
                else if (j == this%atom_idx%xs) then
                    read (atom_parts(j), *) this%coords%coords(1, i)
                else if (j == this%atom_idx%ys) then
                    read (atom_parts(j), *) this%coords%coords(2, i)
                else if (j == this%atom_idx%zs) then
                    read (atom_parts(j), *) this%coords%coords(3, i)
                else if (j == this%atom_idx%xu) then
                    read (atom_parts(j), *) this%coords%coords(1, i)
                else if (j == this%atom_idx%yu) then
                    read (atom_parts(j), *) this%coords%coords(2, i)
                else if (j == this%atom_idx%zu) then
                    read (atom_parts(j), *) this%coords%coords(3, i)
                else if (j == this%atom_idx%id) then
                    read (atom_parts(j), *) this%id(i)
                else if (j == this%atom_idx%mol) then
                    read (atom_parts(j), *) this%mol(i)
                else if (j == this%atom_idx%type) then
                    read (atom_parts(j), *) this%type(i)
                else if (j == this%atom_idx%ix) then
                    read (atom_parts(j), *) this%image_flags(1, i)
                else if (j == this%atom_idx%iy) then
                    read (atom_parts(j), *) this%image_flags(2, i)
                else if (j == this%atom_idx%iz) then
                    read (atom_parts(j), *) this%image_flags(3, i)
                else
                    print *, "Error: Unknown column to substitute variables"
                    print *, "Column: ", j
                    print *, this%atom_idx
                end if
            end do
        end do
    end subroutine read_next_frame

    !> @brief 文字列を空白で区切って配列に分割する
    !> @param[in] str 分割する文字列
    !> @return 分割された文字列の配列
    !> @details 連続した空白は一つの区切りとして扱います
    function split_line(str) result(substrings)
        character(len=*), intent(in) :: str
        character(len=:), allocatable :: substrings(:)
        integer :: i, start, word_end, nWords, strLen
        logical :: inWord

        strLen = len_trim(str)
        nWords = 0
        inWord = .false.
        start = 1

        ! 文字列を走査し、単語の数を数える
        do i = 1, strLen
            if (str(i:i) /= ' ' .and. .not. inWord) then
                inWord = .true.
                start = i
            elseif (str(i:i) == ' ' .and. inWord) then
                inWord = .false.
                nWords = nWords + 1
            end if
        end do
        if (inWord) nWords = nWords + 1 ! 最後の単語を数える

        allocate (character(len=strLen) :: substrings(nWords))
        nWords = 0
        inWord = .false.

        ! 実際に単語を抽出する
        do i = 1, strLen
            if (str(i:i) /= ' ' .and. .not. inWord) then
                inWord = .true.
                start = i
            elseif (str(i:i) == ' ' .and. inWord) then
                inWord = .false.
                word_end = i - 1
                nWords = nWords + 1
                substrings(nWords) = str(start:word_end)
            end if
        end do
        if (inWord) then
            nWords = nWords + 1
            substrings(nWords) = str(start:strLen)
        end if
    end function split_line

    !> @brief トラジェクトリデータをファイルに書き込む
    !> @param[in] this lammpstrjオブジェクト
    !> @details 現在のフレームデータをLAMMPS形式でファイルに書き込みます
    subroutine write_lammpstrj(this)
        class(lammpstrj), intent(in) :: this
        integer :: i

        if (.not. this%is_writing) then
            print *, "エラー: このトラジェクトリは読み込みモードのため書き込めません"
            return
        end if

        ! タイムステップ情報を書き込み
        write (this%unit, '(A)') "ITEM: TIMESTEP"
        write (this%unit, '(I10)') this%timestep

        ! 粒子数を書き込み
        write (this%unit, '(A)') "ITEM: NUMBER OF ATOMS"
        write (this%unit, '(I10)') this%nparticles

        ! ボックス境界情報を書き込み
        if (ubound(this%box_bounds, 1) == 2) then
            write (this%unit, '(A)') "ITEM: BOX BOUNDS pp pp pp"
            write (this%unit, '(2(E20.12,1x))') this%box_bounds(:, 1) ! xlo xhi
            write (this%unit, '(2(E20.12,1x))') this%box_bounds(:, 2) ! ylo yhi
            write (this%unit, '(2(E20.12,1x))') this%box_bounds(:, 3) ! zlo zhi
        else
            write (this%unit, '(A)') "ITEM: BOX BOUNDS xy xz yz pp pp pp"
            write (this%unit, '(3E20.12)') this%box_bounds(:, 1) ! xlo xhi xy
            write (this%unit, '(3E20.12)') this%box_bounds(:, 2) ! ylo yhi xz
            write (this%unit, '(3E20.12)') this%box_bounds(:, 3) ! zlo zhi yz
        end if

        ! 原子データを書き込み
        write (this%unit, '(A)', advance='no') "ITEM: ATOMS id"

        ! ヘッダー情報の書き込み
        if (allocated(this%type)) write (this%unit, '(A)', advance='no') " type"
        if (allocated(this%mol)) write (this%unit, '(A)', advance='no') " mol"
        write (this%unit, '(A)', advance='no') " x y z"
        if (allocated(this%image_flags)) write (this%unit, '(A)', advance='no') " ix iy iz"
        write (this%unit, '(A)') ! 改行

        ! 各原子のデータを書き込み
        do i = 1, this%nparticles
            ! ID を書き込み (必須)
            if (allocated(this%id)) then
                write (this%unit, '(I10)', advance='no') this%id(i)
            else
                write (this%unit, '(I10)', advance='no') i
            end if

            ! タイプを書き込み (オプション)
            if (allocated(this%type)) then
                write (this%unit, '(I5)', advance='no') this%type(i)
            end if

            ! 分子 ID を書き込み (オプション)
            if (allocated(this%mol)) then
                write (this%unit, '(I5)', advance='no') this%mol(i)
            end if

            ! 座標を書き込み (必須)
            write (this%unit, '(3F15.7)', advance='no') this%coords%coords(:, i)

            ! イメージフラグを書き込み (オプション)
            if (allocated(this%image_flags)) then
                write (this%unit, '(3I5)', advance='no') this%image_flags(:, i)
            end if

            write (this%unit, '(A)') ! 改行
        end do
    end subroutine write_lammpstrj

!> @brief 座標データをwrapするサブルーチン
!> @details 傾いたボックスにも対応。イメージフラグを計算または使用します。
!> @param[in] this coordinates型のインスタンス
!> @param[in,out] parent lammpstrj型のインスタンス
!> @return wrapped ラップされた座標
    function wrap_coordinates(this, parent) result(wrapped)
        class(coordinates), intent(in) :: this
        class(lammpstrj), intent(inout) :: parent
        real, allocatable :: fractional_coords(:, :)
        real, allocatable :: wrapped(:, :)
        integer(kind=8), allocatable :: image_flags(:, :)
        integer :: i, j, np
        double precision :: box_len(3), center(3), fractional_center(3)
        double precision :: triclinic_box_len(3), triclinic_center(3)
        double precision :: disp(3)
        double precision :: box_matrix(3, 3)
        double precision :: box_matrix_inv(3, 3)
        logical :: is_triclinic

        call get_box_info(parent%box_bounds, box_len, center)
        call get_box_matrix(parent%box_bounds, box_matrix, box_matrix_inv)
        ! triclinic boxかどうかを判定（box_boundsの3列目に値があるかどうか）
        is_triclinic = .false.
        if (abs(parent%box_bounds(1, 3)) > 1.0e-8 .or. &
            abs(parent%box_bounds(2, 3)) > 1.0e-8 .or. &
            abs(parent%box_bounds(3, 3)) > 1.0e-8) then
            is_triclinic = .true.
        end if

        if (.not. allocated(this%coords)) then
            print *, "Error: coords data is not allocated."
            stop 1
        end if
        np = size(this%coords, 2)
        allocate (wrapped(3, np))
        if (allocated(parent%image_flags)) then
            image_flags = parent%image_flags
        else
            allocate (image_flags(3, np))
            allocate (parent%image_flags(3, np))
        end if

        ! 傾いたボックスの場合
        if (is_triclinic) then
            allocate (fractional_coords(3, np))
            if (.not. this%has_image_flags) then
                ! fractional coordinates
                do i = 1, np
                    fractional_coords(:, i) = matmul(box_matrix_inv, this%coords(:, i) - center)
                end do

                ! wrap fractional coordinates
                image_flags(3, :) = nint(fractional_coords(3, :))
                wrapped(3, :) = fractional_coords(3, :) - image_flags(3, :)
                image_flags(2, :) = nint(fractional_coords(2, :) - image_flags(3, :)*parent%box_bounds(3, 3))
                wrapped(2, :) = fractional_coords(2, :) - image_flags(2, :) - image_flags(3, :)*parent%box_bounds(3, 3)
                image_flags(1, :) = nint(fractional_coords(1, :) - image_flags(3, :)*parent%box_bounds(3, 2) - &
                                          image_flags(2, :)*parent%box_bounds(3, 1))
                wrapped(1, :) = fractional_coords(1, :) - image_flags(1, :) - image_flags(3, :)*parent%box_bounds(3, 2) - &
                                image_flags(2, :)*parent%box_bounds(3, 1)
                parent%image_flags = image_flags

                ! cartesian coordinates
                wrapped = matmul(box_matrix, wrapped)
                do i = 1, np
                    wrapped(:, i) = wrapped(:, i) + center
                end do
            else
                wrapped(3, :) = this%coords(3, :) - box_len(3)*parent%image_flags(3, :)
                wrapped(2, :) = this%coords(2, :) - box_len(2)*parent%image_flags(2, :) - parent%image_flags(3, :)*parent%box_bounds(3, 3)
                wrapped(1, :) = this%coords(1, :) - box_len(1)*parent%image_flags(1, :) - parent%image_flags(3, :)*parent%box_bounds(3, 2) - &
                              parent%image_flags(2, :)*parent%box_bounds(3, 1)
            end if
        else
            ! 直交ボックスの場合（既存のコード）
            if (.not. allocated(parent%image_flags)) then
                allocate (parent%image_flags(3, np))
                do j = 1, np
                    disp = this%coords(:, j) - center
                    do i = 1, 3
                        parent%image_flags(i, j) = nint(disp(i)/box_len(i))
                        wrapped(i, j) = this%coords(i, j) - box_len(i)*parent%image_flags(i, j)
                    end do
                end do
            end if

            do j = 1, np
                do i = 1, 3
                    wrapped(i, j) = this%coords(i, j) - box_len(i)*parent%image_flags(i, j)
                end do
            end do
        end if
    end function wrap_coordinates

!> @brief 座標データをunwrapするサブルーチン
!> @details イメージフラグを使用して周期境界条件を取り除きます。傾いたボックスにも対応。k
!> @param[in] this coordinates型のインスタンス
!> @param[in,out] parent lammpstrj型のインスタンス
!> @return unwrapped アンラップされた座標
    function unwrap_coordinates(this, parent) result(unwrapped)
        class(coordinates), intent(in) :: this
        class(lammpstrj), intent(inout) :: parent
        real, allocatable :: unwrapped(:, :)
        real, allocatable :: fractional_coords(:, :)
        double precision :: box_len(3), center(3)
        double precision :: box_matrix(3, 3), box_matrix_inv(3, 3)

        integer :: i, j, np
        logical :: is_triclinic

        call get_box_matrix(parent%box_bounds, box_matrix, box_matrix_inv)
        call get_box_info(parent%box_bounds, box_len, center)
        ! triclinic boxかどうかを判定（box_boundsの3列目に値があるかどうか）
        is_triclinic = .false.
        if (abs(parent%box_bounds(1, 3)) > 1.0e-8 .or. &
            abs(parent%box_bounds(2, 3)) > 1.0e-8 .or. &
            abs(parent%box_bounds(3, 3)) > 1.0e-8) then
            is_triclinic = .true.
        end if

        if (.not. allocated(parent%image_flags)) then
            print *, "Error: image_flags が定義されていません。unwrapできません。"
            stop 1
        end if
        np = size(this%coords, 2)
        allocate (unwrapped(3, np))

        ! 傾いたボックスの場合
        if (is_triclinic) then
            ! fractional coordinates
            allocate (fractional_coords(3, np))
            do i = 1, np
                fractional_coords(:, i) = matmul(box_matrix_inv, this%coords(:, i) - center)
            end do
            unwrapped(3, :) = fractional_coords(3, :) + parent%image_flags(3, :)
            unwrapped(2, :) = fractional_coords(2, :) + parent%image_flags(2, :) + parent%image_flags(3, :)*parent%box_bounds(3, 3)
            unwrapped(1, :) = fractional_coords(1, :) + parent%image_flags(1, :) + parent%image_flags(3, :)*parent%box_bounds(3, 2) + &
                          parent%image_flags(2, :)*parent%box_bounds(3, 1)
            unwrapped = matmul(box_matrix, unwrapped)
            do i = 1, np
                unwrapped(:, i) = unwrapped(:, i) + center
            end do
        else
            ! 直交ボックスの場合（既存のコード）
            do j = 1, np
                do i = 1, 3
                    unwrapped(i, j) = this%coords(i, j) + parent%image_flags(i, j)*box_len(i)
                end do
            end do
        end if
    end function unwrap_coordinates

    subroutine get_box_info(box_bounds, box_size, center)
        implicit none

        double precision, intent(in) :: box_bounds(:, :) ! LAMMPSのbox_bounds
        ! box_bounds は (1:3, 1:3) の配列で、
        ! xlo, xhi, xy
        ! ylo, yhi, xz
        ! zlo, zhi, yz
        ! の順で格納されている
        double precision, intent(out) :: box_size(3)
        double precision, intent(out) :: center(3)
        double precision :: xlo, xhi, ylo, yhi, zlo, zhi ! 実際のボックスの境界

        ! LAMMPSのbox_boundsから実際のボックスの境界を計算
        xlo = box_bounds(1, 1) - min(0.0, box_bounds(3, 1), box_bounds(3, 2), box_bounds(3, 1) + box_bounds(3, 2))
        xhi = box_bounds(2, 1) - max(0.0, box_bounds(3, 1), box_bounds(3, 2), box_bounds(3, 1) + box_bounds(3, 2))
        ylo = box_bounds(1, 2) - min(0.0, box_bounds(3, 3))
        yhi = box_bounds(2, 2) - max(0.0, box_bounds(3, 3))
        zlo = box_bounds(1, 3)
        zhi = box_bounds(2, 3)
        ! Simulation cell のサイズと中心を計算
        box_size(1) = xhi - xlo
        box_size(2) = yhi - ylo
        box_size(3) = zhi - zlo
        center(1) = (box_bounds(1, 1) + box_bounds(2, 1))/2.0
        center(2) = (box_bounds(1, 2) + box_bounds(2, 2))/2.0
        center(3) = (box_bounds(1, 3) + box_bounds(2, 3))/2.0
    end subroutine get_box_info

    subroutine get_box_matrix(box_bounds, A, A_inv)
        implicit none

        double precision, intent(in) :: box_bounds(:, :) ! LAMMPSのbox_bounds
        double precision, intent(out) :: A(3, 3)
        double precision, intent(out) :: A_inv(3, 3)

        double precision :: xlo, xhi, ylo, yhi, zlo, zhi ! 実際のボックスの境界
        xlo = box_bounds(1, 1) - min(0.0, box_bounds(3, 1), box_bounds(3, 2), box_bounds(3, 1) + box_bounds(3, 2))
        xhi = box_bounds(2, 1) - max(0.0, box_bounds(3, 1), box_bounds(3, 2), box_bounds(3, 1) + box_bounds(3, 2))
        ylo = box_bounds(1, 2) - min(0.0, box_bounds(3, 3))
        yhi = box_bounds(2, 2) - max(0.0, box_bounds(3, 3))
        zlo = box_bounds(1, 3)
        zhi = box_bounds(2, 3)

        A(1, 1) = xhi - xlo
        A(2, 2) = yhi - ylo
        A(3, 3) = zhi - zlo
        A(1, 2) = box_bounds(3, 1)
        A(1, 3) = box_bounds(3, 2)
        A(2, 3) = box_bounds(3, 3)
        A(2, 1) = 0.0
        A(3, 1) = 0.0
        A(3, 2) = 0.0

        A_inv(1, 1) = 1.0/A(1, 1)
        A_inv(2, 2) = 1.0/A(2, 2)
        A_inv(3, 3) = 1.0/A(3, 3)
        A_inv(1, 2) = -A(1, 2)/(A(1, 1)*A(2, 2))
        A_inv(1, 3) = (A(1, 2)*A(2, 3) - A(1, 3)*A(2, 2))/(A(1, 1)*A(2, 2)*A(3, 3))
        A_inv(2, 3) = -A(2, 3)/(A(2, 2)*A(3, 3))
        A_inv(2, 1) = 0.0
        A_inv(3, 1) = 0.0
        A_inv(3, 2) = 0.0
    end subroutine get_box_matrix

end module lammpsio
