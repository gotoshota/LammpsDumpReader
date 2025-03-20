!> @file main.f90
!> @brief LAMMPSトラジェクトリの読み込みと書き込み機能のテストプログラム
!> @details lammpstrjReaderとlammpstrjWriterクラスを使用して、LAMMPSトラジェクトリファイルの
!> 読み込みと書き込み機能をテストします。パフォーマンス測定も行います。
program test_lammpstrj
    use lammpsio
    implicit none

    type(lammpstrjReader) :: reader !< トラジェクトリ読み込み用オブジェクト
    type(lammpstrjWriter) :: writer !< トラジェクトリ書き込み用オブジェクト
    character(len=256) :: input_filename, output_filename
    integer :: i, frame
    integer :: start_time, end_time, read_time, write_time, clock_rate
    logical :: test_write = .true. !< 書き込みテストを行うフラグ
    integer :: count_mol1, j, new_idx !< mol_id = 1のフィルタリングのための変数
    real, allocatable :: wrapped(:,:) !< ラップされた座標用の一時配列

    integer :: debug_i = 0
    
    ! テスト用のファイル名を指定します。適宜変更してください。
    input_filename = "../data/triclinic.lammpstrj"
    output_filename = "../data/output.lammpstrj"

    ! 読み込み用ファイルを開く
    call reader%open(input_filename)
    
    ! 書き込み用ファイルを開く（テスト実行する場合）
    if (test_write) then
        call writer%create(output_filename)
    end if

    frame = 0
    read_time = 0
    write_time = 0
    call system_clock(count_rate=clock_rate)

    ! EOFまで読み込む
    do
        ! 読み込み時間を計測
        call system_clock(start_time)
        call reader%read()
        call system_clock(end_time)
        read_time = read_time + (end_time - start_time)

        if (reader%end_of_file) exit

        frame = frame + 1

        ! 基本情報を表示
        print *, "===== Frame: ", frame, " ====="
        print *, "Timestep: ", reader%timestep
        print *, "Number of particles: ", reader%nparticles
        print *, "Box bounds: "
        do i = 1, 3
            print *, reader%box_bounds(:, i)
        end do

        ! サンプルとして10番目の粒子の座標を表示
        if (reader%nparticles >= 10) then
            print *, "Coords of particle 10: ", reader%coords%data(:, 10)
            if (allocated(reader%id)) print *, "ID of particle 10: ", reader%id(10)
            if (allocated(reader%type)) print *, "Type of particle 10: ", reader%type(10)
        end if

        ! 書き込みテスト
        if (test_write) then
            ! mol_id = 1の粒子数をカウント
            count_mol1 = 0
            
            if (allocated(reader%mol)) then
                do j = 1, reader%nparticles
                    if (reader%mol(j) == 1) count_mol1 = count_mol1 + 1
                end do
            else
                print *, "警告: mol配列が割り当てられていません。全粒子を書き出します。"
                count_mol1 = reader%nparticles
            end if
            
            ! reader からデータをコピー
            writer%timestep = reader%timestep
            writer%nparticles = count_mol1  ! mol_id = 1の粒子数のみ
            writer%box_bounds = reader%box_bounds
            
            ! 配列データのコピー（mol_id = 1の粒子のみ）
            if (.not. allocated(writer%coords%data)) then
                allocate(writer%coords%data(3, count_mol1))
            else if (size(writer%coords%data, 2) /= count_mol1) then
                deallocate(writer%coords%data)
                allocate(writer%coords%data(3, count_mol1))
            end if
            
            ! その他の配列も必要に応じて再割り当て
            if (allocated(reader%id)) then
                if (.not. allocated(writer%id) .or. size(writer%id) /= count_mol1) then
                    if (allocated(writer%id)) deallocate(writer%id)
                    allocate(writer%id(count_mol1))
                end if
            end if
            
            if (allocated(reader%type)) then
                if (.not. allocated(writer%type) .or. size(writer%type) /= count_mol1) then
                    if (allocated(writer%type)) deallocate(writer%type)
                    allocate(writer%type(count_mol1))
                end if
            end if
            
            if (allocated(reader%mol)) then
                if (.not. allocated(writer%mol) .or. size(writer%mol) /= count_mol1) then
                    if (allocated(writer%mol)) deallocate(writer%mol)
                    allocate(writer%mol(count_mol1))
                end if
            end if
            
            if (allocated(reader%image_flags)) then
                if (.not. allocated(writer%image_flags) .or. size(writer%image_flags, 2) /= count_mol1) then
                    if (allocated(writer%image_flags)) deallocate(writer%image_flags)
                    allocate(writer%image_flags(3, count_mol1))
                end if
            end if
            
            ! mol_id = 1の粒子データのみコピー
            new_idx = 0
            do j = 1, reader%nparticles
                if (.not. allocated(reader%mol) .or. reader%mol(j) == 1) then
                    new_idx = new_idx + 1
                    !writer%coords%data(:, new_idx) = reader%coords%data(:, j)
                    ! Use the wrapped coordinates
                    if (.not. allocated(wrapped)) then
                        allocate(wrapped(3, reader%nparticles))
                    end if
                    wrapped = reader%coords%wrap(reader)
                    writer%coords%data(:, new_idx) = wrapped(:, j)
                    if (.not. allocated(writer%image_flags)) then
                        allocate(writer%image_flags(3, count_mol1))
                    end if
                    
                    if (allocated(reader%id)) writer%id(new_idx) = reader%id(j)
                    if (allocated(reader%type)) writer%type(new_idx) = reader%type(j)
                    if (allocated(reader%mol)) writer%mol(new_idx) = reader%mol(j)
                    if (allocated(reader%image_flags)) writer%image_flags(:, new_idx) = reader%image_flags(:, j)
                end if
            end do
            
            ! 書き込み時間を計測
            call system_clock(start_time)
            call writer%write()
            call system_clock(end_time)
            write_time = write_time + (end_time - start_time)
        end if
        
        ! 最初の3フレームだけ詳細に表示
        if (frame >= 3) then
            print *, "... continuing to read frames ..."
            ! 以降は進捗だけ表示
            if (mod(frame, 10) == 0) then
                print *, "Read frame: ", frame
            end if
        end if
    end do

    ! ファイルを閉じる
    call reader%close()
    if (test_write) call writer%close()

    ! パフォーマンス情報を表示
    print *, "============ Performance Summary ============"
    print *, "Total frames read: ", frame
    print *, "Total read time (seconds): ", real(read_time) / real(clock_rate)
    print *, "Average read time per frame (seconds): ", real(read_time) / real(frame) / real(clock_rate)
    
    if (test_write) then
        print *, "Total write time (seconds): ", real(write_time) / real(clock_rate)
        print *, "Average write time per frame (seconds): ", real(write_time) / real(frame) / real(clock_rate)
        print *, "Output file: ", trim(output_filename)
    end if

contains
    subroutine debug_print(i)
        integer, intent(inout) :: i

        print *, "Debug: ", i
        i = i + 1
    end subroutine debug_print
end program test_lammpstrj

