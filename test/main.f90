program test_lammpstrjReader
    use lammpsIO
    implicit none

    type(lammpstrjReader) :: reader
    character(len=256) :: filename
    integer :: i, frame
    integer :: start_time, end_time, total_time, clock_rate

    ! テスト用のファイル名を指定します。適宜変更してください。
    filename = "../data/dump.lammpstrj"

    ! ファイルを開く
    call reader%open(filename)

    frame = 0
    total_time = 0
    call system_clock(count_rate=clock_rate)

    ! EOFまで読み込む
    do
        call system_clock(start_time)
        call reader%read()
        call system_clock(end_time)

        if (reader%end_of_file) exit

        frame = frame + 1
        total_time = total_time + (end_time - start_time)

        print *, "Frame: ", frame
        print *, "Timestep: ", reader%timestep
        print *, "Number of particles: ", reader%nparticles
        print *, "Box bounds: "
        do i = 1, 3
            print *, reader%box_bounds(i, 1), reader%box_bounds(i, 2)
        end do

        print *, "Coordinates: "
        print *, "Particle ", 10, ": ", reader%coords(:, 10)
    end do

    ! ファイルを閉じる
    call reader%close()

    print *, "Total frames read: ", frame
    print *, "Total time (seconds): ", real(total_time) / real(clock_rate)
    print *, "Average time per frame (seconds): ", real(total_time) / real(frame) / real(clock_rate)

end program test_lammpstrjReader

