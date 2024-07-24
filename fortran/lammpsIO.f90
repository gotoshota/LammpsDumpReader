module lammpsIO
    implicit none

    type :: AtomIndex
        INTEGER(KIND=4) :: id   = 0
        INTEGER(KIND=4) :: mol  = 0
        INTEGER(KIND=4) :: type = 0
        INTEGER(KIND=4) :: xu   = 0
        INTEGER(KIND=4) :: yu   = 0
        INTEGER(KIND=4) :: zu   = 0
        INTEGER(KIND=4) :: x    = 0
        INTEGER(KIND=4) :: y    = 0
        INTEGER(KIND=4) :: z    = 0
        INTEGER(KIND=4) :: ix   = 0
        INTEGER(KIND=4) :: iy   = 0
        INTEGER(KIND=4) :: iz   = 0
        integer(KIND=4) :: xs   = 0
        integer(KIND=4) :: ys   = 0
        integer(KIND=4) :: zs   = 0
    end type AtomIndex

    type :: lammpstrjReader
        character(len=256) :: filename
        integer :: unit
        logical :: end_of_file
        logical :: has_atom_idx = .false.
        integer :: timestep
        integer :: nparticles
        double precision :: box_bounds(3, 3) = 0.0d0
        REAL, allocatable :: coords(:, :)
        integer, allocatable :: image_flags(:, :)
        integer, allocatable :: id(:)
        integer, allocatable :: mol(:)
        integer, allocatable :: type(:)
        type(AtomIndex) :: atom_idx
        contains
            procedure :: open => open_trajectory
            procedure :: close => close_trajectory
            procedure :: read => read_next_frame
    end type lammpstrjReader

contains

    subroutine open_trajectory(this, filename)
        class(lammpstrjReader), intent(inout) :: this
        character(len=*), intent(in) :: filename
        integer :: ios

        this%filename = filename
        open(newunit=this%unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error opening file: ", filename
            stop
        end if
        this%end_of_file = .false.
    end subroutine open_trajectory

    subroutine close_trajectory(this)
        class(lammpstrjReader), intent(inout) :: this
        close(this%unit)
    end subroutine close_trajectory

    subroutine read_next_frame(this)
        class(lammpstrjReader), intent(inout) :: this
        character(len=256) :: line
        CHARACTER(LEN=:), allocatable :: atom_header_parts(:), atom_parts(:), box_parts(:)
        integer :: ios, i, j
        integer :: nColums

        if (this%end_of_file) return

        read(this%unit, '(A)', iostat=ios) line ! ITEM: TIMESTEP
        if (ios /= 0) then
            this%end_of_file = .true.
            return
        end if
        read(this%unit, '(I10)', iostat=ios) this%timestep
        read(this%unit, '(A)', iostat=ios) line ! ITEM: NUMBER OF ATOMS
        read(this%unit, '(I10)', iostat=ios) this%nparticles
        read(this%unit, '(A)', iostat=ios) line ! ITEM: BOX BOUNDS
        read(this%unit, '(A)', iostat=ios) line ! xlo xhi (xy)
        box_parts = split_line(line)
        read(box_parts(:), *) this%box_bounds(1:size(box_parts), 1)
        read(this%unit, '(A)', iostat=ios) line ! ylo yhi (yz)
        box_parts = split_line(line)
        read(box_parts(:), *) this%box_bounds(1:size(box_parts), 2)
        read(this%unit, '(A)', iostat=ios) line ! zlo zhi (xz)
        box_parts = split_line(line)
        read(box_parts(:), *) this%box_bounds(1:size(box_parts), 3)
        if (this%has_atom_idx) then
            read(this%unit, '(A)', iostat=ios) line ! Skip Atom header line and use the same atom header
        else if (this%has_atom_idx .eqv. .false.) then
            read(this%unit, '(A)', iostat=ios) line ! Read Atom header line and aplit it 
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
                    case default
                        print *, "Error: Unknown column to read header"
                        print *, "Column: ", i
                end select
            enddo
            this%has_atom_idx = .true.
            allocate(this%coords(3, this%nparticles)) ! x, y, z 
            if (this%atom_idx%id > 0) allocate(this%id(this%nparticles))
            if (this%atom_idx%mol > 0) allocate(this%mol(this%nparticles))
            if (this%atom_idx%type > 0) allocate(this%type(this%nparticles))
            if (this%atom_idx%ix > 0) allocate(this%image_flags(3, this%nparticles))
        end if 

        do i = 1, this%nparticles
            read(this%unit, "(A)") line
            atom_parts = split_line(adjustl(line))
            do j = 1, size(atom_parts)
                if (j == this%atom_idx%x) then
                    read(atom_parts(j), *) this%coords(1, i)
                else if (j == this%atom_idx%y) then
                    read(atom_parts(j), *) this%coords(2, i)
                else if (j == this%atom_idx%z) then
                    read(atom_parts(j), *) this%coords(3, i)
                else if (j == this%atom_idx%xs) then
                    read(atom_parts(j), *) this%coords(1, i)
                else if (j == this%atom_idx%ys) then
                    read(atom_parts(j), *) this%coords(2, i)
                else if (j == this%atom_idx%zs) then
                    read(atom_parts(j), *) this%coords(3, i)
                else if (j == this%atom_idx%xu) then
                    read(atom_parts(j), *) this%coords(1, i)
                else if (j == this%atom_idx%yu) then
                    read(atom_parts(j), *) this%coords(2, i)
                else if (j == this%atom_idx%zu) then
                    read(atom_parts(j), *) this%coords(3, i)
                else if (j == this%atom_idx%id) then
                    read(atom_parts(j), *) this%id(i)
                else if (j == this%atom_idx%mol) then
                    read(atom_parts(j), *) this%mol(i)
                else if (j == this%atom_idx%type) then
                    read(atom_parts(j), *) this%type(i)
                else if (j == this%atom_idx%ix) then
                    read(atom_parts(j), *) this%image_flags(1, i)
                else if (j == this%atom_idx%iy) then
                    read(atom_parts(j), *) this%image_flags(2, i)
                else if (j == this%atom_idx%iz) then
                    read(atom_parts(j), *) this%image_flags(3, i)
                else 
                    print *, "Error: Unknown column to substitute variables"
                    print *, "Column: ", j
                    print *, this%atom_idx
                end if
            end do
        end do
    end subroutine read_next_frame

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
        if (inWord) nWords = nWords + 1  ! 最後の単語を数える

        allocate(character(len=strLen) :: substrings(nWords))
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

end module lammpsIO

