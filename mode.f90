subroutine init_output_dir
    use com_date
    implicit none

    integer :: values(8), ios, unit_id
    logical :: exists
    character(len=256) :: latest_file
    character(len=512) :: cmd

    output_root = 'output'
    output_dir = ''
    restart_input_dir = ''
    latest_file = trim(output_root) // '/latest_run.txt'

    cmd = 'mkdir -p ' // trim(output_root)
    call system(trim(cmd))

    inquire(file=trim(latest_file), exist=exists)
    if (exists) then
        unit_id = 902
        open(unit_id, file=trim(latest_file), status='old', action='read', iostat=ios)
        if (ios .eq. 0) then
            read(unit_id, '(A)', iostat=ios) restart_input_dir
            if (ios .ne. 0) restart_input_dir = ''
            close(unit_id)
        endif
    endif

    call date_and_time(values=values)
    write(output_dir,'(a,"/run_",i4.4,i2.2,i2.2,"_",i2.2,i2.2,i2.2,"_",i3.3)') &
        trim(output_root), values(1), values(2), values(3), values(5), values(6), values(7), values(8)

    cmd = 'mkdir -p ' // trim(output_dir)
    call system(trim(cmd))

    unit_id = 903
    open(unit_id, file=trim(latest_file), status='replace', action='write', iostat=ios)
    if (ios .eq. 0) then
        write(unit_id, '(A)') trim(output_dir)
        close(unit_id)
    endif
end subroutine

subroutine build_output_path(filename, path)
    use com_date
    implicit none

    character(len=*), intent(in) :: filename
    character(len=*), intent(out) :: path

    path = trim(output_dir) // '/' // trim(filename)
end subroutine

subroutine build_restart_input_path(filename, path)
    use com_date
    implicit none

    character(len=*), intent(in) :: filename
    character(len=*), intent(out) :: path

    if (len_trim(restart_input_dir) .gt. 0) then
        path = trim(restart_input_dir) // '/' // trim(filename)
    else
        path = trim(filename)
    endif
end subroutine

subroutine load_motion_table
    use com_date
    implicit none

    integer :: unit_id, ios, count, i
    real(kind=ps) :: row(max_motion_columns)
    character(len=512) :: line

    if (motion_loaded) return

    motion_file = 'motion_params.dat'
    unit_id = 901
    count = 0

    open(unit_id, file=trim(motion_file), status='old', action='read', iostat=ios)
    if (ios .ne. 0) then
        write(*,*) '无法打开运动参数文件: ', trim(motion_file)
        stop 1
    endif

    do
        read(unit_id, '(A)', iostat=ios) line
        if (ios .ne. 0) exit
        if (len_trim(line) .eq. 0) cycle
        row = 0.0d0
        read(line, *, iostat=ios) row(1:motion_column_count)
        if (ios .ne. 0) then
            write(*,*) '运动参数文件格式错误，列数应为: ', motion_column_count, trim(motion_file)
            stop 1
        endif
        count = count + 1
    enddo

    if (count .lt. 1) then
        write(*,*) '运动参数文件为空: ', trim(motion_file)
        stop 1
    endif

    rewind(unit_id)
    allocate(motion_table(count, motion_column_count))

    i = 0
    do
        read(unit_id, '(A)', iostat=ios) line
        if (ios .ne. 0) exit
        if (len_trim(line) .eq. 0) cycle
        i = i + 1
        row = 0.0d0
        read(line, *, iostat=ios) row(1:motion_column_count)
        motion_table(i, :) = row(1:motion_column_count)
        if (ios .ne. 0) then
            write(*,*) '运动参数文件格式错误，列数应为: ', motion_column_count, trim(motion_file)
            stop 1
        endif
    enddo
    close(unit_id)

    motion_points = count
    motion_index = 1
    motion_end_time = motion_table(motion_points, 1)
    motion_loaded = .true.

    do i = 2, motion_points
        if (motion_table(i, 1) .lt. motion_table(i-1, 1)) then
            write(*,*) '运动参数文件时间列必须非递减: ', trim(motion_file)
            stop 1
        endif
    enddo

    if (motion_points .gt. 1) then
        do i = 2, motion_points
            if (motion_table(i, 1) .gt. motion_table(i-1, 1)) exit
        enddo
    endif
end subroutine

subroutine apply_motion_row(row)
    use com_date
    implicit none

    real(kind=ps), intent(in) :: row(motion_column_count)
    integer :: wing, col

    angle_arr = 0.0d0
    omega_arr = 0.0d0
    omgac_arr = 0.0d0
    sa_arr = 0.0d0
    sau_arr = 0.0d0
    saac_arr = 0.0d0
    sb_arr = 0.0d0
    sbv_arr = 0.0d0
    sbac_arr = 0.0d0

    do wing = 1, wing_count
        col = 2 + (wing - 1) * motion_columns_per_wing
        angle_arr(wing) = row(col)
        omega_arr(wing) = row(col + 1)
        omgac_arr(wing) = row(col + 2)
        sa_arr(wing) = row(col + 3)
        sau_arr(wing) = row(col + 4)
        saac_arr(wing) = row(col + 5)
        sb_arr(wing) = row(col + 6)
        sbv_arr(wing) = row(col + 7)
        sbac_arr(wing) = row(col + 8)
    enddo

end subroutine

subroutine interpolate_motion_row(i1, i2, target_time, row)
    use com_date
    implicit none

    integer, intent(in) :: i1, i2
    real(kind=ps), intent(in) :: target_time
    real(kind=ps), intent(out) :: row(motion_column_count)
    real(kind=ps) :: t1, t2, weight

    t1 = motion_table(i1, 1)
    t2 = motion_table(i2, 1)

    if (abs(t2 - t1) .le. 1.0d-14) then
        row = motion_table(i1, :)
        row(1) = target_time
        return
    endif

    weight = (target_time - t1) / (t2 - t1)
    row = (1.0d0 - weight) * motion_table(i1, :) + weight * motion_table(i2, :)
    row(1) = target_time
end subroutine

subroutine lookup_motion_row(target_time, row)
    use com_date
    implicit none

    real(kind=ps), intent(in) :: target_time
    real(kind=ps), intent(out) :: row(motion_column_count)

    if (target_time .le. motion_table(1, 1)) then
        motion_index = 1
        row = motion_table(1, :)
        row(1) = target_time
        return
    endif

    if (target_time .ge. motion_table(motion_points, 1)) then
        motion_index = max(1, motion_points - 1)
        row = motion_table(motion_points, :)
        row(1) = target_time
        return
    endif

    motion_index = max(1, min(motion_index, motion_points - 1))

    do while ((motion_index .gt. 1) .and. (target_time .lt. motion_table(motion_index, 1)))
        motion_index = motion_index - 1
    enddo

    do while ((motion_index .lt. motion_points - 1) .and. (target_time .gt. motion_table(motion_index + 1, 1)))
        motion_index = motion_index + 1
    enddo

    if (abs(target_time - motion_table(motion_index, 1)) .le. 1.0d-12) then
        row = motion_table(motion_index, :)
        row(1) = target_time
        return
    endif

    if (abs(target_time - motion_table(motion_index + 1, 1)) .le. 1.0d-12) then
        row = motion_table(motion_index + 1, :)
        row(1) = target_time
        return
    endif

    call interpolate_motion_row(motion_index, motion_index + 1, target_time, row)
end subroutine

subroutine build_motion_step_table(dt)
    use com_date
    implicit none

    real(kind=ps), intent(in) :: dt
    integer :: i
    real(kind=ps) :: target_time
    real(kind=ps) :: row(max_motion_columns)

    if (.not. motion_loaded) call load_motion_table()
    if (dt .le. 0.0d0) stop '仿真时间步必须大于0'

    motion_dt = dt
    motion_step_count = int((motion_end_time + 1.0d-12) / motion_dt) + 1

    if (allocated(motion_step_table)) deallocate(motion_step_table)
    allocate(motion_step_table(motion_step_count, motion_column_count))

    motion_index = 1
    do i = 1, motion_step_count
        target_time = dble(i - 1) * motion_dt
        call lookup_motion_row(target_time, row)
        motion_step_table(i, :) = row(1:motion_column_count)
    enddo
end subroutine

subroutine mode
    use com_date
    implicit none

    if (.not. motion_loaded) call load_motion_table()
    if (.not. allocated(motion_step_table)) stop '运动参数步表尚未初始化'
    if (motion_dt .le. 0.0d0) stop '运动参数时间步无效'

    motion_step_count = size(motion_step_table, 1)
    motion_index = nint(time / motion_dt) + 1
    motion_index = max(1, min(motion_index, motion_step_count))
    call apply_motion_row(motion_step_table(motion_index, :))
end subroutine

subroutine update_body_extra_from_rel
    use com_date
    implicit none

    integer :: wing
    real(kind=ps) :: accum_v, accum_s, accum_pos

    body_extra_v = 0.0d0
    body_extra_s = 0.0d0
    body_extra_pos = 0.0d0
    accum_v = 0.0d0
    accum_s = 0.0d0
    accum_pos = 0.0d0

    do wing = 1, wing_count
        accum_v = accum_v + body_rel_v(wing)
        accum_s = accum_s + body_rel_s(wing)
        accum_pos = accum_pos + body_rel_pos(wing)
        body_extra_v(wing) = accum_v
        body_extra_s(wing) = accum_s
        body_extra_pos(wing) = accum_pos
    enddo
end subroutine
