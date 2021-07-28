program main

!*****************************************************************************80
!
!!  MAIN is the main program for WAVE_MPI.
!
  use mpi

  integer ( kind = 4 ) id
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) p
  double precision :: wtime

  call MPI_Init ( ierr )

  call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )

  call MPI_Comm_size ( MPI_COMM_WORLD, p, ierr )

  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  WAVE_MPI:'
    write ( *, '(a)' ) '  FORTRAN90/MPI version.'
    write ( *, '(a)' ) '  Solve the 1D time-dependent wave equation.'
  end if
!
!  Record the starting time.
!
  if ( id == 0 ) then
    wtime = MPI_Wtime ( )
  end if

  call solution ( id, p )
!
!  Record the final time.
!
  if ( id == 0 ) then
    wtime = MPI_Wtime ( ) - wtime
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6)' ) '  Wall clock elapsed seconds = ', wtime
  end if
!
!  Terminate MPI.
!
  call MPI_Finalize ( ierr )
!
!  Terminate.
!
  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'WAVE_MPI:'
    write ( *, '(a)' ) '  Normal end of execution.'
  end if

  stop
end




!***********************************************************************************
subroutine solution ( id, p )

!
!! SOLUTION computes the solution of the 1D wave equation.
!

  use mpi

  integer , parameter :: n = 100
  double precision ::  Transition_wave(1:800)
  double precision ::  Exact_sol(1:800)
  double precision ::  u_ex(1:n)
  double precision ::  u(0:n+1)
  double precision ::  u_new(1:n)
  integer :: i
  integer :: id
  integer :: p
  double precision :: initial_condition
  integer :: j
  integer :: j_max
  integer :: j_min
  integer :: k
  integer :: l
  integer :: status(MPI_STATUS_SIZE)
  integer :: tag
  double precision ::  time
  double precision ::  time_delta
  double precision ::  time_max
  double precision ::  time_min
  double precision ::  time_new
  double precision ::  x(0:n+1)
  double precision ::  x_delta
  double precision ::  x_max
  double precision ::  x_min
  double precision ::  Nc
  double precision, parameter :: c=0.5006257822277846d0, kh=0.5d0,alp=8.0d0
  character*11 :: filestr
 

 
  j_max = 400
  j_min = 0
  time_max = 5.0D+00
  time_min = 0.0D+00
  x_max = 5.0D+00
  x_min = 0.0D+00
!
!  Have process 0 print out some information.
!
  if ( id == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Compute an approximate solution to the time dependent'
    write ( *, '(a)' ) '  one dimensional wave equation:'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '    du/dt - c * du/dx = 0 '
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6,a,g14.6)' ) &
    '  for ', x_min, ' = x_min < x < x_max = ', x_max
    write ( *, '(a)' ) ' '
    write ( *, '(a,g14.6,a,g14.6)' ) &
    '  and ', time_min, ' = time_min < t <= t_max = ', time_max
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  Initial condition is specified at time_min.'
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) '  The finite difference method is used to discretize'
    write ( *, '(a)' ) '  internal points of the differential equation.'
    write ( *, '(a)' ) '  and Forward diference for left boundary and Backward for right boundary'
    write ( *, '(a)' ) '  '
    write ( *, '(a,i8,a)' ) '  This uses', p * n, ' equally spaced points in X'
    write ( *, '(a,i8,a)' ) '  and', j_max, ' equally spaced points in time.'
    write ( *, '(a)' ) ' '
    write ( *, '(a,i8,a)' ) &
      '  Parallel execution is done using', p, ' processors.'
    write ( *, '(a)' ) '  Domain decomposition is used.'
    write ( *, '(a,i8,a)' ) '  Each processor works on', n, ' nodes,'
    write ( *, '(a)' ) &
      '  and shares some information with its immediate neighbors.'
  end if
!
!  Set the X coordinates of the N nodes.
!  We don't actually need ghost values of X but we'll throw them in
!  as X(0) and X(N+1).

  do i = 0, n + 1
    x(i) = ( dble ( id * n + i - 1 ) * x_max + dble ( p * n - id * n - i ) * x_min ) / dble (p * n - 1 )
  end do


  time_delta = ( time_max - time_min ) / dble ( j_max - j_min )
  x_delta = ( x_max - x_min ) / dble ( p * n - 1 )

  Nc = c * time_delta/x_delta
  print*, "CFL number is:", Nc
!  Set the values of u at the initial time.
!
  time = time_min

  u(0) = 0.0D+00
  do i = 1, n
    
    u(i) = dexp(-alp * (x(i)-1.0d0)**2) * dsin(kh * (x(i)) / x_delta)
 
  end do
  u(n+1) = 0.0D+00



!  Compute the values of u at the next time, based on current data.

  do j = 1, j_max
    
    

    time_new = ( dble(j - j_min ) * time_max + dble( j_max - j) * time_min ) / dble ( j_max - j_min )

!   Send u(1) to ID-1.

    if ( id > 0 ) then
      tag = 1
      call MPI_Send ( u(1), 1, MPI_DOUBLE_PRECISION, id-1, tag, &
        MPI_COMM_WORLD, ierr )
    end if

!   Receive u(N+1) from ID+1.

    if ( id < p - 1 ) then
      tag = 1
      call MPI_Recv ( u(n+1), 1,  MPI_DOUBLE_PRECISION, id+1, tag, &
        MPI_COMM_WORLD, status, ierr )
    end if

!   Send u(N) to ID+1.

    if ( id < p - 1 ) then
      tag = 2
      call MPI_Send ( u(n), 1, MPI_DOUBLE_PRECISION, id+1, tag, &
        MPI_COMM_WORLD, ierr )
    end if

!   Receive u(0) from ID-1.

    if ( 0 < id ) then
      tag = 2
      call MPI_Recv ( u(0), 1, MPI_DOUBLE_PRECISION, id-1, tag, &
        MPI_COMM_WORLD, status, ierr )
    end if

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

!   Update the temperature based on the four point stencil.

    do i = 1, n
 
      u_new(i) = u(i) - 0.5d0 * Nc * (u(i+1) - u(i-1)) + 0.5d0 * Nc**2 * (u(i+1) -2*u(i) + u(i-1))
      
    end do

!   u at the extreme left and right boundaries was incorrectly computed
!   using the differential equation.  Replace that calculation by
!   the boundary conditions.

    if ( 0 == id ) then
    !Forward differencing
      u_new(1) = u(1) - c * time_delta * (u(2)-u(1)) / x_delta
    end if

    if ( id == p - 1 ) then
    !Backward differencing
      u_new(n) = u(n) - c * time_delta * (u(n)-u(n-1)) / x_delta
     
    end if

    !exact solution
    do k = 1,n
        u_ex(k) = dexp(-alp * ((x(k)-c * time)-1.0d0)**2) * dsin(kh * (x(k)-c * time) / x_delta)
    end do

    !use gather library to send message to rank zero
    call MPI_Gather( u_new, n, MPI_DOUBLE_PRECISION, Transition_wave, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Gather( u_ex, n, MPI_DOUBLE_PRECISION, Exact_sol, n, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

    if (id == 0) then

    filestr(1:1) = 't'
    filestr(2:2) = ACHAR(48+int(time/100000))
    filestr(3:3) = ACHAR(48+int(time/10000)-10*int(time/100000))
    filestr(4:4) = ACHAR(48+int(time/1000)-10*int(time/10000))
    filestr(5:5) = ACHAR(48+int(time/100)-10*int(time/1000))
    filestr(6:6) = ACHAR(48+int(time/10)-10*int(time/100))
    filestr(7:7) = ACHAR(48+int(time/1)-10*int(time/100))
    filestr(8:11) = '.dat'

    open(1,file = filestr, status = 'replace')
    write(1,*) 'Variable = Numerical Solution, Exact Solution'

    do l = 1, 800
    write (1,*)  Transition_wave(l), Exact_sol(l)
    enddo

    close(1)

    end if

    !Update time and temperature.
    time = time_new
    u(1:n) = u_new(1:n)
    u(0) = 0.0D+00
    u(n+1) = 0.0D+00

  end do

  return
end

function initial_condition ( x, time, x_delta )

!*****************************************************************************80
!
!! INITIAL_CONDITION evaluates the initial conditions.
!
!  Parameters:
!
!    Input, real  X, TIME, grid spacing.
!
!    Output, real  INITIAL_CONDITION, the initial condition.
!
  implicit none

  double precision:: initial_condition
  double precision :: time
  double precision :: x
  double precision, parameter :: c=2.5d0, kh=0.5d0,alp=8.0d0 
  double precision :: x_delta
  


  
  initial_condition = dexp(-alp * ((x-c * time)-1.0d0)**2) * dsin(kh * (x-c * time) / x_delta)
    
 
  return 
end

subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.

!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h - 12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end

