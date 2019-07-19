!*****************************************************************************************
!> author: Jacob Williams
!  date: 3/9/2015
!
!  Sample driver program for pikaia.
!  Based on the [xpikaia](http://www.hao.ucar.edu/modeling/pikaia/xpikaia.f) routine.

    program pikaia_test

    use pikaia_module, only: pikaia_class
    use pikaia_binary_module, only: pikaia_binary_class
    use,intrinsic :: iso_fortran_env, wp => real64, I1B => int8 !双精度实数 有改动

    implicit none
    !n为自变量个数
    integer(I1B),parameter  :: n = 2  !! dimension of problem (number of optimization variables) 有改动

    integer                 :: seed !有改动
    integer(I1B)            :: status !有改动
    real(wp),dimension(n)   :: x !自变量数组
    real(wp)                :: f !函数值 适应度
    integer                 :: iunit,istat
    real(wp),dimension(n)   :: xl,xu !自变量下界、上界数组
    !type(pikaia_class) :: p !父类对象
    type(pikaia_binary_class) :: p !子类对象
    logical                 :: header_written
    real                    :: tstart,tend !起始时间 结束时间
    !输出文件名
    character(len=*),parameter :: filename = 'pikaia_test.txt'

    seed = 999  ! specify the random seed

    !output file:
    open(newunit=iunit,file=filename,iostat=istat)
    if (istat/=0) stop 'error opening output file.'
    header_written = .false. !表头打印关闭

    !initial guess:
    write(output_unit,'(A)') ''
    write(output_unit,'(A)') ' TWOD Example'
    write(output_unit,'(A)') ''

    x  = 0.0_wp !猜测值[xl,xu]
    xl = 0.0_wp !自变量下界
    xu = 1.0_wp !自变量上界

    !initialize the class:
    call p%init(n,xl,xu,twod,status,&
                iter_f              = report_iteration,&
                ngen                = 1000,&
                nd                  = 9_I1B,&
                ivrb                = 1_I1B,&    !0,1,2
                convergence_tol     = 1.0e-6_wp,&
                convergence_window  = 200,&
                iseed               = seed)
    !call p%SetIsGray(IsGray=.true.) !F:Binary T:Gray
    !imut  1   2  3    4   5   6
    !ngen  51 75 203  28  22  24 superclass
    !ngen 103 85 154  38 143 213 subclass method1
    !ngen  51 32 305 206  98 343 subclass method2
    !Now call pikaia:
    call cpu_time(tstart)
    call p%solve(x,f,status) ![xl,xu]
    call cpu_time(tend)

    !Print the results:
    write(output_unit,'(A)') ''
    write(output_unit,'(A,1X,*(I4))')    '  status: ',status
    write(output_unit,'(A,1X,*(F12.6))') '       x: ',x !最强适应度对应自变量
    write(output_unit,'(A,1X,*(F12.6))') '       f: ',f !最强适应度
    write(output_unit,'(A)') ''
    write(output_unit,'(A,1X,F12.6,A)')  'run time: ',tend-tstart,' sec'
    write(output_unit,'(A)') ''

    !-------------------------

    !initial guess:
    write(output_unit,'(A)') ''
    write(output_unit,'(A)') ' ROSENBROCK Example'
    write(output_unit,'(A)') ''

    x  = 0.5_wp !猜测值[xl,xu]
    xl = 0.0_wp !自变量下界
    xu = 2.0_wp !自变量上界

    !initialize the class:
    call p%init(n,xl,xu,rosenbrock,status,&
                np                  = 500,&        !try a larger population for this one
                ngen                = 1000,&
                nd                  = 9_I1B,&
                convergence_tol     = 1.0e-10_wp,& !tighter tolerance also
                convergence_window  = 200,&
                iseed               = seed)
    !call p%SetIsGray(IsGray=.true.) !F:Binary T:Gray
    !Now call pikaia:
    call cpu_time(tstart)
    call p%solve(x,f,status) ![xl,xu]
    call cpu_time(tend)

    !Print the results:
    write(output_unit,'(A)') ''
    write(output_unit,'(A,1X,*(I4))')    '  status: ',status
    write(output_unit,'(A,1X,*(F12.6))') '       x: ',x !最强适应度对应自变量
    write(output_unit,'(A,1X,*(F12.6))') '       f: ',f !最强适应度
    write(output_unit,'(A)') ''
    write(output_unit,'(A,1X,F12.6,A)')  'run time: ',tend-tstart,' sec'
    write(output_unit,'(A)') ''

    close(iunit,iostat=istat)

    contains !抽象接口的具体实现

        subroutine twod(me,x,f)

        !! Compute sample fitness function
        !! (a smooth 2-d landscape)

        implicit none

        !Input/output:
        class(pikaia_class),intent(inout) :: me
        real(wp),dimension(:),intent(in)  :: x
        real(wp),intent(out)              :: f

        !Constant
        real(wp),parameter  :: pi=acos(-1.0_wp)
        real(wp),parameter  :: sigma2=0.15_wp
        integer,parameter   :: nn = 9

        !Local
        real(wp) :: rr

        if (x(1)>1.0_wp .or. x(2)>1.0_wp) then
            write(output_unit,*) 'Error in function twod: invalid inputs.'
            write(output_unit,*) 'x(1)=',x(1),'>1.0'
            write(output_unit,*) 'x(2)=',x(2),'>1.0'
            stop
        else
            rr=sqrt( (0.5_wp-x(1))**2+ (0.5_wp-x(2))**2)
            f=cos(rr*nn*pi)**2 *exp(-rr**2/sigma2) !默认是求最大值
        end if

        end subroutine twod

        subroutine rosenbrock(me,x,f)

        !! Rosenbrock function for testing the algorithm.
        !! The minimum is at f(1,1) = 0.
        !!
        !!#Reference
        !!   * [Rosenbrock function](http://en.wikipedia.org/wiki/Rosenbrock_function)

        implicit none

        class(pikaia_class),intent(inout)   :: me
        real(wp),dimension(:), intent(in)   :: x
        real(wp),intent(out)                :: f

        real(wp),parameter :: one     = 1.0_wp
        real(wp),parameter :: hundred = 100.0_wp

        !the rosenbrock function:
        f = (one-x(1))**2 + hundred*(x(2)-x(1)**2)**2

        f = -f    !since pikaia maximizes !求最小值

        end subroutine rosenbrock

        subroutine report_iteration(me,iter,x,f) !历次迭代结果

        !! A simple iteration reporting function.
        !! Writes `iter,x,f` to the output file.

        implicit none

        class(pikaia_class),intent(inout)  :: me
        integer,intent(in)                 :: iter
        real(wp),dimension(:),intent(in)   :: x
        real(wp),intent(in)                :: f

        character(len=10),dimension(n) :: xheader !表头
        integer :: i

        !the first time it is called, also write a header:
        if (.not. header_written) then
            do i=1,n
                write(xheader(i),'(I10)') i !整数转换为字符串
                xheader(i) = 'X'//trim(adjustl(xheader(i)))
                xheader(i) = repeat(' ',10-len_trim(xheader(i)))//xheader(i)
            end do
            write(iunit,'(A5,1X,*(A10,1X))') 'ITER',xheader,'F'
            header_written = .true.
        end if

        write(iunit,'(I5,1X,*(F10.6,1X))') iter,x,f

        end subroutine report_iteration

    end program pikaia_test
!*****************************************************************************************
