!*****************************************************************************************
!>
!  PIKAIA is a general purpose unconstrained optimization
!  method based on a genetic algorithm.
!  This is an object-oriented version of the algorithm for Fortran 2003/2008.
!
!# See also
!  * [Original description page](http://www.hao.ucar.edu/modeling/pikaia/pikaia.php)
!  * [Original sourcecode](http://download.hao.ucar.edu/archive/pikaia/)
!
!# License
!
!    Copyright (c) 2015, Jacob Williams
!
!    http://github.com/jacobwilliams/pikaia
!
!    All rights reserved.
!
!    Redistribution and use in source and binary forms, with or without
!    modification, are permitted provided that the following conditions are met:
!    * Redistributions of source code must retain the above copyright notice, this
!      list of conditions and the following disclaimer.
!    * Redistributions in binary form must reproduce the above copyright notice,
!      this list of conditions and the following disclaimer in the documentation
!      and/or other materials provided with the distribution.
!    * Neither the name of pikaia nor the names of its
!      contributors may be used to endorse or promote products derived from
!      this software without specific prior written permission.
!
!    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!    AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!    IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!    FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!    DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!    CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!    OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!    ------------------------------------------------------------------------------
!
!    The original version of the PIKAIA software is public domain software
!    and is available electronically from the High Altitude Observatory.
!    http://www.hao.ucar.edu/modeling/pikaia/pikaia.php
!
!
!# History
!  * Jacob Williams : 3/8/2015 : Significant refactoring of original PIKAIA code.
!    Converted to free-form source, double precision real variables, added various
!    new features, and an object-oriented interface.
!
!*****************************************************************************************

    module pikaia_module

    use,intrinsic :: iso_fortran_env

    implicit none

    private

    integer,parameter :: wp = real64  !! Default real kind [8 bytes]. Ĭ��ʵ��Ϊ˫����
    integer,parameter :: I1B = int8   !Ⱦɫ���������8�ֽھ��㹻�� κ������ [-128,127] �иĶ�
    integer,parameter :: IB = int32   !����Ⱦɫ������Ʊ���   �иĶ� ������

    !*********************************************************
    type,public :: pikaia_class !pikaia��

        !! Main class for using the Pikaia algorithm.
        !! INIT and SOLVE are the only public methods. ����������init �� solve

        private

        integer(I1B) :: n = 0  !number of solution variables �Ա������� �иĶ�
        real(wp),dimension(:),allocatable :: xl    !! lower bounds of x ��Ա����Ϊ��̬����
        real(wp),dimension(:),allocatable :: xu    !! upper bound of x ��Ա����Ϊ��̬����
        real(wp),dimension(:),allocatable :: del !�Ա��������򳤶� ��Ա����Ϊ��̬����

        !other solution inputs (with default values):
        integer  :: np                 = 100 !��Ⱥ�� ����Ϊż��
        integer  :: ngen               = 500 !������
        integer(I1B)  :: nd            = 5   !��Ч���� �иĶ�
        integer(I1B)  :: nb            = 30  !��Ч���� �иĶ� ������
        real(wp)      :: z             = 2**30  !�иĶ� ������
        real(wp) :: pcross             = 0.85_wp !�������
        integer(I1B)  :: imut          = 2 !����ģʽ �иĶ�
        real(wp) :: pmuti              = 0.005_wp  !������ʳ�ʼֵ
        real(wp) :: pmutmn             = 0.0005_wp !��С�������
        real(wp) :: pmutmx             = 0.25_wp   !���������
        real(wp) :: fdif               = 1.0_wp !�����Ӧ�Ȳ��� �����ӳ���select_parents
        integer(I1B)  :: irep          = 1 !��ֳ�ƻ� �иĶ�
        integer(I1B)  :: ielite        = 1 !�Ƿ��Ŵ���Ӣ����,0��1�� �����ӳ���newpop  �иĶ�
        integer(I1B)  :: ivrb          = 0 !��ӡģʽ 0��1��С3�߳�  �иĶ�
        real(wp) :: convergence_tol    = 0.0001_wp !�����б� ���ε�����Ӧ�Ȳ�ֵ
        integer :: convergence_window  = 20 !��������  �иĶ�
        integer  :: iseed              = 999 !���������
        real(wp) :: initial_guess_frac = 0.1_wp !��ʼ��Ⱥ�²��� �����ӳ���pikaia ��10%����Ⱥ�Ա�����ʼֵ���ò²�ֵ

        !used internally:
        real(wp) :: pmut   = -huge(1.0_wp)
        real(wp) :: bestft = huge(1.0_wp) !��ǿ��Ӧ��
        real(wp) :: pmutpv = huge(1.0_wp) !��ǿ��Ӧ�ȱ�����

        !user-supplied procedures: ����ָ��user_f����ָ��������pikaia_func��ͬ�β��б�Ĺ���
        procedure(pikaia_func),pointer :: user_f => null()  !! fitness function
        procedure(iter_func),pointer   :: iter_f => null()  !! reporting function (best member of population)

    contains

        !public routines:
        procedure,non_overridable,public :: init   => set_inputs !���캯�� ��������
        procedure,non_overridable,public :: solve  => solve_with_pikaia !������� �������� [xl,xu]

        !private routines:
        procedure,non_overridable :: ff  => func_wrapper  !! internal pikaia function (x:[0,1]) ��Ӧ�Ⱥ��� ������Ϊ[0,1]
        procedure,non_overridable :: newpop !����Ⱥ�滻����Ⱥ irep=1
        procedure,non_overridable :: stdrep !�㶨״̬��ֳ irep=2��3
        procedure,non_overridable :: genrep !ȫ������滻 irep=1
        procedure,non_overridable :: adjmut !�����ʵĶ�̬���� imut=2,3,5,6
        procedure,non_overridable :: cross  !����
        procedure,non_overridable :: encode !����
        procedure,non_overridable :: mutate !���� imut>=4��<4
        procedure,non_overridable :: decode !����
        procedure,non_overridable :: select_parents !���̶��㷨ѡ��ĸ
        procedure,non_overridable :: report !��ӡ������Ϣ ivrb>0
        procedure,non_overridable :: rnkpop !��Ⱥ��Ӧ������
        procedure,non_overridable :: pikaia !������� [0,1]

    end type pikaia_class
    !*********************************************************

    abstract interface !����ӿ�

        subroutine pikaia_func(me,x,f)  !��Ӧ�Ⱥ����ӿ� ����ʵ�ּ�������
            !! The interface for the function that pikaia will be maximizing.
        import :: wp,pikaia_class !����ģ��������� pikaia_classΪ����
        implicit none
        class(pikaia_class),intent(inout)  :: me    !! pikaia class
        real(wp),dimension(:),intent(in)   :: x     !! optimization variable vector
        real(wp),intent(out)               :: f     !! fitness value
        end subroutine pikaia_func

        subroutine iter_func(me,iter,x,f) !���ε������ ����ʵ�ּ�������
            !! The interface for the function that user can specify
            !! to report each iteration when pikaia is running.
            !! The best (fittest) population member is passed to
            !! this routine in each generation.
        import :: wp,pikaia_class !����ģ��������� pikaia_classΪ����
        implicit none
        class(pikaia_class),intent(inout)  :: me    !! pikaia class
        integer,intent(in)                 :: iter  !! iteration number
        real(wp),dimension(:),intent(in)   :: x     !! optimization variable vector
        real(wp),intent(in)                :: f     !! fitness value
        end subroutine iter_func

    end interface

    contains
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!
!  Constructor for the [[pikaia_class]].
!  The routine must be called before the solve routine can be used.
!
!  The following inputs are required: n, f, xl, xu.
!  For the others, if they are not present, then
!  the default values are used
!
!@note Based on setctl in the original code.

    subroutine set_inputs(me,&
                            n,xl,xu,f,status,&
                            iter_f,&
                            np,ngen,nd,pcross,pmutmn,pmutmx,pmut,imut,&
                            fdif,irep,ielite,ivrb,&
                            convergence_tol,convergence_window,initial_guess_frac,&
                            iseed) !���캯�� ȷ����Ա������ֵ

    implicit none

    class(pikaia_class),intent(out)    :: me        !! pikaia class
    integer(I1B),intent(in)            :: n         !! the parameter space dimension, i.e., the number �иĶ�
                                                    !! of adjustable parameters (size of the x vector).
    real(wp),dimension(n),intent(in)   :: xl        !! vector of lower bounds for x ���뾲̬����
    real(wp),dimension(n),intent(in)   :: xu        !! vector of upper bounds for x ���뾲̬����
    procedure(pikaia_func)             :: f         !! user-supplied scalar function of n variables, �κ���pikaia_func�ӿ�һ���������̶��ɴ���
                                                    !! which must have the [[pikaia_func]] procedure interface.
                                                    !! By convention, f should return higher values for more optimal
                                                    !! parameter values (i.e., individuals which are more "fit").
                                                    !! For example, in fitting a function through data points, f
                                                    !! could return the inverse of chi**2.
    integer(I1B),intent(out)           :: status    !! status output flag (0 if there were no errors) ��������  �иĶ�
    procedure(iter_func),optional      :: iter_f    !! user-supplied subroutine that will report the �κ���iter_func�ӿ�һ���������̶��ɴ���
                                                    !! best solution for each generation.
                                                    !! It must have the [[iter_func]] procedure interface.  If not present,
                                                    !! then it is not used.  (note: this is independent of ivrb).
    integer,intent(in),optional        :: np        !! number of individuals in a population (default is 100)
    integer,intent(in),optional        :: ngen      !! maximum number of iterations
    integer(I1B),intent(in),optional   :: nd        !! number of significant digits (i.e., number of �иĶ�
                                                    !! genes) retained in chromosomal encoding (default is 6).
    real(wp),intent(in),optional       :: pcross    !! crossover probability; must be  <= 1.0 (default
                                                    !! is 0.85). If crossover takes place, either one
                                                    !! or two splicing points are used, with equal
                                                    !! probabilities
    real(wp),intent(in),optional       :: pmutmn    !! minimum mutation rate; must be >= 0.0 (default is 0.0005)
    real(wp),intent(in),optional       :: pmutmx    !! maximum mutation rate; must be <= 1.0 (default is 0.25)
    real(wp),intent(in),optional       :: pmut      !! initial mutation rate; should be small (default
                                                    !! is 0.005) (Note: the mutation rate is the probability
                                                    !! that any one gene locus will mutate in
                                                    !! any one generation.)
    integer(I1B),intent(in),optional   :: imut      !! mutation mode; 1/2/3/4/5 (default is 2). �иĶ�
                                                    !!  1=one-point mutation, fixed rate.
                                                    !!  2=one-point, adjustable rate based on fitness.
                                                    !!  3=one-point, adjustable rate based on distance.
                                                    !!  4=one-point+creep, fixed rate.
                                                    !!  5=one-point+creep, adjustable rate based on fitness.
                                                    !!  6=one-point+creep, adjustable rate based on distance.
    real(wp),intent(in),optional       :: fdif      !! relative fitness differential; range from 0 ѡ��ĸʱʹ��
                                                    !! (none) to 1 (maximum).  (default is 1.0)
    integer(I1B),intent(in),optional   :: irep      !! reproduction plan; 1/2/3=Full generational �иĶ�
                                                    !! replacement/Steady-state-replace-random/Steady-
                                                    !! state-replace-worst (default is 3)
    integer(I1B),intent(in),optional   :: ielite    !! elitism flag; 0/1=off/on (default is 0) ��Ӣ���� �иĶ�
                                                    !! (Applies only to reproduction plans 1 and 2)
    integer(I1B),intent(in),optional   :: ivrb      !! printed output 0/1/2=None/Minimal/Verbose  0��1��С3�߳� �иĶ�
                                                    !! (default is 0)
    real(wp),intent(in),optional       :: convergence_tol    !! convergence tolerance; must be > 0.0 (default is 0.0001)
    integer,intent(in),optional        :: convergence_window !! convergence window; must be >= 0  �иĶ�
                                                             !! This is the number of consecutive solutions
                                                             !! within the tolerance for convergence to
                                                             !! be declared (default is 20)
    real(wp),intent(in),optional       :: initial_guess_frac !! fraction of the initial population
                                                             !! to set equal to the initial guess.  Range from 0
                                                             !! (none) to 1.0 (all). (default is 0.1 or 10%).
    integer,intent(in),optional        :: iseed              !! random seed value; must be > 0 (default is 999)

    me%n = n

    if (allocated(me%xl)) deallocate(me%xl)
    allocate(me%xl(n)) !��̬����ȷ���ڴ�
    me%xl = xl

    if (allocated(me%xu)) deallocate(me%xu)
    allocate(me%xu(n)) !��̬����ȷ���ڴ�
    me%xu = xu

    if (allocated(me%del)) deallocate(me%del)
    allocate(me%del(n)) !��̬����ȷ���ڴ�
    me%del = me%xu - me%xl !�Ա��������򳤶�

    me%user_f => f !����ָ��ָ��������

    if (present(iter_f)) me%iter_f => iter_f !����ָ��ָ��������

    if (present(np                 )) me%np                 = np
    if (present(ngen               )) me%ngen               = ngen
    if (present(nd                 )) me%nd                 = nd
    if (present(pcross             )) me%pcross             = pcross
    if (present(imut               )) me%imut               = imut
    if (present(pmut               )) me%pmuti              = pmut  !initial value
    if (present(pmutmn             )) me%pmutmn             = pmutmn
    if (present(pmutmx             )) me%pmutmx             = pmutmx
    if (present(fdif               )) me%fdif               = fdif
    if (present(irep               )) me%irep               = irep
    if (present(ielite             )) me%ielite             = ielite
    if (present(ivrb               )) me%ivrb               = ivrb
    if (present(convergence_tol    )) me%convergence_tol    = convergence_tol
    if (present(convergence_window )) me%convergence_window = convergence_window
    if (present(initial_guess_frac )) me%initial_guess_frac = initial_guess_frac
    if (present(iseed              )) me%iseed              = iseed
    me%nb=CEILING(log(1.0_wp*10**me%nd)/log(2.0_wp)) !ʮ������Ч����ת��Ϊ��������Ч����  �иĶ� ������
    me%z=2**me%nb-1 !�иĶ� ������
    !check for errors:

    !initialize error flag:
    status = 0

    !Print a header
    if (me%ivrb>0) then !��ӡ��ͷ
        write(output_unit,'(A)') '------------------------------------------------------------'
        write(output_unit,'(A)') '              PIKAIA Genetic Algorithm Report               '
        write(output_unit,'(A)') '------------------------------------------------------------'
        write(output_unit,'(A,I4)')    ' Number of Generations evolving: ',me%ngen
        write(output_unit,'(A,I4)')    '     Individuals per generation: ',me%np
        write(output_unit,'(A,I4)')    '  Number of Chromosome segments: ',me%n
        write(output_unit,'(A,I4)')    '  Length of Chromosome segments: ',me%nd
        write(output_unit,'(A,I4)')    '    Binary Length of Chromosome: ',me%nb !�иĶ� ������
        write(output_unit,'(A,E11.4)') '          Crossover probability: ',me%pcross
        write(output_unit,'(A,E11.4)') '          Initial mutation rate: ',me%pmuti
        write(output_unit,'(A,E11.4)') '          Minimum mutation rate: ',me%pmutmn
        write(output_unit,'(A,E11.4)') '          Maximum mutation rate: ',me%pmutmx
        write(output_unit,'(A,E11.4)') '  Relative fitness differential: ',me%fdif
        write(output_unit,'(A,E11.4)') '         Initial guess fraction: ',me%initial_guess_frac
        write(output_unit,'(A,E11.4)') '          Convergence tolerance: ',me%convergence_tol
        write(output_unit,'(A,I4)')    '             Convergence window: ',me%convergence_window
        select case (me%imut) !����ģʽ
        case(1); write(output_unit,'(A)') '                  Mutation Mode: Uniform, Constant Rate'
        case(2); write(output_unit,'(A)') '                  Mutation Mode: Uniform, Variable Rate (F)'
        case(3); write(output_unit,'(A)') '                  Mutation Mode: Uniform, Variable Rate (D)'
        case(4); write(output_unit,'(A)') '                  Mutation Mode: Uniform+Creep, Constant Rate'
        case(5); write(output_unit,'(A)') '                  Mutation Mode: Uniform+Creep, Variable Rate (F)'
        case(6); write(output_unit,'(A)') '                  Mutation Mode: Uniform+Creep, Variable Rate (D)'
        end select
        select case (me%irep)  !��ֳģʽ
        case(1); write(output_unit,'(A)') '              Reproduction Plan: Full generational replacement'
        case(2); write(output_unit,'(A)') '              Reproduction Plan: Steady-state-replace-random'
        case(3); write(output_unit,'(A)') '              Reproduction Plan: Steady-state-replace-worst'
        end select
        write(output_unit,'(A)') '------------------------------------------------------------'
    end if

    !Check some control values ������ģʽ
    if (me%imut/=1 .and. me%imut/=2 .and. me%imut/=3 .and. &
          me%imut/=4 .and. me%imut/=5 .and. me%imut/=6) then
       write(output_unit,'(A)') ' ERROR: illegal value for Mutation Mode.'
       status = 5
    end if

    if (me%fdif>1.0_wp) then
       write(output_unit,'(A)') ' ERROR: illegal value for Relative fitness differential.'
       status = 9
    end if

    if (me%irep/=1 .and. me%irep/=2 .and. me%irep/=3) then !��鷱ֳ����ģʽ
       write(output_unit,'(A)') ' ERROR: illegal value for Reproduction plan.'
       status = 10
    end if

    if (me%pcross>1.0_wp .or. me%pcross<0.0_wp) then !��齻�����
       write(output_unit,'(A)') ' ERROR: illegal value for Crossover probability.'
       status = 4
    end if

    if (me%ielite/=0 .and. me%ielite/=1) then
       write(output_unit,'(A)') ' ERROR: illegal value for Elitism flag.'
       status = 11
    end if

    if (me%convergence_tol<=0.0_wp) then
        write(output_unit,'(A)') ' ERROR: illegal value for Convergence tolerance.'
        status = 101
    end if

    if (me%convergence_window<=0) then
        write(output_unit,'(A)') ' ERROR: illegal value for Convergence window.'
        status = 102
    end if

    if (me%iseed<=0) then
        write(output_unit,'(A)') ' ERROR: illegal value for iseed.'
        status = 103
    end if

    if (me%nd>9 .or. me%nd<1) then !�����Ч����
        write(output_unit,'(A)') ' ERROR: illegal value for Chromosome length.'
        status = 104
    end if

    if (mod(me%np,2)>0) then  !�����Ⱥ���Ƿ�Ϊż��
       write(output_unit,'(A)') ' ERROR: population size must be an even number.'
       status = 105
    end if

    if (me%initial_guess_frac<0.0_wp .or. me%initial_guess_frac>1.0_wp) then
       write(output_unit,'(A)') ' ERROR: illegal value for Initial guess fraction.'
       status = 106
    end if
    !�����Ǿ���
    if (me%irep==1 .and. me%imut==1 .and. me%pmuti>0.5_wp .and. me%ielite==0) then
       write(output_unit,'(A)') &
        ' WARNING: dangerously high value for Initial mutation rate; '//&
        '(Should enforce elitism with ielite=1.)'
    end if

    if (me%irep==1 .and. me%imut==2 .and. me%pmutmx>0.5_wp .and. me%ielite==0) then
       write(output_unit,'(A)') &
       ' WARNING: dangerously high value for Maximum mutation rate; '//&
       '(Should enforce elitism with ielite=1.)'
    end if

    if (me%fdif<0.33_wp .and. me%irep/=3) then
       write(output_unit,'(A)') &
       ' WARNING: dangerously low value of Relative fitness differential.'
    end if

    block !�иĶ� ������
        integer(IB),pointer::a
        allocate(a)
        if (me%nb > bit_size(a) - 1) then
            write(output_unit,'(A)') &
            ' WARNING: IB is too small,shoud be int64.' 
        end if
        deallocate(a)
    end block

    end subroutine set_inputs
!*****************************************************************************************

!*****************************************************************************************
!>
!  Optimization (maximization) of user-supplied "fitness" function
!  over n-dimensional parameter space x using a basic genetic
!  algorithm method.
!
!  Genetic algorithms are heuristic search techniques that
!  incorporate in a computational setting, the biological notion
!  of evolution by means of natural selection.  This subroutine
!  implements the three basic operations of selection, crossover,
!  and mutation, operating on "genotypes" encoded as strings.
!
!  Version 1.2 differs from version 1.0 (December 1995) in that
!  it includes (1) two-point crossover, (2) creep mutation, and
!  (3) dynamical adjustment of the mutation rate based on metric
!  distance in parameter space.
!
!# Authors
!   * Paul Charbonneau & Barry Knapp
!     (High Altitude Observatory, National Center for Atmospheric Research)
!     Version 1.2 [ 2002 April 3 ]
!   * Jacob Williams : 3/8/2015 : Refactoring and some new features.
!
!# References
!   * Charbonneau, Paul. "An introduction to genetic algorithms for
!     numerical optimization", NCAR Technical Note TN-450+IA
!     (April 2002)
!   * Charbonneau, Paul. "Release Notes for PIKAIA 1.2",
!     NCAR Technical Note TN-451+STR (April 2002)
!   * Charbonneau, Paul, and Knapp, Barry. "A User's Guide
!     to PIKAIA 1.0" NCAR Technical Note TN-418+IA
!     (December 1995)
!   * Goldberg, David E.  Genetic Algorithms in Search, Optimization,
!     & Machine Learning.  Addison-Wesley, 1989.
!   * Davis, Lawrence, ed.  Handbook of Genetic Algorithms.
!     Van Nostrand Reinhold, 1991.

    subroutine pikaia(me,x,f,status) !x[0,1] �иĶ� ������

    implicit none

    !subroutine arguments:
    class(pikaia_class),intent(inout)      :: me
    real(wp),dimension(:),intent(inout)    :: x      !! Input - initial guess for solution vector. [0,1]
                                                     !! Output - the "fittest" (optimal) solution found, [0,1]
                                                     !! i.e., the solution which maximizes the fitness function.
    real(wp),intent(out)                   :: f      !! the (scalar) value of the fitness function at x
    integer(I1B),intent(out)               :: status !! an indicator of the success or failure �иĶ�
                                                     !! of the call to pikaia (0=success; non-zero=failure)

    !Local variables
    integer  :: k,ip,ig,ip1,ip2,new,newtot,istart,i_window
    real(wp) :: current_best_f, last_best_f, fguess
    logical  :: convergence
    real(wp),dimension(me%n,2)     :: ph    !��ĸ��Ⱥ�Ա��� [0,1]
    real(wp),dimension(me%n,me%np) :: oldph !����Ⱥ  [0,1]
    real(wp),dimension(me%n,me%np) :: newph !����Ⱥ  [0,1]
    integer(IB),dimension(me%n)   :: gn1 !ĸ��ȺȾɫ�����  �иĶ� ������
    integer(IB),dimension(me%n)   :: gn2 !����ȺȾɫ�����  �иĶ� ������
    integer,dimension(me%np)       :: ifit !ifit�����±�
    integer,dimension(me%np)       :: jfit !jfit��������
    real(wp),dimension(me%np)      :: fitns !��Ⱥ��Ӧ��
    real(wp),dimension(me%n)       :: xguess !��Ⱥ�Ա����²�ֵ

    real(wp),parameter :: big = huge(1.0_wp)    !! a large number

    !initialize:
    call rninit(me%iseed)
    me%bestft   = -big !����Ϊ��Сֵ
    me%pmutpv   = -big !����Ϊ��Сֵ
    me%pmut     = me%pmuti  !set initial mutation rate (it can change)
    i_window    = 0
    last_best_f = -big !����Ϊ��Сֵ
    convergence = .false. !������
    status      = 0

    !Handle the initial guess:
    if (me%initial_guess_frac==0.0_wp) then !��ʱ�²�ֵû�� ��Ⱥȫ���������ֵ

        !initial guess not used (totally random population)

        istart = 1  !index to start random population members �����Ⱥ������

    else

        !use the initial guess:

        xguess = x !�Ա����²�ֵ[0,1]
        do k=1,me%n    !make sure they are all within the [0,1] bounds
            xguess(k) = max( 0.0_wp, min(1.0_wp,xguess(k)) ) !�²�ֵ�߽���
        end do
        call me%ff(xguess,fguess) !����²���Ⱥ��Ӧ�Ⱥ���

        !how many elements in the population to set to xguess?:
        ! [at least 1, at most n]
        istart = max(1, min(me%np, int(me%np * me%initial_guess_frac)))

        do k=1,istart !ǰistart����Ⱥ���ò²�ֵ
            oldph(:,k) = xguess
            fitns(k)   = fguess
        end do

        istart = istart + 1  !index to start random population members

    end if
    !ǰistart����Ⱥ֮��Ĳ������ֵ
    !Compute initial (random but bounded) phenotypes
    do ip=istart,me%np
        do k=1,me%n
            oldph(k,ip)=urand()  !from [0,1]
        end do
        call me%ff(oldph(:,ip),fitns(ip)) !���������Ⱥ��Ӧ�Ⱥ���
    end do

    !Rank initial population by fitness order ��ʼ��Ⱥ��Ӧ������
    call me%rnkpop(fitns,ifit,jfit)  !ifit�����±� jfit��������

    !Main Generation Loop
    do ig=1,me%ngen !����

        !Main Population Loop
        newtot=0
        do ip=1,me%np/2 !һ����Ⱥ

            !1. pick two parents ѡ��ĸ
            call me%select_parents(jfit,ip1,ip2) !jfit��������
            !ip1ĸ��Ⱥ  ip2����Ⱥ
            !2. encode parent phenotypes
            call me%encode(oldph(:,ip1),gn1) !ʵ��ת��Ϊʮ��������
            call me%encode(oldph(:,ip2),gn2) !ʵ��ת��Ϊʮ��������

            !3. breed ����
            call me%cross(gn1,gn2) !����
            call me%mutate(gn1) !ĸ��ȺȾɫ�����
            call me%mutate(gn2) !����ȺȾɫ�����

            !4. decode offspring genotypes
            call me%decode(gn1,ph(:,1)) !ʮ��������ת��Ϊʵ��
            call me%decode(gn2,ph(:,2)) !ʮ��������ת��Ϊʵ��

            !5. insert into population
            if (me%irep==1) then !��ֳģʽ Full generational replacement
                call me%genrep(ip,ph,newph)
            else  !ifit�����±� jfit��������
                call me%stdrep(ph,oldph,fitns,ifit,jfit,new) !������Ⱥoldph ��������
                newtot = newtot+new !newtot����Ⱥ���¸�����Ŀ
            end if

        end do    !End of Main Population Loop

        !if running full generational replacement: swap populations
        if (me%irep==1) call me%newpop(oldph,newph,ifit,jfit,fitns,newtot) !��ֳģʽ Full generational replacement
        !newtot����Ⱥ���¸�����Ŀ ��������Ⱥ��Ӧ������ ��ʱoldph�Ѿ��Ǹ��º���Ⱥ oldph=newph
        !adjust mutation rate? �����ʶ�̬����
        if (any(me%imut==[2,3,5,6])) call adjmut(me,oldph,fitns,ifit)

        !report this iteration:
        if (me%ivrb>0) call me%report(oldph,fitns,ifit,ig,newtot) !��ϸ��ӡ������Ϣ

        !report (unscaled) x:
        if (associated(me%iter_f)) & !����ָ��me%iter_f��ָ��
            call me%iter_f(ig,me%xl+me%del*oldph(1:me%n,ifit(me%np)),fitns(ifit(me%np))) !��ӡÿ�ε�����ǿ��Ⱥ������Ӧ��

        !JW additions: add a convergence criteria
        ! [stop if the last convergence_window iterations are all within the convergence_tol]
        current_best_f = fitns(ifit(me%np))    !current iteration best fitness
        if (abs(current_best_f-last_best_f)<=me%convergence_tol) then
            !this solution is within the tol from the previous one
            i_window = i_window + 1    !number of solutions within the convergence tolerance
        else
            i_window = 0    !a significantly better solution was found, reset window
        end if
        !����convergence_window�ε����� ��ǰ���ε�����Ӧ��������С��convergence_tol����Ϊ������
        if (i_window>=me%convergence_window) then
            convergence = .true.
            exit !exit main loop -> convergence
        end if
        last_best_f = current_best_f    !to compare with next iteration ��һ�ε�����ǿ��Ӧ��

    end do    !End of Main Generation Loop

    !JW additions:
    if (me%ivrb>0) then !��ϸ��ӡ������Ϣ
        if (convergence) then
            write(output_unit,'(A)') 'Solution Converged'
        else
            write(output_unit,'(A)') 'Iteration Limit Reached'
        end if
    end if

    !Return best phenotype and its fitness
    x = oldph(1:me%n,ifit(me%np)) !��ǿ��Ӧ�ȶ�Ӧ�Ա���[0,1]
    f = fitns(ifit(me%np)) !��ǿ��Ӧ��

    end subroutine pikaia
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!
!  Main pikaia wrapper used by the class.

    subroutine solve_with_pikaia(me,x,f,status) ![xl,xu]

    implicit none

    class(pikaia_class),intent(inout)   :: me
    real(wp),dimension(:),intent(inout) :: x !���봫����Ϊ[xl,xu]
    real(wp),intent(out)                :: f
    integer(I1B),intent(out)            :: status !�иĶ�

    if (associated(me%user_f)) then !����ȷ��Ӧ�Ⱥ���

        !scale input initial guess to be [0,1]:
        x = (x-me%xl)/me%del !���Ա�����Χת��Ϊ[0,1]

        !call the main routine, using the wrapper function:
        call me%pikaia(x,f,status) ![0,1]

        !unscale output to be [xl,xu]:
        x = me%xl + me%del*x !���Ա�����Χ��[0,1]�ָ���ԭʼ

    else !��Ӧ�Ⱥ���û�й��������庯��

        write(output_unit,'(A)') 'Error: pikaia class not initialized.'
        status = -1

    end if

    end subroutine solve_with_pikaia
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!
!  Wrapper for the user's function that is used by the main pikaia routine
!  The x input to this function comes from pikaia, and will be between [0,1].

    subroutine func_wrapper(me,x,f) !������Ϊ[0,1]��Ӧ��Ӧ�Ⱥ���

    implicit none

    class(pikaia_class),intent(inout) :: me   ! pikaia class
    real(wp),dimension(:),intent(in)  :: x    ! optimization variable vector [0,1]
    real(wp),intent(out)              :: f    ! fitness value

    real(wp),dimension(me%n) :: xp    !unscaled x vector: [xu,xl]

    !map each x variable from [0,1] to [xl,xu]:
    xp = me%xl + me%del*x !���Ա�����Χ��[0,1]�ָ���ԭʼ

    !call the user's function with xp:
    call me%user_f(xp,f)

    end subroutine func_wrapper
!*****************************************************************************************

!*****************************************************************************************
!> author: B. G. Knapp
!  date: 86/12/23
!
!  Return integer array p which indexes array a in increasing order.
!  Array a is not disturbed.  The Quicksort algorithm is used.
!
!# Reference
!   * N. Wirth, "Algorithms and Data Structures", Prentice-Hall, 1986

    subroutine rqsort(n,a,p) !��������p��ʵ������a�������е��±�

    implicit none

    integer,intent(in)               :: n
    real(wp),dimension(n),intent(in) :: a
    integer,dimension(n),intent(out) :: p

    !Constants
    integer,parameter :: LGN = 32      !! log base 2 of maximum n
    integer,parameter :: Q   = 11      !! smallest subfile to use quicksort on

    !Local:
    integer,dimension(LGN) :: stackl,stackr
    real(wp) :: x
    integer :: s,t,l,m,r,i,j

    !Initialize the stack
    stackl(1)=1
    stackr(1)=n
    s=1

    !Initialize the pointer array
    p = [(i, i=1,n)]

    do while (s>0)

        l=stackl(s)
        r=stackr(s)
        s=s-1

3       if ((r-l)<Q) then

            !Use straight insertion
            do i=l+1,r
                t = p(i)
                x = a(t)
                do j=i-1,l,-1
                    if (a(p(j))<=x) goto 5
                    p(j+1) = p(j)
                end do
                j=l-1
5               p(j+1) = t
            end do

        else

            !Use quicksort, with pivot as median of a(l), a(m), a(r)
            m=(l+r)/2
            t=p(m)
            if (a(t)<a(p(l))) then
                p(m)=p(l)
                p(l)=t
                t=p(m)
            end if
            if (a(t)>a(p(r))) then
                p(m)=p(r)
                p(r)=t
                t=p(m)
                if (a(t)<a(p(l))) then
                    p(m)=p(l)
                    p(l)=t
                    t=p(m)
                end if
            end if

            !Partition
            x=a(t)
            i=l+1
            j=r-1
            do while (i<=j)

                do while (a(p(i))<x)
                    i=i+1
                end do

                do while (x<a(p(j)))
                    j=j-1
                end do

                if (i<=j) then
                    t=p(i)
                    p(i)=p(j)
                    p(j)=t
                    i=i+1
                    j=j-1
                end if

            end do

            !Stack the larger subfile
            s=s+1
            if ((j-l)>(r-i)) then
                stackl(s)=l
                stackr(s)=j
                l=i
            else
                stackl(s)=i
                stackr(s)=r
                r=j
            end if

            goto 3
        end if

    end do

    end subroutine rqsort
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 3/8/2015
!
!  Return the next pseudo-random deviate from a sequence which is
!  uniformly distributed in the interval [0,1]
!
!@note This is now just a wrapper for the intrinsic random_number function.

    function urand() result(r) !�������������װ

    implicit none

    real(wp) :: r

    call random_number(r)

    end function urand
!*****************************************************************************************

!*****************************************************************************************
!> author: Jacob Williams
!  date: 3/8/2015
!
!  Initialize the random number generator with the input seed value.
!
!@note This is now just a wrapper for the intrinsic random_seed function.

    subroutine rninit(iseed) !������������װ

    implicit none

    integer,intent(in) :: iseed

    integer,dimension(:),allocatable :: seed
    integer :: n

    call random_seed(size=n)
    allocate(seed(n))
    seed = iseed
    call random_seed(put=seed)
    deallocate(seed)

    end subroutine rninit
!*****************************************************************************************

!*****************************************************************************************
!>
!  Write generation report to standard output

    subroutine report(me,oldph,fitns,ifit,ig,nnew) !��ӡ������Ϣ

    implicit none

    class(pikaia_class),intent(inout)         :: me
    real(wp),dimension(me%n,me%np),intent(in) :: oldph
    real(wp),dimension(me%np),intent(in)      :: fitns
    integer,dimension(me%np),intent(in)       :: ifit
    integer,intent(in)                        :: ig   !�������
    integer,intent(in)                        :: nnew !newtot����Ⱥ���¸�����Ŀ

    integer :: ndpwr,k
    logical :: rpt

    rpt=.false.

    if (me%pmut/=me%pmutpv) then !����ʱ�����ʷ����仯
       me%pmutpv=me%pmut
       rpt=.true.
    end if

    if (fitns(ifit(me%np))/=me%bestft) then !����ʱ��ǿ��Ӧ�ȷ����仯
       me%bestft=fitns(ifit(me%np))
       rpt=.true.
    end if

    if (rpt .or. me%ivrb>=2) then

        !Power of 10 to make integer genotypes for display
        ndpwr = 10**me%nd

        write(output_unit,'(/I6,I6,F10.6,4F10.6)') & !��ӡ��ǿ��Ӧ�� ��ǿ��Ӧ�� ��λ����Ӧ��
            ig,nnew,me%pmut,fitns(ifit(me%np)),&
            fitns(ifit(me%np-1)),fitns(ifit(me%np/2))

        do k=1,me%n
            write(output_unit,'(22X,3I10)') &         !��ӡ��ӦȾɫ�����
                    nint(ndpwr*oldph(k,ifit(me%np))),&
                    nint(ndpwr*oldph(k,ifit(me%np-1))),&
                    nint(ndpwr*oldph(k,ifit(me%np/2)))
        end do

    end if

    end subroutine report
!*****************************************************************************************

!*****************************************************************************************
!>
!  Encode phenotype parameters into integer genotype
!  ph(k) are x,y coordinates [ 0 < x,y < 1 ]

    subroutine encode(me,ph,gn) !ʵ��ת��Ϊ����������  �иĶ� ������

    implicit none

    class(pikaia_class),intent(inout)         :: me
    real(wp),dimension(me%n),intent(in)       :: ph
    integer(IB),dimension(me%n),intent(out)   :: gn  !�иĶ� ������

    gn=0_IB !��ֵΪ0
    gn=ph*me%z

    end subroutine encode
!*****************************************************************************************

!*****************************************************************************************
!>
!  decode genotype into phenotype parameters
!  ph(k) are x,y coordinates [ 0 < x,y < 1 ]

    subroutine decode(me,gn,ph) !����������ת��Ϊʵ��  �иĶ� ������

    implicit none

    class(pikaia_class),intent(inout)        :: me
    integer(IB),dimension(me%n),intent(in)   :: gn !�иĶ� ������
    real(wp),dimension(me%n),intent(out)     :: ph

    ph=gn/me%z !�иĶ� ������

    end subroutine decode
!*****************************************************************************************

!*****************************************************************************************
!>
!  breeds two parent chromosomes into two offspring chromosomes.
!  breeding occurs through crossover. If the crossover probability
!  test yields true (crossover taking place), either one-point or
!  two-point crossover is used, with equal probabilities.
!
!@note Compatibility with version 1.0: To enforce 100% use of one-point
!      crossover, un-comment appropriate line in source code below

    subroutine cross(me,gn1,gn2) !���� �иĶ� ������

    implicit none

    class(pikaia_class),intent(inout)           :: me
    integer(IB),dimension(me%n),intent(inout)   :: gn1  !�иĶ� ������
    integer(IB),dimension(me%n),intent(inout)   :: gn2  !�иĶ� ������

    integer :: i, ispl, ispl2, itmp
    integer(IB) :: t

    !Use crossover probability to decide whether a crossover occurs
    if (urand()<me%pcross) then

        !Compute first crossover point
        ispl=int(urand()*me%nb)+1 ![1,me%nb] ����λ�����

        !Now choose between one-point and two-point crossover
        if (urand()<0.5_wp) then !����λ���յ�
            ispl2=me%nb
        else
            ispl2=int(urand()*me%nb)+1 ![1,me%nb]
            !Un-comment following line to enforce one-point crossover
            if (ispl2<ispl) then !������ȷ��ispl1<isp2
                itmp=ispl2
                ispl2=ispl
                ispl=itmp
            end if
        end if

        !Swap genes from ispl to ispl2
        do i=1,me%n !��ĸ����Ⱦɫ��
            !CALL MVBITS(FROM, FROMPOS, LEN, TO, TOPOS)
            !������FROM��FROMPOSλ��FROMPOS+LEN-1λ��Ϣ���Ƶ�����TO��TOPOSλ��TOPOS+LEN-1 LENΪ����λ�ĳ���
            t=0_IB
            CALL MVBITS(gn2(i), ispl-1, ispl2-ispl+1, t, ispl-1)
            CALL MVBITS(gn1(i), ispl-1, ispl2-ispl+1, gn2(i), ispl-1)
            CALL MVBITS(t, ispl-1, ispl2-ispl+1, gn1(i), ispl-1)
        end do

    end if

    end subroutine cross
!*****************************************************************************************

!*****************************************************************************************
!>
!  Introduces random mutation in a genotype.
!  Mutations occur at rate pmut at all gene loci.
!
!# Input
!   * imut=1    Uniform mutation, constant rate
!   * imut=2    Uniform mutation, variable rate based on fitness
!   * imut=3    Uniform mutation, variable rate based on distance
!   * imut=4    Uniform or creep mutation, constant rate
!   * imut=5    Uniform or creep mutation, variable rate based on fitness
!   * imut=6    Uniform or creep mutation, variable rate based on distance

    subroutine mutate(me,gn) !���� �иĶ� ������

    implicit none

    class(pikaia_class),intent(inout)           :: me
    integer(IB),dimension(me%n),intent(inout)   :: gn  !�иĶ� ������

    integer :: i,j,inc

    !Decide which type of mutation is to occur
    if (me%imut>=4 .and. urand()<=0.5_wp) then !�������+���

        !CREEP MUTATION OPERATOR ��伴Ⱦɫ����������1���1 �����С
        !Subject each locus to random +/- 1 increment at the rate pmut
        do i=1,me%n
            do j=1,me%nb
                if (urand()<me%pmut) then
                    inc=nint( urand() )*2-1 !-1��1
                    gn(i)=gn(i)+shiftl(1_IB,j-1)*inc
                end if
                !���Ʒ�Χ
                if (gn(i)<0_IB) gn(i)=0_IB
                if (gn(i)>me%z) gn(i)=me%z
            end do
        end do

    else !�������

        !UNIFORM MUTATION OPERATOR ͳһ������� ����ϴ� ���ܳ���0��7��7��0
        !Subject each locus to random mutation at the rate pmut
        do i=1,me%n !Ⱦɫ��ÿ��λ��
            do j=1,me%nb,3
                if (urand()<me%pmut) then
                    !int(urand()*8) ![0,7] or [000,111]
                    CALL MVBITS(int(urand()*8), 0, 3, gn(i), j-1)
                    !���Ʒ�Χ
                    if (gn(i)<0_IB) gn(i)=0_IB
                    if (gn(i)>me%z) gn(i)=me%z
                end if
            end do
        end do

    end if

    end subroutine mutate
!*****************************************************************************************

!*****************************************************************************************
!>
!  Dynamical adjustment of mutation rate:
!
!   * imut=2 or imut=5 : adjustment based on fitness differential
!     between best and median individuals
!   * imut=3 or imut=6 : adjustment based on metric distance
!     between best and median individuals

    subroutine adjmut(me,oldph,fitns,ifit) !�����ʵĶ�̬����

    implicit none

    class(pikaia_class),intent(inout)         :: me
    integer,dimension(me%np),intent(in)       :: ifit
    real(wp),dimension(me%n,me%np),intent(in) :: oldph
    real(wp),dimension(me%np),intent(in)      :: fitns

    integer  :: i
    real(wp) :: rdif

    real(wp),parameter :: rdiflo = 0.05_wp
    real(wp),parameter :: rdifhi = 0.25_wp
    real(wp),parameter :: delta  = 1.5_wp

    if (me%imut==2 .or. me%imut==5) then !����Ӧ�ȵ���������

        !Adjustment based on fitness differential
        rdif = abs(fitns(ifit(me%np)) - &
               fitns(ifit(me%np/2)))/(fitns(ifit(me%np)) + &
               fitns(ifit(me%np/2)))

    else if (me%imut==3 .or. me%imut==6) then !���������������

        !Adjustment based on normalized metric distance
        rdif=0.0_wp
        do i=1,me%n
            rdif=rdif+( oldph(i,ifit(me%np))-oldph(i,ifit(me%np/2)) )**2
        end do
        rdif=sqrt( rdif ) / real(me%n,wp)

    end if

    if (rdif<=rdiflo) then
        me%pmut=min(me%pmutmx,me%pmut*delta)
    else if (rdif>=rdifhi) then
        me%pmut=max(me%pmutmn,me%pmut/delta)
    end if

    end subroutine adjmut
!*****************************************************************************************

!*****************************************************************************************
!>
!  Selects two parents from the population, using roulette wheel
!  algorithm with the relative fitnesses of the phenotypes as
!  the "hit" probabilities.
!
!# Reference
!  * Davis 1991, chap. 1.
!
!# History
!  * Jacob Williams : 3/10/2015 : rewrote this routine to return both parents,
!    and also protect against the loop exiting without selecting a parent.
    !���̶��㷨ѡ��ĸ
    subroutine select_parents(me,jfit,imom,idad) !jfit��������

    implicit none

    class(pikaia_class),intent(inout)   :: me
    integer,dimension(me%np),intent(in) :: jfit
    integer,intent(out)                 :: imom
    integer,intent(out)                 :: idad

    integer  :: np1,i,j
    real(wp) :: dice,rtfit
    integer,dimension(2) :: parents

    !initialize:
    np1 = me%np+1
    parents = -99

    !get two (unequal) parents:
    do j=1,2
        main: do
            dice = urand()*me%np*np1 ![0,np**2+np)
            rtfit = 0.0_wp
            do i=1,me%np
                rtfit = rtfit+np1+me%fdif*(np1-2*jfit(i))
                if (rtfit>=dice) then
                    parents(j) = i
                    if (parents(1)/=parents(2)) exit main
                end if
            end do
        end do main
    end do
    !��Ӧ��Խ�� jfitԽС(������) rtfit��Խ��
    imom = parents(1) !ĸ��Ⱥλ��
    idad = parents(2) !����Ⱥλ��
    !me%fdif=1ʱ ��1�� 2*np1-2 ��2�� 2*np1-4 ���� ���1�� 2
    end subroutine select_parents
!*****************************************************************************************

!*****************************************************************************************
!>
!  Ranks initial population.
!  Calls external sort routine to produce key index and rank order
!  of input array arrin (which is not altered).

    subroutine rnkpop(me,arrin,indx,rank) !indx�����±� rank��������
    !��arrin(indx(1))��С ����rankΪme%np arrin(indx(me%np))��� ����rankΪ1
    implicit none

    class(pikaia_class),intent(inout)     :: me
    real(wp),dimension(me%np),intent(in)  :: arrin
    integer,dimension(me%np),intent(out)  :: indx
    integer,dimension(me%np),intent(out)  :: rank

    integer :: i

    !Compute the key index
    call rqsort(me%np,arrin,indx)  !��������indx��ʵ������arrin�������е��±�

    !and the rank order
    do i=1,me%np
        rank(indx(i)) = me%np-i+1 !��������
    end do

    end subroutine rnkpop
!*****************************************************************************************

!*****************************************************************************************
!>
!  Full generational replacement: accumulate offspring into new
!  population array

    subroutine genrep(me,ip,ph,newph) !ȫ������滻 irep=1

    implicit none

    class(pikaia_class),intent(inout)          :: me
    integer,intent(in)                         :: ip !��Ⱥλ��
    real(wp),dimension(me%n,2),intent(in)      :: ph !��ĸ��Ⱥ�Ա��� [0,1]
    real(wp),dimension(me%n,me%np),intent(out) :: newph !����Ⱥ [0,1]

    integer :: i1,i2,k

    !Insert one offspring pair into new population
    i1=2*ip-1   !ip 1 2 3 ... np/2
    i2=i1+1     !i1 1 3 5 ... np-1
    do k=1,me%n !i2 2 4 6 ... np
        newph(k,i1)=ph(k,1)
        newph(k,i2)=ph(k,2)
    end do

    end subroutine genrep
!*****************************************************************************************

!*****************************************************************************************
!>
!  Steady-state reproduction: insert offspring pair into population
!  only if they are fit enough (replace-random if irep=2 or
!  replace-worst if irep=3).

    subroutine stdrep(me,ph,oldph,fitns,ifit,jfit,nnew) !�㶨״̬��ֳ irep=2��3

    implicit none

    class(pikaia_class),intent(inout)             :: me
    real(wp),dimension(me%n,2),intent(in)         :: ph
    real(wp),dimension(me%n,me%np),intent(inout)  :: oldph
    real(wp),dimension(me%np),intent(inout)       :: fitns
    integer,dimension(me%np),intent(inout)        :: ifit !�����±�
    integer,dimension(me%np),intent(inout)        :: jfit !��������
    integer,intent(out)                           :: nnew !����Ⱥ���¸�����Ŀ

    integer  :: i,j,k,i1,if1
    real(wp) :: fit

    nnew = 0

    main_loop : do j=1,2

        !1. compute offspring fitness (with caller's fitness function)
        call me%ff(ph(:,j),fit) !���㸸ĸ��Ⱥ��Ӧ��

        !2. if fit enough, insert in population
        do i=me%np,1,-1 !����

            if (fit>fitns(ifit(i))) then !��ĸ��Ⱥ����

                !make sure the phenotype is not already in the population
                if (i<me%np) then
                    if (all(oldph(1:me%n,ifit(i+1))==ph(1:me%n,j))) cycle main_loop
                end if

                !offspring is fit enough for insertion, and is unique

                !(i) insert phenotype at appropriate place in population
                if (me%irep==3) then
                    i1=1 !ֻ���������Ⱥ
                else if (me%ielite==0 .or. i==me%np) then !me%irep==2
                    i1=int(urand()*me%np)+1 ![1,me%np]
                else  !me%irep==2
                    i1=int(urand()*(me%np-1))+1 ![1,me%np-1]
                end if
                if1 = ifit(i1) !�����±� ����Ⱥλ��
                fitns(if1)=fit
                do k=1,me%n
                    oldph(k,if1)=ph(k,j) !������Ⱥ
                end do

                !(ii) shift and update ranking arrays
                if (i<i1) then

                    !shift up ����
                    jfit(if1)=me%np-i !��������
                    do k=i1-1,i+1,-1
                        jfit(ifit(k))=jfit(ifit(k))-1
                        ifit(k+1)=ifit(k) !�����±�
                    end do
                    ifit(i+1)=if1

                else !i>=i1

                    !shift down ����
                    jfit(if1)=me%np-i+1 !��������
                    do k=i1+1,i
                        jfit(ifit(k))=jfit(ifit(k))+1
                        ifit(k-1)=ifit(k) !�����±�
                    end do
                    ifit(i)=if1

                end if

                nnew = nnew+1 !����Ⱥ��
                cycle main_loop

            end if

        end do

    end do main_loop

    end subroutine stdrep
!*****************************************************************************************

!*****************************************************************************************
!>
!  Replaces old population by new; recomputes fitnesses & ranks
!
!# History
!  * Jacob Williams : 3/9/2015 : avoid unnecessary function evaluation if `ielite/=1`.

    subroutine newpop(me,oldph,newph,ifit,jfit,fitns,nnew) !irep==1

    implicit none

    class(pikaia_class),intent(inout)            :: me
    real(wp),dimension(me%n,me%np),intent(inout) :: oldph
    real(wp),dimension(me%n,me%np),intent(inout) :: newph
    integer,dimension(me%np),intent(out)         :: ifit !�����±�
    integer,dimension(me%np),intent(out)         :: jfit !��������
    real(wp),dimension(me%np),intent(out)        :: fitns
    integer,intent(out)                          :: nnew !����Ⱥ���¸�����Ŀ

    integer  :: i
    real(wp) :: f

    nnew = me%np

    if (me%ielite==1) then !��Ӣ��Ⱥ����

        !if using elitism, introduce in new population fittest of old
        !population (if greater than fitness of the individual it is
        !to replace)
        call me%ff(newph(:,1),f) !��������Ⱥ��1��������Ӧ��

        if (f<fitns(ifit(me%np))) then !����Ӧ�ȱ�������Ӧ�Ȳ�
            newph(:,1)=oldph(:,ifit(me%np)) !����Ⱥ��1�������滻�ɾ�Ӣ
            nnew = nnew-1
        end if

    end if

    !replace population
    do i=1,me%np

        oldph(:,i)=newph(:,i) !������Ⱥ

        !get fitness using caller's fitness function
        call me%ff(oldph(:,i),fitns(i)) !������Ӧ��

    end do

    !compute new population fitness rank order
    call me%rnkpop(fitns,ifit,jfit)  !������Ⱥ����

    end subroutine newpop
!*****************************************************************************************

!*****************************************************************************************
    end module pikaia_module
!*****************************************************************************************
