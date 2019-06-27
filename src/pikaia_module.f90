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

    public !Ϊ�˱�������̳и������й��е�ʵ�������ͷ���

    integer,parameter,private :: wp  = real64 !! Default real kind [8 bytes]. Ĭ��ʵ��Ϊ˫����
    integer,parameter,private :: I1B = int8   !Ⱦɫ��ʮ���Ʊ�������8�ֽھ��㹻�� [-127,128] �иĶ�
    integer,parameter,private :: IB  = int32  !����Ⱦɫ������Ʊ���  ��������
    
    !*********************************************************
    type,public :: pikaia_class !pikaia_class��
        !���๫�е�ʵ�������뷽�����ܱ�����̳�
        !! Main class for using the Pikaia algorithm.
        !! INIT and SOLVE are the only public methods. ����������init �� solve

        integer(I1B) :: n = 0_I1B  !number of solution variables �Ա������� �иĶ�
        real(wp),dimension(:),allocatable :: xl    !! lower bounds of x ��Ա����Ϊ��̬����
        real(wp),dimension(:),allocatable :: xu    !! upper bound of x ��Ա����Ϊ��̬����
        real(wp),dimension(:),allocatable :: del   !�Ա��������򳤶� ��Ա����Ϊ��̬����

        !other solution inputs (with default values):
        integer  :: np                 = 100 !��Ⱥ�� ����Ϊż��
        integer  :: ngen               = 500 !������
        integer(I1B)  :: nd            = 5_I1B   !��Ч���� �иĶ� min:4 max:9
        real(wp) :: pcross             = 0.85_wp !�������
        integer(I1B)  :: imut          = 2_I1B !����ģʽ �иĶ�
        real(wp) :: pmuti              = 0.005_wp  !������ʳ�ʼֵ
        real(wp) :: pmutmn             = 0.0005_wp !��С�������
        real(wp) :: pmutmx             = 0.25_wp   !���������
        real(wp) :: fdif               = 1.0_wp !�����Ӧ�Ȳ��� �����ӳ���select_parents
        integer(I1B)  :: irep          = 1_I1B !��ֳ�ƻ� �иĶ�
        integer(I1B)  :: ielite        = 1_I1B !�Ƿ��Ŵ���Ӣ����,0��1�� �����ӳ���newpop  �иĶ�
        integer(I1B)  :: ivrb          = 0_I1B !��ӡģʽ 0��1��С3�߳�  �иĶ�
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
        procedure,public :: init   => set_inputs !���캯�� ��������
        procedure,non_overridable,public :: solve  => solve_with_pikaia !������� �������� [xl,xu]

        !private routines:
        procedure,non_overridable :: ff  => func_wrapper  !! internal pikaia function (x:[0,1]) ��Ӧ�Ⱥ��� ������Ϊ[0,1]
        procedure,non_overridable :: newpop !����Ⱥ�滻����Ⱥ irep=1
        procedure,non_overridable :: stdrep !�㶨״̬��ֳ irep=2��3
        procedure,non_overridable :: genrep !ȫ������滻 irep=1
        procedure,non_overridable :: adjmut !�����ʵĶ�̬���� imut=2,3,5,6
        procedure,non_overridable :: select_parents !���̶��㷨ѡ��ĸ
        procedure,non_overridable :: report !��ӡ������Ϣ ivrb>0
        procedure,non_overridable :: rnkpop !��Ⱥ��Ӧ������
        !�󶨹��̵�ͨ������ ���ݽӿڲ�ͬ�������÷��� �������decimal �������binary
        generic :: cross  => cross_decimal ,cross_binary  !���� ʮ���� ������
        generic :: encode => encode_decimal,encode_binary !���� ʮ���� ������
        generic :: mutate => mutate_decimal,mutate_binary !���� ʮ���� ������ imut>=4��<4
        generic :: decode => decode_decimal,decode_binary !���� ʮ���� ������
        !�����踲�����ع��� ��Ϊgn�Ķ���ı�
        procedure :: pikaia !������� [0,1]
        procedure :: cross_decimal ,cross_binary
        procedure :: encode_decimal,encode_binary
        procedure :: mutate_decimal,mutate_binary 
        procedure :: decode_decimal,decode_binary 

    end type pikaia_class
    !*********************************************************
    
    include "interface.f90" !abstract interface ���������๲��ӿڴ���
    !ʵ�÷�����:rninit(����ʵ������pikaia����).urand,rqsort(����ʵ������rnkpop����)
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
    integer(I1B),intent(in),optional   :: nd        !! number of significant digits (i.e., number of �иĶ� min:4 max:9
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

    !check for errors:

    !initialize error flag:
    status = 0

    !Print a header
    if (me%ivrb>0_I1B) then !��ӡ��ͷ
        write(output_unit,'(A)') '------------------------------------------------------------'
        write(output_unit,'(A)') '              PIKAIA Genetic Algorithm Report               '
        write(output_unit,'(A)') '------------------------------------------------------------'
        write(output_unit,'(A,I4)')    ' Number of Generations evolving: ',me%ngen
        write(output_unit,'(A,I4)')    '     Individuals per generation: ',me%np
        write(output_unit,'(A,I4)')    '  Number of Chromosome segments: ',me%n
        write(output_unit,'(A,I4)')    '  Length of Chromosome segments: ',me%nd
        write(output_unit,'(A,E11.4)') '          Crossover probability: ',me%pcross
        write(output_unit,'(A,E11.4)') '          Initial mutation rate: ',me%pmuti
        write(output_unit,'(A,E11.4)') '          Minimum mutation rate: ',me%pmutmn
        write(output_unit,'(A,E11.4)') '          Maximum mutation rate: ',me%pmutmx
        write(output_unit,'(A,E11.4)') '  Relative fitness differential: ',me%fdif
        write(output_unit,'(A,E11.4)') '         Initial guess fraction: ',me%initial_guess_frac
        write(output_unit,'(A,E11.4)') '          Convergence tolerance: ',me%convergence_tol
        write(output_unit,'(A,I4)')    '             Convergence window: ',me%convergence_window
        select case (me%imut) !����ģʽ
        case(1_I1B); write(output_unit,'(A)') '                  Mutation Mode: Uniform, Constant Rate'
        case(2_I1B); write(output_unit,'(A)') '                  Mutation Mode: Uniform, Variable Rate (F)'
        case(3_I1B); write(output_unit,'(A)') '                  Mutation Mode: Uniform, Variable Rate (D)'
        case(4_I1B); write(output_unit,'(A)') '                  Mutation Mode: Uniform+Creep, Constant Rate'
        case(5_I1B); write(output_unit,'(A)') '                  Mutation Mode: Uniform+Creep, Variable Rate (F)'
        case(6_I1B); write(output_unit,'(A)') '                  Mutation Mode: Uniform+Creep, Variable Rate (D)'
        end select
        select case (me%irep)  !��ֳģʽ
        case(1_I1B); write(output_unit,'(A)') '              Reproduction Plan: Full generational replacement'
        case(2_I1B); write(output_unit,'(A)') '              Reproduction Plan: Steady-state-replace-random'
        case(3_I1B); write(output_unit,'(A)') '              Reproduction Plan: Steady-state-replace-worst'
        end select
        write(output_unit,'(A)') '------------------------------------------------------------'
    end if

    !Check some control values ������ģʽ
    if (me%imut/=1_I1B .and. me%imut/=2_I1B .and. me%imut/=3_I1B .and. &
          me%imut/=4_I1B .and. me%imut/=5_I1B .and. me%imut/=6_I1B) then
       write(output_unit,'(A)') ' ERROR: illegal value for Mutation Mode.'
       status = 5
    end if

    if (me%fdif>1.0_wp) then
       write(output_unit,'(A)') ' ERROR: illegal value for Relative fitness differential.'
       status = 9
    end if

    if (me%irep/=1_I1B .and. me%irep/=2_I1B .and. me%irep/=3_I1B) then !��鷱ֳ����ģʽ
       write(output_unit,'(A)') ' ERROR: illegal value for Reproduction plan.'
       status = 10
    end if

    if (me%pcross>1.0_wp .or. me%pcross<0.0_wp) then !��齻�����
       write(output_unit,'(A)') ' ERROR: illegal value for Crossover probability.'
       status = 4
    end if

    if (me%ielite/=0_I1B .and. me%ielite/=1_I1B) then
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

    if (me%nd>9_I1B .or. me%nd<4_I1B) then !�����Ч���� �иĶ�
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
    if (me%irep==1_I1B .and. me%imut==1_I1B .and. me%pmuti>0.5_wp .and. me%ielite==0_I1B) then
       write(output_unit,'(A)') &
        ' WARNING: dangerously high value for Initial mutation rate; '//&
        '(Should enforce elitism with ielite=1.)'
    end if

    if (me%irep==1_I1B .and. me%imut==2_I1B .and. me%pmutmx>0.5_wp .and. me%ielite==0_I1B) then
       write(output_unit,'(A)') &
       ' WARNING: dangerously high value for Maximum mutation rate; '//&
       '(Should enforce elitism with ielite=1.)'
    end if

    if (me%fdif<0.33_wp .and. me%irep/=3_I1B) then
       write(output_unit,'(A)') &
       ' WARNING: dangerously low value of Relative fitness differential.'
    end if

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

    subroutine pikaia(me,x,f,status) !x[0,1] �������д�÷���

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
    !�������д��� ���������಻ͬ  ����ʮ���Ʊ���
    integer(I1B),dimension(me%n*me%nd)  :: gn1 !ĸ��ȺȾɫ����� �иĶ� ʮ����
    integer(I1B),dimension(me%n*me%nd)  :: gn2 !����ȺȾɫ����� �иĶ� ʮ����
    
    include "pikaia_sub.f90" !�˶δ��븸�������๲��
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

    subroutine encode_decimal(me,ph,gn) !ʵ��[0,1]ת��Ϊʮ�������� �������

    implicit none

    class(pikaia_class),intent(inout)         :: me
    real(wp),dimension(me%n),intent(in)       :: ph      ![0,1]
    integer(I1B),dimension(me%n*me%nd),intent(out) :: gn !�иĶ�

    integer  :: ip,i,j,ii
    real(wp) :: z

    z=10**me%nd !ndΪ��Ч���� �иĶ�  z=10.0_wp**me%nd
    ii=0
    do i=1,me%n
        ip=int(ph(i)*z)
        do j=me%nd,1,-1 !��λ����λ˳��
            gn(ii+j)=mod(ip,10)
            ip=ip/10
        end do
        ii=ii+me%nd
    end do

    end subroutine encode_decimal
    
    subroutine encode_binary(me,ph,gn) !ʵ��[0,1]ת��Ϊ���������� ������д�÷���

    implicit none

    class(pikaia_class),intent(inout)         :: me
    real(wp),dimension(me%n),intent(in)       :: ph ![0,1]
    integer(IB),dimension(me%n),intent(out)   :: gn !�иĶ�

    gn=0_IB !�иĶ�

    end subroutine encode_binary
!*****************************************************************************************

!*****************************************************************************************
!>
!  decode genotype into phenotype parameters
!  ph(k) are x,y coordinates [ 0 < x,y < 1 ]

    subroutine decode_decimal(me,gn,ph) !ʮ��������ת��Ϊʵ�� �������

    implicit none

    class(pikaia_class),intent(inout)             :: me
    integer(I1B),dimension(me%n*me%nd),intent(in) :: gn !�иĶ�
    real(wp),dimension(me%n),intent(out)          :: ph ![0,1]

    integer  :: ip,i,j,ii
    real(wp) :: z

    !z=10.0_wp**(-me%nd)
    z=10**me%nd !�иĶ�
    ii=0
    do i=1,me%n
        ip=0
        do j=1,me%nd !��λ����λ˳��
            ip=10*ip+gn(ii+j)
        end do
        !ph(i)=ip*z
        ph(i)=ip/z !�иĶ�
        ii=ii+me%nd
    end do

    end subroutine decode_decimal
    
    subroutine decode_binary(me,gn,ph) !����������ת��Ϊʵ�� ������д�÷���

    implicit none

    class(pikaia_class),intent(inout)        :: me
    integer(IB),dimension(me%n),intent(in)   :: gn !�иĶ�
    real(wp),dimension(me%n),intent(out)     :: ph ![0,1]

    ph=0.0_wp !�иĶ�

    end subroutine decode_binary
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

    subroutine cross_decimal(me,gn1,gn2) !���� �������

    implicit none

    class(pikaia_class),intent(inout)           :: me
    integer(I1B),dimension(me%n*me%nd),intent(inout) :: gn1 !�иĶ�
    integer(I1B),dimension(me%n*me%nd),intent(inout) :: gn2 !�иĶ�

    integer :: i, ispl, ispl2, itmp
    integer(I1B) :: t !�иĶ�

    !Use crossover probability to decide whether a crossover occurs
    if (urand()<me%pcross) then

        !Compute first crossover point
        ispl=int(urand()*me%n*me%nd)+1 ![1,me%n*me%nd] ����λ�����

        !Now choose between one-point and two-point crossover
        if (urand()<0.5_wp) then !����λ���յ�
            ispl2=me%n*me%nd !�൱�ڵ��㽻��
        else
            ispl2=int(urand()*me%n*me%nd)+1 !�൱�����㽻��
            !Un-comment following line to enforce one-point crossover
            !ispl2=me%n*me%nd !ǿ�Ƶ��㽻��
            if (ispl2<ispl) then !������ȷ��ispl1<isp2
                itmp=ispl2
                ispl2=ispl
                ispl=itmp
            end if
        end if

        !Swap genes from ispl to ispl2
        do i=ispl,ispl2 !��ĸ����Ⱦɫ��
            t=gn2(i)
            gn2(i)=gn1(i)
            gn1(i)=t
        end do

    end if

    end subroutine cross_decimal
    
    subroutine cross_binary(me,gn1,gn2) !���� ������д�÷���

    implicit none

    class(pikaia_class),intent(inout)           :: me
    integer(IB),dimension(me%n),intent(inout)   :: gn1 !�иĶ�
    integer(IB),dimension(me%n),intent(inout)   :: gn2 !�иĶ�

    gn1=0_IB !�иĶ�
    gn2=0_IB !�иĶ�
    end subroutine cross_binary
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

    subroutine mutate_decimal(me,gn) !���� �������

    implicit none

    class(pikaia_class),intent(inout)                :: me
    integer(I1B),dimension(me%n*me%nd),intent(inout) :: gn !�иĶ�

    integer :: i,j,k,l,ist,inc,loc
    logical :: fix

    !Decide which type of mutation is to occur
    if (me%imut>=4_I1B .and. urand()<=0.5_wp) then !�������+���

        !CREEP MUTATION OPERATOR ��伴Ⱦɫ����������1���1 �����С
        !Subject each locus to random +/- 1 increment at the rate pmut
        do i=1,me%n
            do j=1,me%nd

                if (urand()<me%pmut) then

                    !Construct integer
                    loc=(i-1)*me%nd+j !Ⱦɫ��λ��
                    inc=nint( urand() )*2-1 !-1��1
                    ist=(i-1)*me%nd+1 !���Ա�����1��Ⱦɫ��λ��
                    gn(loc)=gn(loc)+inc

                    !This is where we carry over the one (up to two digits)
                    !first take care of decrement below 0 case
                    if (inc<0 .and. gn(loc)<0) then !0��1����
                        if (j==1) then !��1��Ⱦɫ��λ�ò����� ���λ
                            gn(loc)=0
                        else !���ǵ�1��Ⱦɫ��λ�ñ���
                            fix = .true.
                            do k=loc,ist+1,-1 !�൱�ڸ�λ�ü�1
                                gn(k)=9
                                gn(k-1)=gn(k-1)-1 !��λ��1����10
                                if ( gn(k-1)>=0 ) then
                                    fix = .false.
                                    exit
                                end if
                            end do
                            if (fix) then
                                !we popped under 0.00000 lower bound; fix it up
                                if ( gn(ist)<0) then !����000000
                                    do l=ist,loc
                                        gn(l)=0
                                    end do
                                end if
                            end if
                        end if
                    end if

                    if (inc>0 .and. gn(loc)>9) then !9��1���λ
                        if (j==1) then !��1��Ⱦɫ��λ�ò����� ���λ
                            gn(loc)=9
                        else !���ǵ�1��Ⱦɫ��λ�ñ���
                            fix = .true.
                            do k=loc,ist+1,-1 !�൱�ڸ�λ�ü�1
                                gn(k)=0
                                gn(k-1)=gn(k-1)+1!��λ��һ
                                if ( gn(k-1)<=9 ) then
                                    fix = .false.
                                    exit
                                end if
                            end do
                            if (fix) then
                                !we popped over 9.99999 upper bound; fix it up
                                if ( gn(ist)>9 ) then !����999999
                                    do l=ist,loc
                                        gn(l)=9
                                    end do
                                end if
                            end if
                        end if
                    end if

                end if

            end do
        end do

    else !�������

        !UNIFORM MUTATION OPERATOR ͳһ������� ����ϴ� ���ܳ���0��9��9��0
        !Subject each locus to random mutation at the rate pmut
        do i=1,me%n*me%nd !Ⱦɫ��ÿ��λ��
            if (urand()<me%pmut) then
                gn(i)=int(urand()*10.0_wp)
            end if
        end do

    end if

    end subroutine mutate_decimal
    
    subroutine mutate_binary(me,gn) !���� ������д�÷���

    implicit none

    class(pikaia_class),intent(inout)           :: me
    integer(IB),dimension(me%n),intent(inout)   :: gn !�иĶ�

    gn=0_IB
    end subroutine mutate_binary
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
