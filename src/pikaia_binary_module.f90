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

    module pikaia_binary_module

    use pikaia_module, only: pikaia_class,urand,rninit !rqsort所在实例方法rnkpop直接继承父类
    use,intrinsic :: iso_fortran_env

    implicit none

    private

    integer,parameter,private :: wp  = real64 !! Default real kind [8 bytes]. 默认实数为双精度
    integer,parameter,private :: I1B = int8   !染色体十进制编码数组8字节就足够了 [-127,128] 有改?
    integer,parameter,private :: IB  = int32  !用于染色体二进制编码  用于子类
    !IB  int16   int32  int64
    !nd    4       9     18
    !nb   14      30     60
    !*********************************************************
    type, extends(pikaia_class),public :: pikaia_binary_class !pikaia_class的子类
        !子类继承父类所有公有的实例变量和方法
        !! Main class for using the Pikaia algorithm.
        !! INIT and SOLVE are the only public methods. 公开方法：init 和 solve

        private
        !nd 11  10  9  8  7  6  5  4  3
        !nb 37  34  30 27 24 20 17 14 10
        integer(I1B)  :: nb     = 30_I1B                           !有效数字 有改动 二进制
        real(wp)      :: z      = real(ibset(0_IB,30_IB)-1_IB,wp)  !有改动 二进制 2**30-1  fast
        logical       :: IsGray = .false.                          !F:Binary T:Gray

    contains

        !public routines:
        procedure,non_overridable,public :: init   => set_inputs !构造函数 公开方法
        procedure,non_overridable,public :: SetIsGray
        !private routines:
        procedure,non_overridable :: cross_binary  !交叉
        procedure,non_overridable :: encode_binary !编码
        procedure,non_overridable :: mutate_binary !变异 imut>=4或<4
        procedure,non_overridable :: decode_binary !解码
        procedure,non_overridable :: pikaia        !计算程序 [0,1]

    end type pikaia_binary_class
    !*********************************************************
    include "interface.f90" !abstract interface 父类与子类共享接口代码
    !子类新的实用方法有:binary2gray,gray2binary
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
                            iseed) !构造函数 确定实例变量初值
    !若子类构造函数与父类同名覆盖 则接口(形参列表)必须相同需另写方法进行新的实例变量初始化
    !若子类构造函数与父类不同名 则接口不相同 新的实例变量初始化可在子类构造函数进行
    implicit none

    class(pikaia_binary_class),intent(out) :: me    !! pikaia binary class
    integer(I1B),intent(in)            :: n         !! the parameter space dimension, i.e., the number 有改动
                                                    !! of adjustable parameters (size of the x vector).
    real(wp),dimension(n),intent(in)   :: xl        !! vector of lower bounds for x 传入静态数组
    real(wp),dimension(n),intent(in)   :: xu        !! vector of upper bounds for x 传入静态数组
    procedure(pikaia_func)             :: f         !! user-supplied scalar function of n variables, 任何与pikaia_func接口一样的子例程都可传入
                                                    !! which must have the [[pikaia_func]] procedure interface.
                                                    !! By convention, f should return higher values for more optimal
                                                    !! parameter values (i.e., individuals which are more "fit").
                                                    !! For example, in fitting a function through data points, f
                                                    !! could return the inverse of chi**2.
    integer(I1B),intent(out)           :: status    !! status output flag (0 if there were no errors) 错误类型  有改动
    procedure(iter_func),optional      :: iter_f    !! user-supplied subroutine that will report the 任何与iter_func接口一样的子例程都可传入
                                                    !! best solution for each generation.
                                                    !! It must have the [[iter_func]] procedure interface.  If not present,
                                                    !! then it is not used.  (note: this is independent of ivrb).
    integer,intent(in),optional        :: np        !! number of individuals in a population (default is 100)
    integer,intent(in),optional        :: ngen      !! maximum number of iterations
    integer(I1B),intent(in),optional   :: nd        !! number of significant digits (i.e., number of 有改动
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
    integer(I1B),intent(in),optional   :: imut      !! mutation mode; 1/2/3/4/5 (default is 2). 有改动
                                                    !!  1=one-point mutation, fixed rate.
                                                    !!  2=one-point, adjustable rate based on fitness.
                                                    !!  3=one-point, adjustable rate based on distance.
                                                    !!  4=one-point+creep, fixed rate.
                                                    !!  5=one-point+creep, adjustable rate based on fitness.
                                                    !!  6=one-point+creep, adjustable rate based on distance.
    real(wp),intent(in),optional       :: fdif      !! relative fitness differential; range from 0 选择父母时使用
                                                    !! (none) to 1 (maximum).  (default is 1.0)
    integer(I1B),intent(in),optional   :: irep      !! reproduction plan; 1/2/3=Full generational 有改动
                                                    !! replacement/Steady-state-replace-random/Steady-
                                                    !! state-replace-worst (default is 3)
    integer(I1B),intent(in),optional   :: ielite    !! elitism flag; 0/1=off/on (default is 0) 精英开关 有改动
                                                    !! (Applies only to reproduction plans 1 and 2)
    integer(I1B),intent(in),optional   :: ivrb      !! printed output 0/1/2=None/Minimal/Verbose  0无1最小3冗长 有改动
                                                    !! (default is 0)
    real(wp),intent(in),optional       :: convergence_tol    !! convergence tolerance; must be > 0.0 (default is 0.0001)
    integer,intent(in),optional        :: convergence_window !! convergence window; must be >= 0  有改动
                                                             !! This is the number of consecutive solutions
                                                             !! within the tolerance for convergence to
                                                             !! be declared (default is 20)
    real(wp),intent(in),optional       :: initial_guess_frac !! fraction of the initial population
                                                             !! to set equal to the initial guess.  Range from 0
                                                             !! (none) to 1.0 (all). (default is 0.1 or 10%).
    integer,intent(in),optional        :: iseed              !! random seed value; must be > 0 (default is 999)
    !call me%父类名%父类实例方法
    call me%pikaia_class%init(n,xl,xu,f,status,& !调用父类的构造函数 节省代码工作量
                            iter_f,&
                            np,ngen,nd,pcross,pmutmn,pmutmx,pmut,imut,&
                            fdif,irep,ielite,ivrb,&
                            convergence_tol,convergence_window,initial_guess_frac,&
                            iseed)
    !以下为子类独有实例变量赋值 依赖于父类实例变量
    me%nb=CEILING(me%nd/log10(2.0_wp),I1B) !十进制有效数字转换为二进制有效数字  有改动 二进制
    me%z=real(ibset(0_IB,me%nb)-1_IB,wp)   !有改动 二进制 2**me%nb-1

    !Print a header
    if (me%ivrb>0_I1B) then !打印表头
        write(output_unit,'(A,I4)')    '    Binary Length of Chromosome: ',me%nb     !有改动 二进制
    end if

    block !有改动 二进制
        integer(IB),pointer::a
        allocate(a)

        if (me%nb > bit_size(a) - 1) then !bit_size(a)=32 最高位为符号位 编码不得占用
            write(output_unit,'(A)') &
            ' ERROR: IB is too small, shoud be int64.'
            stop  !有改动
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

    subroutine pikaia(me,x,f,status) !x[0,1] 有改动 二进制 子类覆盖父类同名同接口方法

    implicit none

    !subroutine arguments:
    class(pikaia_binary_class),intent(inout) :: me
    real(wp),dimension(:),intent(inout)    :: x      !! Input - initial guess for solution vector. [0,1]
                                                     !! Output - the "fittest" (optimal) solution found, [0,1]
                                                     !! i.e., the solution which maximizes the fitness function.
    real(wp),intent(out)                   :: f      !! the (scalar) value of the fitness function at x
    integer(I1B),intent(out)               :: status !! an indicator of the success or failure 有改动
                                                     !! of the call to pikaia (0=success; non-zero=failure)

    !Local variables 
    !以下两行代码 父类与子类不同 子类二进制编码
    integer(IB),dimension(me%n)    :: gn1 !母种群染色体编码  有改动 二进制
    integer(IB),dimension(me%n)    :: gn2 !父种群染色体编码  有改动 二进制

    include "pikaia_sub.f90" !此段代码父类与子类共享
    end subroutine pikaia
!*****************************************************************************************

!*****************************************************************************************
!>
!  Encode phenotype parameters into integer genotype
!  ph(k) are x,y coordinates [ 0 < x,y < 1 ]

    subroutine encode_binary(me,ph,gn) !实数转换为二进制整数  有改动 二进制 子类覆盖父类同名同接口方法

    implicit none

    class(pikaia_binary_class),intent(inout)  :: me
    real(wp),dimension(me%n),intent(in)       :: ph  ![0,1]
    integer(IB),dimension(me%n),intent(out)   :: gn  !有改动 二进制

    gn=int(ph*me%z,IB)  !实数[0,1]转换为二进制编码
    if(me%IsGray)call binary2gray(gn,gn) !二进制编码转换为格雷编码

    end subroutine encode_binary
!*****************************************************************************************

!*****************************************************************************************
!>
!  decode genotype into phenotype parameters
!  ph(k) are x,y coordinates [ 0 < x,y < 1 ]

    subroutine decode_binary(me,gn,ph) !二进制整数转换为实数  有改动 二进制 子类覆盖父类同名同接口方法

    implicit none

    class(pikaia_binary_class),intent(inout) :: me
    integer(IB),dimension(me%n),intent(in)   :: gn !有改动 二进制
    real(wp),dimension(me%n),intent(out)     :: ph ![0,1]

    if(me%IsGray)then
        block
            integer(IB),dimension(me%n) :: gn0
            call gray2binary(gn,gn0) !格雷编码转换为二进制编码
            ph=real(gn0,wp)/me%z !有改动 二进制 二进制编码转换为实数[0,1]
        end block
    else
        ph=gn/me%z !有改动 二进制 二进制编码转换为实数[0,1]
    end if

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

    subroutine cross_binary(me,gn1,gn2) !交叉 有改动 二进制 子类覆盖父类同名同接口方法

    implicit none

    class(pikaia_binary_class),intent(inout)    :: me
    integer(IB),dimension(me%n),intent(inout)   :: gn1  !有改动 二进制
    integer(IB),dimension(me%n),intent(inout)   :: gn2  !有改动 二进制

    integer :: i, ispl, ispl2, itmp
    integer(IB) :: t=0_IB
    integer(I1B),parameter :: method = 1_I1B ! speed:method1 > method2 

    !Use crossover probability to decide whether a crossover occurs
    if (urand()<me%pcross) then
        select case(method)
        case(1_I1B)  !method1:all vectors have same ispl and ispl2.
            !Compute first crossover point
            ispl=int(urand()*me%nb)+1 ![1,me%nb] 交叉位置起点

            !Now choose between one-point and two-point crossover
            if (urand()<0.5_wp) then !交叉位置终点
                ispl2=me%nb !相当于单点交叉
            else
                ispl2=int(urand()*me%nb)+1 ![1,me%nb] 相当于两点交叉
                !Un-comment following line to enforce one-point crossover
                !ispl2=me%nb !强制单点交叉
                if (ispl2<ispl) then !交换以确保ispl1<=isp2
                    itmp=ispl2
                    ispl2=ispl
                    ispl=itmp
                end if
            end if

            itmp=ispl2-ispl+1 !避免重复计算 交叉基因长度
            !Swap genes from ispl to ispl2
            do i=1,me%n !每个自变量均交叉
                !CALL MVBITS(FROM, FROMPOS, LEN, TO, TOPOS)
                !将整数FROM的FROMPOS位到FROMPOS+LEN-1位信息复制到整数TO的TOPOS位到TOPOS+LEN-1 LEN为复制位的长度
                t=IBITS(gn2(i),ispl-1,itmp) !t储存gn2(i)基因信息
                CALL MVBITS(gn1(i), ispl-1, itmp, gn2(i), ispl-1) !gn2=gn1
                CALL MVBITS(t     , 0     , itmp, gn1(i), ispl-1) !gn1=t
            end do
        case(2_I1B)  !method2:the same method as surperclass
            !Compute first crossover point
            ispl=int(urand()*me%nb*me%n)+1 ![1,me%nb*me%n] 交叉位置起点

            !Now choose between one-point and two-point crossover
            if (urand()<0.5_wp) then !交叉位置终点
                ispl2=me%nb*me%n !相当于单点交叉
            else
                ispl2=int(urand()*me%nb*me%n)+1 ![1,me%nb*me%n] 相当于两点交叉
                !Un-comment following line to enforce one-point crossover
                !ispl2=me%nb*me%n !强制单点交叉
                if (ispl2<ispl) then !交换以确保ispl1<=isp2
                    itmp=ispl2
                    ispl2=ispl
                    ispl=itmp
                end if
            end if
            
            block
                integer,dimension(:),allocatable :: nn1,nn2 !自变量染色体中的起点、终点位置
                integer :: i1,i2
                i1=(ispl -1)/me%nb+1  ![1,me%n] !起点所在自变量  ispl  ispl2
                i2=(ispl2-1)/me%nb+1  ![1,me%n] !终点所在自变量   30     31
                ALLOCATE(nn1(i1:i2),nn2(i1:i2)) !分配内存         i1     i2
                !                                                  1      2
                if (i1/=i2)then !ispl2、ispl不在同一个自变量
                    do i=i1,i2
                        if (i==i1) then !最左边自变量
                            nn1(i)=mod(ispl -1,me%nb)+1
                            nn2(i)=me%nb
                        else if (i==i2) then !最右边自变量
                            nn1(i)=1
                            nn2(i)=mod(ispl2-1,me%nb)+1
                        else !中间自变量全部交叉
                            nn1(i)=1
                            nn2(i)=me%nb
                        end if
                        
                        itmp=nn2(i)-nn1(i)+1 !避免重复计算
                        t=IBITS(gn2(i),nn1(i)-1,itmp) !t储存gn2(i)基因信息
                        CALL MVBITS(gn1(i), nn1(i)-1, itmp, gn2(i), nn1(i)-1)  !gn2=gn1
                        CALL MVBITS(t     , 0       , itmp, gn1(i), nn1(i)-1)  !gn1=t
                    end do
                else !i1==i2  ispl2、ispl在同一个自变量
                    itmp=ispl2-ispl+1 !避免重复计算
                    nn1(i1)=mod(ispl -1,me%nb)+1 !自变量染色体中的起点位置
                    !nn2(i1)=mod(ispl2-1,me%nb)+1 !自变量染色体中的终点位置
                    nn2(i1)=nn1(i1)+itmp-1 !自变量染色体中的终点位置
                    
                    !itmp==ispl2-ispl+1==nn(2,i1)-nn(1,i1)+1              
                    t=IBITS(gn2(i1),nn1(i1)-1,itmp) !t储存gn2(i1)基因信息
                    CALL MVBITS(gn1(i1), nn1(i1)-1, itmp, gn2(i1), nn1(i1)-1)  !gn2=gn1
                    CALL MVBITS(t      , 0        , itmp, gn1(i1), nn1(i1)-1)  !gn1=t
                end if
                
                DEALLOCATE(nn1,nn2) !释放内存
            end block
        end select
    end if

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

    subroutine mutate_binary(me,gn) !变异 有改动 二进制 子类覆盖父类同名同接口方法

    implicit none

    class(pikaia_binary_class),intent(inout)    :: me
    integer(IB),dimension(me%n),intent(inout)   :: gn  !有改动 二进制

    integer(I1B) :: i,j,inc

    !Decide which type of mutation is to occur
    if (me%imut>=4_I1B .and. urand()<=0.5_wp) then !单点变异+蠕变

        !CREEP MUTATION OPERATOR 蠕变即染色体编码随机加1或减1 变异较小
        !Subject each locus to random +/- 1 increment at the rate pmut
        do i=1_I1B,me%n
            do j=1_I1B,me%nb
                if (urand()<me%pmut) then
                    if (me%IsGray) then !Gray
                        gn(i) = ieor(gn(i),ibset(0_IB,j-1))   !取反  
                        !another method shiftl(1_IB,j-1)
                    else !Binary
                        inc=nint( urand() )*2-1 !-1或1
                        gn(i)=gn(i)+ibset(0_IB,j-1)*inc   !ibset(0_IB,j-1) =2**(j-1)  
                        !another method shiftl(1_IB,j-1)
                        !限制范围
                        if (gn(i)<0_IB) gn(i)=0_IB
                        if (gn(i)>me%z) gn(i)=me%z
                    end if
                end if
            end do
        end do

    else !单点变异

        !UNIFORM MUTATION OPERATOR 统一变异操作 变异较大 可能出现0变7或7变0
        !Subject each locus to random mutation at the rate pmut
        do i=1_I1B,me%n !染色体每个位置
            do j=1_I1B,me%nb,3_I1B ! 三个位置分别为j j+1 j+2
                !nd 9  8  7  6  5  4 
                !nb 30 27 24 20 17 14
                if (urand()<me%pmut) then ! 4<=nd<=9
                    !int(urand()*4) ![0,3] or [00,11]
                    if(me%nd<=6_I1B)then !mod(me%nb,3)==2
                        if(j==me%nb-1_I1B)then !最后一个位置的位起点 19或16或13
                            !int32 最低位为0(右边第1位) 最高位为31(右边第32位)
                            !0为随机数位起点 2为位长度 j-1为gn(i)位起点
                            CALL MVBITS(int(urand()*4), 0, 2, gn(i), j-1)
                            exit !退出内循环
                        end if
                    end if
                    !int(urand()*8) ![0,7] or [000,111]
                    !0为随机数位起点 3为位长度 j-1为gn(i)位起点
                    CALL MVBITS(int(urand()*8), 0, 3, gn(i), j-1)
                end if
            end do
        end do

    end if

    end subroutine mutate_binary
!*****************************************************************************************
    !二进制编码转换为格雷码 适用于无符号整数  增加
    pure elemental subroutine binary2gray(b,g) !子类独有实用方法
    implicit none

    integer(kind=IB), intent(in)  :: b
    integer(kind=IB), intent(out) :: g
    !ieor位异或 同0异1 shiftr位右移
    g = ieor(b,shiftr(b,1))  ! 参考c语言x^(x>>1) F08
    return
    end subroutine binary2gray
    !格雷码转换为二进制编码 适用于无符号整数  增加
    pure elemental subroutine gray2binary(g,b) !子类独有实用方法
    implicit none

    integer(kind=IB), value       :: g ! 传值
    integer(kind=IB), intent(out) :: b

    b=g                 !unsigned int y = x;

    do while (g > 0_IB) !while(x>>=1) !x>>=1等价于x=x>>1 非0真0假
      g=shiftr(g,1)     !位右移
      b=ieor(b,g)       !y ^= x; 同0异1
    end do

    return
    end subroutine gray2binary

    subroutine SetIsGray(me,IsGray) !子类独有实例方法 新实例变量IsGray初始化
    implicit none

    class(pikaia_binary_class),intent(inout) :: me
    logical,intent(in)                       :: IsGray

    me%IsGray=IsGray !T Gray,F Binary
    !Print a header
    if (me%ivrb>0_I1B) then !打印表头
        write(output_unit,'(A,L4)')    'Chromosome Is Encoded with Gray: ',me%IsGray !有改动 二进制
    end if
    end subroutine SetIsGray
!*****************************************************************************************
    end module pikaia_binary_module
!*****************************************************************************************
