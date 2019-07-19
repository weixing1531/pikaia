    integer  :: k,ip,ig,ip1,ip2,new,newtot,istart,i_window
    real(wp) :: current_best_f, last_best_f, fguess
    logical  :: convergence
    real(wp),dimension(me%n,2)     :: ph    !父母种群自变量 [0,1]
    real(wp),dimension(me%n,me%np) :: oldph !旧种群  [0,1]
    real(wp),dimension(me%n,me%np) :: newph !新种群  [0,1]
    integer,dimension(me%np)       :: ifit !ifit升序下标
    integer,dimension(me%np)       :: jfit !jfit升序排名
    real(wp),dimension(me%np)      :: fitns !种群适应度
    real(wp),dimension(me%n)       :: xguess !种群自变量猜测值

    real(wp),parameter :: big = huge(1.0_wp)    !! a large number

    !initialize:
    call rninit(me%iseed)  !种子初始化
    me%bestft   = -big     !先设为极小值
    me%pmutpv   = -big     !先设为极小值
    me%pmut     = me%pmuti !set initial mutation rate (it can change)
    i_window    = 0
    last_best_f = -big     !先设为极小值
    convergence = .false.  !不收敛
    status      = 0

    !Handle the initial guess:
    if (me%initial_guess_frac==0.0_wp) then !此时猜测值没用 种群全部采用随机值

        !initial guess not used (totally random population)

        istart = 1  !index to start random population members 随机种群序号起点

    else

        !use the initial guess:

        xguess = x !自变量猜测值[0,1]
        do k=1,me%n    !make sure they are all within the [0,1] bounds
            xguess(k) = max( 0.0_wp, min(1.0_wp,xguess(k)) ) !猜测值边界检查
        end do
        call me%ff(xguess,fguess) !计算猜测种群适应度函数

        !how many elements in the population to set to xguess?:
        ! [at least 1, at most n]
        istart = max(1, min(me%np, int(me%np * me%initial_guess_frac)))

        do k=1,istart !前istart个种群采用猜测值
            oldph(:,k) = xguess
            fitns(k)   = fguess
        end do

        istart = istart + 1  !index to start random population members

    end if
    !前istart个种群之外的采用随机值
    !Compute initial (random but bounded) phenotypes
    do ip=istart,me%np
        do k=1,me%n
            oldph(k,ip)=urand()  !from [0,1]
        end do
        call me%ff(oldph(:,ip),fitns(ip)) !计算随机种群适应度函数
    end do

    !Rank initial population by fitness order 初始种群适应度排名
    call me%rnkpop(fitns,ifit,jfit)  !ifit升序下标 jfit升序排名

    !Main Generation Loop
    do ig=1,me%ngen !迭代

        !Main Population Loop
        newtot=0
        do ip=1,me%np/2 !一半种群

            !1. pick two parents 选择父母
            call me%select_parents(jfit,ip1,ip2) !jfit升序排名
            !ip1母种群  ip2父种群
            !2. encode parent phenotypes
            call me%encode(oldph(:,ip1),gn1) !实数转换为十进制整数
            call me%encode(oldph(:,ip2),gn2) !实数转换为十进制整数

            !3. breed 培育
            call me%cross(gn1,gn2) !交叉
            call me%mutate(gn1) !母种群染色体变异
            call me%mutate(gn2) !父种群染色体变异

            !4. decode offspring genotypes
            call me%decode(gn1,ph(:,1)) !十进制整数转换为实数
            call me%decode(gn2,ph(:,2)) !十进制整数转换为实数

            !5. insert into population
            if (me%irep==1) then !繁殖模式 Full generational replacement
                call me%genrep(ip,ph,newph)
            else  !ifit升序下标 jfit升序排名
                call me%stdrep(ph,oldph,fitns,ifit,jfit,new) !更新种群oldph 更新排名
                newtot = newtot+new !newtot新种群中新个体数目
            end if

        end do    !End of Main Population Loop

        !if running full generational replacement: swap populations
        if (me%irep==1) call me%newpop(oldph,newph,ifit,jfit,fitns,newtot) !繁殖模式 Full generational replacement
        !newtot新种群中新个体数目 更新了种群适应度排名 此时oldph已经是更新后种群 oldph=newph
        !adjust mutation rate? 变异率动态调整
        if (any(me%imut==[2,3,5,6])) call me%adjmut(oldph,fitns,ifit) !有改动

        !report this iteration:
        if (me%ivrb>0) call me%report(oldph,fitns,ifit,ig,newtot) !详细打印迭代信息

        !report (unscaled) x:
        if (associated(me%iter_f)) & !过程指针me%iter_f有指向
            call me%iter_f(ig,me%xl+me%del*oldph(1:me%n,ifit(me%np)),fitns(ifit(me%np))) !打印每次迭代最强种群及其适应度

        !JW additions: add a convergence criteria
        ! [stop if the last convergence_window iterations are all within the convergence_tol]
        current_best_f = fitns(ifit(me%np))    !current iteration best fitness
        if (abs(current_best_f-last_best_f)<=me%convergence_tol) then
            !this solution is within the tol from the previous one
            i_window = i_window + 1    !number of solutions within the convergence tolerance
        else
            i_window = 0    !a significantly better solution was found, reset window
        end if
        !连续convergence_window次迭代后 后前两次迭代适应度相对误差小于convergence_tol即认为已收敛
        if (i_window>=me%convergence_window) then
            convergence = .true.
            exit !exit main loop -> convergence
        end if
        last_best_f = current_best_f    !to compare with next iteration 上一次迭代最强适应度

    end do    !End of Main Generation Loop

    !JW additions:
    if (me%ivrb>0) then !详细打印迭代信息
        if (convergence) then
            write(output_unit,'(A)') 'Solution Converged'
        else
            write(output_unit,'(A)') 'Iteration Limit Reached'
        end if
    end if

    !Return best phenotype and its fitness
    x = oldph(1:me%n,ifit(me%np)) !最强适应度对应自变量[0,1]
    f = fitns(ifit(me%np)) !最强适应度
