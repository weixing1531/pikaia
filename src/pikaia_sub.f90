    integer  :: k,ip,ig,ip1,ip2,new,newtot,istart,i_window
    real(wp) :: current_best_f, last_best_f, fguess
    logical  :: convergence
    real(wp),dimension(me%n,2)     :: ph    !��ĸ��Ⱥ�Ա��� [0,1]
    real(wp),dimension(me%n,me%np) :: oldph !����Ⱥ  [0,1]
    real(wp),dimension(me%n,me%np) :: newph !����Ⱥ  [0,1]
    integer,dimension(me%np)       :: ifit !ifit�����±�
    integer,dimension(me%np)       :: jfit !jfit��������
    real(wp),dimension(me%np)      :: fitns !��Ⱥ��Ӧ��
    real(wp),dimension(me%n)       :: xguess !��Ⱥ�Ա����²�ֵ

    real(wp),parameter :: big = huge(1.0_wp)    !! a large number

    !initialize:
    call rninit(me%iseed)  !���ӳ�ʼ��
    me%bestft   = -big     !����Ϊ��Сֵ
    me%pmutpv   = -big     !����Ϊ��Сֵ
    me%pmut     = me%pmuti !set initial mutation rate (it can change)
    i_window    = 0
    last_best_f = -big     !����Ϊ��Сֵ
    convergence = .false.  !������
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
        if (any(me%imut==[2,3,5,6])) call me%adjmut(oldph,fitns,ifit) !�иĶ�

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
