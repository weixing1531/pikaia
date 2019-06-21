    abstract interface !抽象接口

        subroutine pikaia_func(me,x,f)  !适应度函数接口 具体实现见主程序
        !! The interface for the function that pikaia will be maximizing.
        import :: wp,pikaia_class !导入模块变量定义 pikaia_class为类名
        implicit none
        class(pikaia_class),intent(inout)  :: me    !! pikaia class
        real(wp),dimension(:),intent(in)   :: x     !! optimization variable vector
        real(wp),intent(out)               :: f     !! fitness value
        end subroutine pikaia_func

        subroutine iter_func(me,iter,x,f) !历次迭代结果 具体实现见主程序
        !! The interface for the function that user can specify
        !! to report each iteration when pikaia is running.
        !! The best (fittest) population member is passed to
        !! this routine in each generation.
        import :: wp,pikaia_class !导入模块变量定义 pikaia_class为类名
        implicit none
        class(pikaia_class),intent(inout)  :: me    !! pikaia class
        integer,intent(in)                 :: iter  !! iteration number
        real(wp),dimension(:),intent(in)   :: x     !! optimization variable vector
        real(wp),intent(in)                :: f     !! fitness value
        end subroutine iter_func

    end interface