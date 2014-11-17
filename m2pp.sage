min_m = 3
max_m = 5

max_step = 4
max_abs = 4
max_period = 3

var('a0,a_n,b_n')

flist = ["Random", "Piecewise Linear", "Other"]
llist = ["Random", "From Graph", "From Formula"]
ilist = ["Random", "$[-L,L]$", "$[0,L]$"]
rlist = ["Random", "Coefficients", "Series", "Partial Sum"]
plist = ["Fourier", "Inhomogeneous Equation"]
elist = ["Random", "Real Distinct", "Complex Conjugate Pair", "Repeated"]

@interact
def _f(psel = Selector(plist, label = 'Problem Topic:', selector_type='button')):
    if psel == plist[0]:
        @interact
        def _f(fsel = Selector(flist, label = 'Function type:', selector_type = 'button'), isel = Selector(ilist, label = 'Interval:', selector_type = 'button'), rsel = Selector(rlist, label = 'Question:', selector_type = 'button'), regen = Button(text="Regenerate Problem", default=True, value=True, label = "")):

            ftype = flist.index(fsel)
            itype = ilist.index(isel)
            rtype = rlist.index(rsel)

            if ftype == 0:
                ftype = randint(1,len(flist)-1)
            if rtype == 0:
                rtype = randint(1,len(rlist)-1)
            if itype == 0:
                itype = randint(1,len(ilist)-1)

            if ftype == 1:
                @interact
                def _f(lsel = Selector(llist, label = 'Function given:', selector_type = 'button'), regen1 = Button(text="Regenerate Problem", default=True, value=True, label = "")):
                    ltype = llist.index(lsel)

                    if ltype == 0:
                        ltype = randint(1,len(llist)-1)
                    generate_problem(ftype, ltype, itype, rtype)
                    say_problem(ftype, ltype, itype, rtype)
                    @interact
                    def _f(solution = Checkbox(label = "Show Solution", default = False),plots = Checkbox(label = "Draw Plots", default = False)):
                        if solution:
                            say_solution(ftype, ltype, itype, rtype)
                        if plots:
                            @interact
                            def _f(plot_m=(5,(1..50))):
                                if f.flist[0][0][0] == 0:
                                    show(plot(f.extension(s=s).plot_function(thickness=2)+f.extension(s=s).plot_psum(plot_m,thickness=2)))
                                else:
                                    show(plot(f.plot_function(thickness=2)+f.plot_psum(plot_m,thickness=2)))

            else:
                ltype = 0
                generate_problem(ftype, ltype, itype, rtype)
                say_problem(ftype, ltype, itype, rtype)
                @interact
                def _f(solution = Checkbox(label = "Show Solution", default = False),plots = Checkbox(label = "Draw Plots", default = False)):
                    if solution:
                        say_solution(ftype, ltype, itype, rtype)
                    if plots:
                        @interact
                        def _f(plot_m=(5,(1..50))):
                            if f.flist[0][0][0] == 0:
                                show(plot(f.extension(s=s).plot_function(thickness=2)+f.extension(s=s).plot_psum(plot_m,thickness=2)))
                            else:
                                show(plot(f.plot_function(thickness=2)+f.plot_psum(plot_m,thickness=2)))
    else:
        @interact
        def _f(esel = Selector(elist, label = 'Eigenvalues:', selector_type='button'), regen2 = Button(text="Regenerate Problem", default=True, value=True, label = "")):
            var('alpha,t,u,x,y,x1,x2,c1,c2,y1,y2,s');x1f=function('x1f',t); x2f=function('x2f',t); cv = vector([c1,c2]); yx = vector([y1,y2])
            etype = elist.index(esel)
            A = HDE2d(generate_mx(evtype = etype))
            A = IHDE2d(A.A, A.generate_g())
            html(r"Find a particular solution of the inhomogeneous equation $$\mathbf x'(t)={A}\mathbf x(t)+{g}$$".format(A=latex(A.A), g=latex(A.g.column())))
            @interact
            def _f(solution = Checkbox(label = "Show Solution", default = False)):
                if solution:
                    if A.disc() < 0:
                        html(A.say_u())
                    else:
                        html(A.say_y())