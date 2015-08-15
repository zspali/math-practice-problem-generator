min_m = 3
max_m = 5

max_step = 4
max_abs = 4
max_period = 3

var('a0,a_n,b_n')

load('https://github.com/zspali/math-practice-problem-generator/raw/master/pc.sage')
load('https://github.com/zspali/math-practice-problem-generator/raw/master/fourier-problems.sage')
load('https://github.com/zspali/math-practice-problem-generator/raw/master/pdes-with-fourier.sage')

plist = ["Fourier", "PDE"]


@interact
def _f(psel = Selector(plist, label = 'Problem Topic:', selector_type='button')):
    flist = ["Random", "Piecewise Linear", "Other"]
    llist = ["Random", "From Graph", "From Formula"]
    ilist = ["Random", "$[-L,L]$", "$[0,L]$"]
    rlist = ["Random", "Coefficients", "Series", "Partial Sum"]

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
        min_m = 3
        max_m = 5

        max_step = 4
        max_abs = 4

        p2list = ["Random", "Heat Equation", "Wave Equation", "Laplace Equation"]

        @interact
        def _f(psel = Selector(p2list, label = "PDE type:", selector_type = "button"), regen = Button(text="Regenerate Problem", default=True, value=True, label = "")):
            if psel == p2list[0]:
                psel = choice(plist[1:])
            if psel == p2list[1]:

                blist = ["Random", "Temperatures fixed at ends", "Insulated ends"]
                rlist = ["Random", "Formal Solution", "Partial Sum"]
                llist = ["Random", "From Graph", "From Formula"]

                @interact
                def _f(bsel = Selector(blist, label = "BC type:", selector_type = "button"), rsel = Selector(rlist, label = "problem type:", selector_type = "button"), lsel = Selector(llist, label = "NHBC represented by:", selector_type = "button"), regen1 = Button(text="Regenerate Problem", default=True, value=True, label = "")):
                    if bsel == blist[0]:
                        bsel = choice(blist[1:])
                    if rsel == rlist[0]:
                        rsel = choice(rlist[1:])
                    if lsel == llist[0]:
                        lsel = choice(llist[1:])

                    f = generate_pc()
                    alpha2 = 2^randint(-3,3)

                    bc = [randint(-max_abs,max_abs), randint(-max_abs,max_abs)]

                    if bsel == blist[1]:


                        if bc == [0,0]:
                            problem = HHE1d(f = f, bc = bc, alpha2 = alpha2)
                        else:
                            problem = NHHE1d(f = f, bc = bc, alpha2 = alpha2)

                        if rsel == rlist[1]:
                            html(r"Find the formal solution to the problem")

                        else:
                            m = randint(min_m, max_m)
                            html(r"Find the partial sum $s_{{{m}}}(t,x)$ to the problem".format(m = latex(m)))

                        html(problem.say_eqs())

                        if lsel == llist[1]:
                            html(r"Where $f(x)$ has graph<br>")
                            show(f.plot_function())
                        else:
                            html(r"Where $f(x)$ has formula" + f.say_function())

                        @interact
                        def _f(solution = Checkbox(label = "Show Solution", default = False), cplot = Checkbox(label = r"Draw Contour Plot of $s_{50}(t,x)$ ", default = False)):
                            if solution:
                                html(problem.say_fseries())
                                if rsel == rlist[2]:
                                    html(problem.say_psum(m))

                            if cplot:
                                show(problem.plot_psum(50))

                    if bsel == blist[2]:

                        bc = [0,0]

                        problem = IHE1d(f = f, alpha2 = alpha2)

                        if rsel == rlist[1]:
                            html(r"Find the formal solution to the problem")

                        else:
                            m = randint(min_m, max_m)
                            html(r"Find the partial sum $s_{{{m}}}(t,x)$ to the problem".format(m = latex(m)))

                        html(problem.say_eqs())

                        if lsel == llist[1]:
                            html(r"Where $f(x)$ has graph<br>")
                            show(f.plot_function())
                        else:
                            html(r"Where $f(x)$ has formula" + f.say_function())

                        @interact
                        def _f(solution = Checkbox(label = "Show Solution", default = False), cplot = Checkbox(label = r"Draw Contour Plot of $s_{50}(t,x)$ ", default = False)):
                            if solution:
                                html(problem.say_fseries())
                                if rsel == rlist[2]:
                                    html(problem.say_psum(m))

                            if cplot:
                                show(problem.plot_psum(50))

            if psel == p2list[2]:

                blist = ["Random", "Zero Initial Velocity", "Zero Initial Displacement"]
                rlist = ["Random", "Formal Solution", "Partial Sum"]
                llist = ["Random", "From Graph", "From Formula"]

                @interact
                def _f(bsel = Selector(blist, label = "BC type:", selector_type = "button"), rsel = Selector(rlist, label = "problem type:", selector_type = "button"), lsel = Selector(llist, label = "NHBC represented by:", selector_type = "button"), regen1 = Button(text="Regenerate Problem", default=True, value=True, label = "")):
                    if bsel == blist[0]:
                        bsel = choice(blist[1:])
                    if rsel == rlist[0]:
                        rsel = choice(rlist[1:])
                    if lsel == llist[0]:
                        lsel = choice(llist[1:])

                    f = generate_pc()
                    a = 2^randint(-3,3)

                    if bsel == blist[1]:
                        problem = ZVWE1d(f = f, a = a)
                    else:
                        problem = ZDWE1d(f = f, a = a)

                    if rsel == rlist[1]:
                        html(r"Find the formal solution to the problem")

                    else:
                        m = randint(min_m, max_m)
                        html(r"Find the partial sum $s_{{{m}}}(t,x)$ to the problem".format(m = latex(m)))

                    html(problem.say_eqs())

                    if lsel == llist[1]:
                        html(r"Where $f(x)$ has graph<br>")
                        show(f.plot_function())
                    else:
                        html(r"Where $f(x)$ has formula" + f.say_function())

                    @interact
                    def _f(solution = Checkbox(label = "Show Solution", default = False), cplot = Checkbox(label = r"Draw Contour Plot of $s_{50}(t,x)$ ", default = False)):
                        if solution:
                            html(problem.say_fseries())
                            if rsel == rlist[2]:
                                html(problem.say_psum(m))

                        if cplot:
                            show(problem.plot_psum(50))

            if psel == p2list[3]:

                blist = ["Random", "Dirichlet Problem", "Neumann Problem"]
                rlist = ["Random", "Formal Solution", "Partial Sum"]
                llist = ["Random", "From Graph", "From Formula"]

                @interact
                def _f(bsel = Selector(blist, label = "BC type:", selector_type = "button"), rsel = Selector(rlist, label = "problem type:", selector_type = "button"), lsel = Selector(llist, label = "NHBC represented by:", selector_type = "button"), regen1 = Button(text="Regenerate Problem", default=True, value=True, label = "")):
                    if bsel == blist[0]:
                        bsel = choice(blist[1:])
                    if rsel == rlist[0]:
                        rsel = choice(rlist[1:])
                    if lsel == llist[0]:
                        lsel = choice(llist[1:])

                    xin = randint(0,1)
                    vari = [x,y][xin]

                    zero_int = bool( bsel == blist[2] )
                    f = generate_pc(vari=vari, zero_int = zero_int)

                    a = f.L()
                    b = randint(1,max_abs)

                    if randint(0,1) == 0:
                        bc = [[0,f]]
                    else:
                        bc = [[f,0]]

                    bc.insert(xin,[0,0])

                    if xin == 0:
                        ab = [a,b]
                    else:
                        ab = [b,a]

                    if bsel == blist[1]:
                        problem = DP2d(bc = bc, ab = ab)
                    else:
                        problem = NP2d(bc = bc, ab = ab)

                    if rsel == rlist[1]:
                        html(r"Find the formal solution to the problem")

                    else:
                        m = randint(min_m, max_m)
                        html(r"Find the partial sum $s_{{{m}}}(x,y)$ to the problem".format(m = latex(m)))

                    html(problem.say_eqs())

                    if lsel == llist[1]:
                        html(r"Where ${f}$ has graph<br>".format(f=latex(f)))
                        show(f.plot_function())
                    else:
                        html(r"Where ${f}$ has formula ".format(f=latex(f)) + f.say_function())

                    @interact
                    def _f(solution = Checkbox(label = "Show Solution", default = False), cplot = Checkbox(label = r"Draw Contour Plot of $s_{50}(x,y)$ ", default = False)):
                        if solution:
                            html(problem.say_fseries())
                            if rsel == rlist[2]:
                                html(problem.say_psum(m))

                        if cplot:
                            show(problem.plot_psum(50))      
