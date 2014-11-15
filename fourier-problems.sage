load('https://github.com/zspali/math-practice-problem-generator/raw/master/pc.sage')

var('a0,a_n,b_n')

flist = ["Random", "Piecewise Linear", "Other"]
llist = ["Random", "From Graph", "From Formula"]
ilist = ["Random", "$[-L,L]$", "$[0,L]$"]
rlist = ["Random", "Coefficients", "Series", "Partial Sum"]

def generate_problem():
    global per, f, s, m, p, pname, rname, lname, fdef

    pname = ["", "sine","cosine"][p]
    rname = ["Fourier {} coefficients".format(pname), "Fourier {} series".format(pname), "partial sum $s_{{{}}}(x)$ of the Fourier {} series".format(latex(m),pname)][rtype-1]
    lname = ["piecewise linear",""][ftype-1]

    if ftype == 1:
        f = generate_pc(max_step = max_step, max_abs = max_abs)
        if itype == 1:
            s = randint(0,2)
            if s == 0:
                f=f.translate(-f.L()/2)
                if f.is_even():
                    s = 2
                elif f.is_odd():
                    s = 1
            else:
                f=f.extension(s=s)
            per = randint(1,max_period)
            f=f.periods(per)
    if ftype == 2:
            L = randint(1,max_step)
            fun = choice([randint(1,max_abs)*(-1)^randint(0,1)*x^2,cos(per*pi*x/L),sin(per*pi*x/L),exp(per*pi*x/L)])
            f = pc([[[(-L)*(2-itype),L],fun]])
            per = 1

    if ftype == 2 or ltype == 2:
        fdef = "with formula {}".format(f.say_function())
    if ftype == 1 and ltype == 1:
        fdef = "with graph"

    m = randint(min_m,max_m)
    p = (itype-1)*randint(1,2)
    if itype == 2:
        s = p

def say_problem():
    html("Compute the {} of the {} function {} <br>".format(rname,lname,fdef))
    if ftype == 1 and ltype == 1:
        f.plot_function(thickness=2).show()
        html("<br>")

def say_solution():
    global f
    print s
    
    if itype == 1:
        f=f.slice((-f.L())*(1-min(s,1))/per,f.L()/per)
        say = "we can let $L={}$. <br>".format(latex(f.L()))
        if per > 1:
            say = "Since ${}$ has period {}, ".format(latex(f), latex(2*f.L())) + say
        html(say[0].upper()+say[1:])
    
    cosf = cos(n*pi*x/f.L())
    sinf = sin(n*pi*x/f.L())
    mult = (1+min(s,1))/f.L()
    cterm = ["",latex(a0/2)+"+"][(s-1) % 2]
    sterm = latex((a_n*cosf)*((s-1) % 2) + b_n*sinf*(min(1,2-s)))
        
    if rtype == 2:
        say = ""
        if itype == 1:
            if s == 1:
                say = "Since ${}$ is odd, ".format(latex(f))
            elif s == 2:
                say = "Since ${}$ is even, ".format(latex(f))
        say = say + r"the {} is ${}\sum_{{n=1}}^\infty {}$ <br>".format(rname, cterm, sterm)
        html(say[0].upper()+say[1:])
        
    if ftype == 1 and ltype == 1:
        html("The formula for ${}$ is {}".format(latex(f),f.say_function()))
        
    html("The Fourier {} coefficients are: <br>".format(pname))
    if s == 1 and rtype != 2 and itype == 1:
        html(r"Since ${}$ is odd, we have $$a_n=0,\quad n\ge0$$".format(latex(f)))
    if s != 1:
        html(f.say_integral(mult=mult, fsymbol="a_0"))
        html(f.say_integral(mult=mult, fsymbol="a_n", fmult = cosf))
    if s == 2 and rtype != 2 and itype == 1:
        html(r"Since ${}$ is even, we have $$b_n=0,\quad n>0$$".format(latex(f)))
    if s != 2:
        html(f.say_integral(mult=mult, fsymbol="b_n", fmult = sinf))
    if rtype == 3:
        html("The {} is".format(rname))
        if s > 0:
            f = f.extension(s=s)
        html("$$s_{{{}}}({})={}\sum_{{n=1}}^{{{}}} {}={}$$".format(latex(m),latex(f.fvar), cterm,latex(m), sterm,latex(f.partial_sum(m=m))))