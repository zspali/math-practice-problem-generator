var('n,x')
class pc:

    def __init__(self, flist, fsymbol = function("f",x)):

        self.flist = [[e[0],SR(e[1])] for e in flist]
        self.fsymbol = fsymbol
        self.fvar = [fvar for fvar in fsymbol.args() if fvar != n][0]

    def _latex_(self):
        return latex(self.fsymbol)
    
    def split_flist(self, x0 = 0, side = 1):
        old_fl = list(self.flist)
        left_fl = []
        right_fl=[]
        for e in old_fl:
            if e[0][1] <= x0:
                left_fl.append(e)
            elif e[0][0] < x0:
                left_fl.append([[e[0][0],x0],e[1]])
                right_fl.append([[x0,e[0][1]],e[1]])
            else:
                right_fl.append(e)
        return [pc(left_fl),pc(right_fl)][side]

    def optimized_flist(self):
        old_fl = list(self.flist)
        new_fl = []
        c1 = old_fl.pop(0)
        while not old_fl == []:
            c2 = old_fl.pop(0)
            if c1[1] == c2[1]:
                c1[0] = (c1[0][0],c2[0][1])
            else:
                new_fl.append(c1)
                c1 = c2
        new_fl.append(c1)
        return pc(new_fl)
    
    def translate(self, delta):
        return pc([[[e[0][0]+delta,e[0][1]+delta],e[1].substitute(self.fvar == self.fvar-delta)] for e in self.flist])
    
    def periods(self, n_periods, centered = true):
        f = pc(reduce(lambda x,y: x+y, [self.translate(k*2*self.L()).flist for k in range(n_periods)]))
        if centered:
            f = f.translate(-(n_periods-1)*self.L())
        return f

    def L(self):
        return self.flist[-1][0][1]

    def lb(self):
        return self.flist[0][0][0]

    def sin_coeff(self, m = n):
        return 1/self.L()*(sum([integrate(c[1]*sin(m*pi*self.fvar/self.L()),self.fvar,c[0][0],c[0][1]) for c in self.flist])-self.fvar+self.fvar).full_simplify().expand().subs_expr(sin(n*pi) == 0)

    def cos_coeff(self, m = n):
        return 1/self.L()*(sum([integrate(c[1]*cos(m*pi*self.fvar/self.L()),self.fvar,c[0][0],c[0][1]) for c in self.flist])-self.fvar+self.fvar).full_simplify().expand().subs_expr(sin(n*pi) == 0)

    def say_function(self):

        if len(self.flist) == 1:
            return r"{}".format(latex(self.flist[0][1]))
        else:
            return r"\begin{cases} " + reduce(lambda x,y: x + r" \\ " + y, [r"{} & {} \le ".format(latex(c[1]),latex(c[0][0])) + latex(self.fvar) + " < {}".format(latex(c[0][1])) for c in self.flist]) + r" \end{cases}"

    def say_integral(self, lb, ub, mult):

        fl_comp = list(self.flist)

        while fl_comp[0][0][1] <= lb:
            fl_comp.pop(0)

        while fl_comp[-1][0][0] >= ub:
            fl_comp.pop()

        fl_comp[0][0] = (lb,fl_comp[0][0][1])
        fl_comp[-1][0] = (fl_comp[-1][0][0],ub)

        say = "$" +latex(mult) + r"\left(" + r"\int_{0}^{1}{2}\,\mathrm d{3}".format(blatex(lb),blatex(ub),latex(self.fsymbol),latex(self.fvar))  + r"\right)"
        fl_comp = [c for c in fl_comp if c[1] != 0]



        if fl_comp == []:
            return say +r"=0$ <br>"
        else:
            say += "$ <br> $=" + latex(mult) + r"\left(" + reduce(lambda x,y: x + "+" + y, [r"\int_{}^{}{}\,\mathrm d{}".format(blatex(c[0][0]),blatex(c[0][1]),blatex(c[1]),latex(self.fvar)) for c in fl_comp]) + r"\right)"
            say += "$ <br> $=" + latex(mult) + r"\left(" + reduce(lambda x,y: x + "+" + y, [r"\left[{}\right]_{}^{}".format(latex(c[1]),blatex(c[0][0]),blatex(c[0][1])) for c in [[c[0],(integrate(c[1],self.fvar)-self.fvar+self.fvar).full_simplify().expand()] for c in fl_comp]]) + r"\right)"
            return say + r"$ <br> $={}$ <br>".format(latex((sum([integrate(c[1]*mult,self.fvar,c[0][0],c[0][1]) for c in fl_comp])-self.fvar+self.fvar).full_simplify().expand().subs_expr(sin(n*pi) == 0)))

    def extension(self, s = 0):
        split_fl = self.split_flist().flist
        mirror = [[(-c[0][1],-c[0][0]),(-1)^s*(c[1]).subs(self.fvar==-self.fvar)] for c in split_fl]
        mirror.reverse()
        new_fl = mirror + split_fl

        return pc(new_fl)

    def fval(self, ival):

        if ival == self.flist[0][0][0] or ival == self.flist[-1][0][1]:
            return 0
        else:
            fiter = iter(self.flist)
            c = fiter.next()
            while c[0][1] < ival:
                c = fiter.next()
            if ival == c[0][1]:
                return (c[1] + fiter.next()[1]+self.fvar-self.fvar).subs(self.fvar==ival)/2
            else:
                return (c[1]-self.fvar+self.fvar).subs(self.fvar==ival)

    def say_sin_coeff(self):

        return r"We have $b_n=$ " + pc([[c[0],c[1]*sin(n*pi*self.fvar/self.L())] for c in self.flist], fsymbol = self.fsymbol * sin(n*pi*self.fvar/self.L())).say_integral(0,self.L(),2/self.L()) + "<br>"

    def say_cos_coeff(self, cterm = True):

        say = r"We have "
        if cterm:
            say += "$a_0=$ " + self.say_integral(0,self.L(),2/self.L()) + ", <br> and for $n>0:$ "
        return say + "$a_n=$ " +pc([[c[0],c[1]*cos(n*pi*self.fvar/self.L())] for c in self.flist], fsymbol = self.fsymbol * cos(n*pi*self.fvar/self.L())).say_integral(0,self.L(),2/self.L())

    def plot_function(self, color='blue', thickness=1, legend_label=""):

        return sum([plot(c[1],(self.fvar,c[0][0],c[0][1]), color=color, thickness=thickness) for c in self.flist[:-1]]+[plot(c[1],(self.fvar,c[0][0],c[0][1]), color=color, thickness=thickness, legend_label=legend_label) for c in self.flist[-1:]])
