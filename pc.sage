var('n,x')
class pc:

    def __init__(self, flist, fsymbol = function("f",x)):

        self.flist = flist
        self.fsymbol = fsymbol
        self.fvar = [fvar for fvar in fsymbol.args() if fvar != n][0]

    def _latex_(self):
        return latex(self.fsymbol)

    def optimize_flist(self):
        old_fl = list(self.flist_init)
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

        self.flist = list(new_fl)
        return self.flist

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

    def odd_ext(self, n_periods = 1, s = 0):

        new_fl = list(self.flist)
        mirror = [[(-c[0][1],-c[0][0]),(-1)*(c[1]-self.fvar+self.fvar).subs(self.fvar==-self.fvar)] for c in new_fl]
        mirror.reverse()
        new_fl = mirror + new_fl

        step = new_fl[0][0][1] - new_fl[0][0][0]

        l = len(new_fl);
        for i in range(1, n_periods):
            for j in range(0,l):
                c = new_fl[j]
                new_fl.append([(c[0][0]+i*l*step,c[0][1]+i*l*step),(c[1]-self.fvar+self.fvar).subs(self.fvar==self.fvar-i*l*step)])

        new_fl = [[(c[0][0]+s,c[0][1]+s),(c[1]-self.fvar+self.fvar).subs(self.fvar==self.fvar-s)] for c in new_fl]

        return new_fl

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