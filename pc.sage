var('n,x')

def blatex(toLatex):
    return '{' + latex(toLatex) + '}'

def latex_mult_sum(exp_list, mult = 1):
    subexpr = reduce(lambda x,y: x + "+" + y, exp_list)
    if mult == 0:
        return latex(0)
    elif mult == 1:
        return subexpr
    elif len(exp_list) == 1:
        return latex(mult) + subexpr
    else:
        return latex(mult) + r"\left(" + subexpr + r"\right)"

class pc:

    def __init__(self, flist, fsymbol = function("f",x)):

        self.flist = [[e[0],SR(e[1])] for e in flist]
        self.fsymbol = fsymbol
        self.fvar = [fvar for fvar in fsymbol.args() if fvar != n][0]

    def _latex_(self):
        return latex(self.fsymbol)
    
    def changed_flist(self, flist):
        return pc(flist, fsymbol = self.fsymbol)
    
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
        return [self.changed_flist(left_fl),self.changed_flist(right_fl)][side]
    
    def slice(self, x0, x1):
        old_fl = list(self.flist)
        new_fl=[]
        if x0 < x1:
            for e in old_fl:
                if e[0][0] < x1 and e[0][1] > x0:
                    new_fl.append([[max(e[0][0],x0),min(e[0][1],x1)],e[1]])
        return self.changed_flist(new_fl)

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
        return self.changed_flist(new_fl)
    
    def translate(self, delta):
        return self.changed_flist([[[e[0][0]+delta,e[0][1]+delta],e[1].substitute(self.fvar == self.fvar-delta)] for e in self.flist])
    
    def periods(self, n_periods, centered = true):
        f = self.changed_flist(reduce(lambda x,y: x+y, [self.translate(k*2*self.L()).flist for k in range(n_periods)]))
        if centered:
            f = f.translate(-(n_periods-1)*self.L())
        return f

    def L(self):
        return self.flist[-1][0][1]

    def lb(self):
        return self.flist[0][0][0]

    def is_even(self):
        flist = self.split_flist(side=0).flist+self.split_flist().flist
        return not False in [bool(SR(flist[len(flist)/2 + k][1])==SR(flist[len(flist)/2 - k - 1][1](x=-x))) for k in range(0,len(flist)/2)]
    
    def is_odd(self):
        flist = self.split_flist(side=0).flist+self.split_flist().flist
        return not False in [bool(SR(flist[len(flist)/2 + k][1])==-SR(flist[len(flist)/2 - k - 1][1](x=-x))) for k in range(0,len(flist)/2)]
    
    def sin_coeff(self, m = n):
        if self.flist[0][0][0]==0:
            d=2
        else:
            d=1
        return d/self.L()*(sum([integrate(c[1]*sin(m*pi*self.fvar/self.L()),self.fvar,c[0][0],c[0][1]) for c in self.flist])-self.fvar+self.fvar).full_simplify().expand().subs_expr(sin(n*pi) == 0)

    def cos_coeff(self, m = n):
        if self.flist[0][0][0]==0:
            d=2
        else:
            d=1
        return d/self.L()*(sum([integrate(c[1]*cos(m*pi*self.fvar/self.L()),self.fvar,c[0][0],c[0][1]) for c in self.flist])-self.fvar+self.fvar).full_simplify().expand().subs_expr(sin(n*pi) == 0)
    
    def partial_sum(self, m):
        x = self.fvar
        return self.cos_coeff(m=0)/2 + sum([SR(self.cos_coeff())(n=k)*cos(k*pi*x/self.L())+SR(self.sin_coeff())(n=k)*sin(k*pi*x/self.L()) for k in range(1,m+1)])

    def say_function(self):
        
        flist = self.optimized_flist().flist

        if len(flist) == 1:
            return r"$${}={},\quad {}\le {}< {}$$".format(latex(self),latex(flist[0][1]),latex(flist[0][0][0]),latex(self.fvar),latex(flist[0][0][1]))
        else:
            return r"$${}=\begin{{cases}} ".format(latex(self)) + reduce(lambda x,y: x + r" \\ " + y, [r"{} & {} \le ".format(latex(c[1]),latex(c[0][0])) + latex(self.fvar) + " < {}".format(latex(c[0][1])) for c in flist]) + r" \end{cases}$$"

    def say_integral(self, lb=0, ub=0, sb = true, mult = 1, fsymbol = "", fmult = 1):

        if sb:
            lb = self.flist[0][0][0]
            ub = self.flist[-1][0][1]
        
        fl_comp = list(self.optimized_flist().flist)

        while fl_comp[0][0][1] <= lb:
            fl_comp.pop(0)

        while fl_comp[-1][0][0] >= ub:
            fl_comp.pop()

        fl_comp[0][0] = (lb,fl_comp[0][0][1])
        fl_comp[-1][0] = (fl_comp[-1][0][0],ub)
        
        fl_comp = [[e[0],e[1]*fmult] for e in fl_comp]
        
        if fmult == 1:
            flatex = latex(self)
        else:
            flatex = latex(self) + latex(fmult)
        
        say = r"$$"
        if fsymbol != "":
            say += "{}=".format(fsymbol)
        say += latex_mult_sum([r"\int_{{{}}}^{{{}}}{}\,\mathrm d{}".format(latex(lb),latex(ub),flatex,latex(self.fvar))], mult=mult) + "="
        if fl_comp == []:
            return say +r"0 $$"
        else:
            say += latex_mult_sum([r"\int_{}^{}{}\,\mathrm d{}".format(blatex(c[0][0]),blatex(c[0][1]),blatex(c[1]),latex(self.fvar)) for c in fl_comp], mult=mult) + "$$"
            say += "$$=" + latex_mult_sum([r"\left[{}\right]_{}^{}".format(latex(c[1]),blatex(c[0][0]),blatex(c[0][1])) for c in [[c[0],SR(integrate(c[1],self.fvar)).full_simplify().expand()] for c in fl_comp]], mult=mult) + "$$"
            say += r"$$={}$$".format(latex(mult*sum([SR(integrate(c[1],self.fvar,c[0][0],c[0][1])).subs_expr(sin(n*pi) == 0) for c in fl_comp])))
        return say

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

        return r"We have $b_n=$ " + pc([[c[0],c[1]*sin(n*pi*self.fvar/self.L())] for c in self.flist], fsymbol = self.fsymbol * sin(n*pi*self.fvar/self.L())).say_integral(lb=0,ub=self.L(),mult=2/self.L()) + "<br>"

    def say_cos_coeff(self, cterm = True):

        say = r"We have "
        if cterm:
            say += "$a_0=$ " + self.say_integral(lb=0,ub=self.L(),mult=2/self.L()) + ", <br> and for $n>0:$ "
        return say + "$a_n=$ " +pc([[c[0],c[1]*cos(n*pi*self.fvar/self.L())] for c in self.flist], fsymbol = self.fsymbol * cos(n*pi*self.fvar/self.L())).say_integral(lb=0,ub=self.L(),mult=2/self.L())

    def plot_function(self, color='blue', thickness=1, legend_label="", default_label=true):
        if default_label:
            legend_label=r"${}$".format(latex(self))
        return sum([plot(c[1],(self.fvar,c[0][0],c[0][1]), color=color, thickness=thickness) for c in self.flist[:-1]]+[plot(c[1],(self.fvar,c[0][0],c[0][1]), color=color, thickness=thickness, legend_label=legend_label) for c in self.flist[-1:]])
    
    def plot_psum(self, m, color='red', thickness=1, legend_label="", default_label=true):
        if default_label:
            legend_label=r"$s_{{{}}}({})$".format(latex(m),latex(self.fvar))
        return plot(self.partial_sum(m), (self.fvar,-self.L(),self.L()), legend_label=legend_label, color=color, thickness=thickness)

def generate_pc(max_step = 2, max_abs = 4, zero_int=false, vari=x):
    
    r = range(-max_abs,max_abs + 1)
    r.pop(r.index(0))
    i = randint(0,2*max_abs - 1)
    c1 = r.pop(i)
    j = randint(0,2*max_abs - 2)
    c2 = r.pop(j)
    step = randint(1,max_step)

    if randint(0,1) == 0: # 2 pieces or not
        if randint(0,1) == 0: # lin start or not
            list_base = [[(0,step),c1*vari],[(step,2*step),(c1*step)*randint(0,1)]]
        else:
            list_base = [[(0,step),c1],[(step,2*step),c2]]
    else:
        list_base = [[(0,step),c1+c2*vari]]
    
    f=pc(list_base, fsymbol=function("f",vari))
    if zero_int:
        a0 = f.extension(s=0).cos_coeff(m=0)
        f = f.changed_flist([[c[0],c[1]-a0/2] for c in f.flist])
    
    return f

