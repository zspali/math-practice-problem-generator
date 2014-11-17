R.<a>=QQ[]
var('alpha,t,u,x,y,x1,x2,c1,c2,y1,y2,s');x1f=function('x1f',t); x2f=function('x2f',t); cv = vector([c1,c2]); yx = vector([y1,y2])
  
def generate_trafo(max_abs = 8):
    T = [[1,0],[0,1]]
    d = randint(0,1)
    r = randint(0,1)
    if d == 0:
        T[r][r] = T[r][r] * randint(1,max_abs)*(-1)^randint(0,1)
    if d == 1:
        T[r][1-r] = randint(-max_abs,max_abs)
    return matrix(T)

def randomly_conjugate(A, max_abs = 8):
    T = generate_trafo(max_abs)
    A = T^(-1)*A*T
    return A

def OK_mx(A, max_abs = 8, max_denom = 2):
    return max(map(abs,A.list())) <= max_abs and lcm(map(denominator,A.list())) <= max_denom

def generate_mx(max_abs = 8, evtype = 0, parametric = False, max_conj = 50):
    if evtype == 0:
        evtype = randint(1,3)

    if evtype == 1:
        evlist = range(-max_abs,max_abs+1)
        J = matrix([[evlist.pop(choice(range(len(evlist)))),0],[0,evlist.pop(choice(range(len(evlist))))]])

    if evtype == 2:
        ip = randint(1,max_abs)*(-1)^randint(0,1)
        mabs = floor(sqrt(max_abs^2-ip^2))
        rp = randint(-mabs,mabs)
        J = matrix([[0,1],[-ip^2-rp^2,2*rp]])

    if evtype == 3:
        ev = randint(-max_abs,max_abs)
        if randint(0,5) == 0:
            J = matrix([[ev,0],[0,ev]])
        else:
            J = matrix([[ev,1],[0,ev]])

    A = J
    for i in range(max_conj):
        Ao = A
        A = randomly_conjugate(A, max_abs = max_abs)
        if not OK_mx(A, max_abs = 2*max_abs):
            A = Ao
        elif not 0 in A.list():
            break

    if parametric:
        mlist = [list(r) for r in list(A)]
        i, j, ptype = randint(0,1), randint(0,1), randint(0,1)
        if ptype == 0:
            mlist[i][j] = a
        else:
            mlist[i][j], mlist[1-i][1-j] = a, a
        A = matrix(mlist)
        
    return matrix(A)
    
def sol_row(B):
        if B[0] != vector([0,0]):
            row = B[0];rowin=0
        else:
            row = B[1];rowin=1

        if real(row[1]) >= 0 and real(row[0]) <= 0:
            evect = vector([row[1], -row[0]])
        else:
            evect = vector([-row[1], row[0]])
        return evect, rowin
    
def simplify_v(v, expr = None):
    if expr != None:
        return vector([SR(c).simplify().subs_expr(expr).expand() for c in v])
    else:
        return vector([SR(c).simplify().expand() for c in v])

def simplify_mx(A, expr = None):
    if expr != None:
        return matrix([[SR(c).simplify().subs_expr(expr).expand() for c in r] for r in A])
    else:
        return matrix([[SR(c).simplify().expand() for c in r] for r in A])
    
def latex_vmatrix(A):
    return r"\begin{{vmatrix}} {} \end{{vmatrix}}".format(reduce(lambda x,y: reduce(lambda z,w: latex(z) + " & " + latex(w), x) + r" \\ " + reduce(lambda z,w: latex(z) + " & " + latex(w), y), list(A)))
    
class HDE2d:
    def __init__(self, A, cx_exp=False, vari=t):
        self.A = A
        self.cx_exp = cx_exp
        self.t = vari

    def p(self):
        return (self.A).trace()

    def q(self):
        return (self.A).det()

    def disc(self):
        return (self.p())^2-4*self.q()

    def evals(self):
        if self.disc() == 0:
            return [self.p()/2,self.p()/2]
        else:
            if self.disc() > 0:
                return [self.p()/2 + (-1)^i*sqrt(self.disc())/2 for i in range(1,3)]
            else:
                return [self.p()/2 + (-1)^i*I*sqrt(-self.disc())/2 for i in range(0,2)]
    
    def la(self):
        return self.evals()[0].real()
    
    def mu(self):
        return self.evals()[0].imag()
    
    def ivects(self, cx_exp = False):
        if self.disc() > 0:
            return [sol_row(self.A-eval)[0] for eval in self.evals()]
        elif self.disc() < 0:
            evect = sol_row(self.A-self.evals()[0])[0]
            if cx_exp:
                return [evect, conjugate(evect)]
            else:
                return [vector([z.real() for z in evect]),vector([z.imag() for z in evect])]
        else:
            if self.A-self.p()/2 == 0:
                return [vector([1,0]),vector([0,1])]
            else:
                B=self.A-self.evals()[0]
                evect,rowin = sol_row(B)
                gvect_l=[0,0]
                gvect_l[rowin-1]=B[rowin][rowin-1]/evect[rowin]
                gvect_l[rowin]=0
                gvect=vector(gvect_l)
                return [evect,gvect]
    
    def T(self):
        return column_matrix(self.ivects())
    
    def is_defective(self):
        return bool(self.disc() == 0 and self.A-self.p()/2 != 0)
    
    def J(self):
        return diagonal_matrix(self.evals())+SR(self.is_defective())*matrix([[0,1],[0,0]])
    
    def fsols(self):
        if self.disc() > 0:
            return [e^(self.evals()[i]*self.t)*self.ivects()[i] for i in range(2)]
        elif self.disc() < 0:
            if self.cx_exp:
                return [e^(self.evals()[i]*self.t)*self.ivects(cx_exp=self.cx_exp)[i] for i in range(2)]
            else:
                la = self.evals()[0].real()
                mu = self.evals()[0].imag()
                a = self.ivects()[0]
                b = self.ivects()[1]
                return [e^(la*self.t)*(cos(mu*self.t)*a - sin(mu*self.t)*b), e^(la*self.t)*(cos(mu*self.t)*b + sin(mu*self.t)*a)]
        else:
            if self.A-self.p()/2 == 0:
                return [e^(self.evals()[0]*self.t)*self.ivects()[i] for i in range(2)]
            else:
                r = self.evals()[0]
                u = self.ivects()[0]
                v = self.ivects()[1]
                return [e^(r*self.t)*u,e^(r*self.t)*(v+self.t*u)]
            
    def Phi(self):
        return column_matrix(self.fsols())
    
    def PhiI(self):
        if self.disc() < 0:
            return simplify_mx(exp(-2*self.evals()[0].real()*self.t)/self.T().det()*self.Phi().adjoint())
        else:
            return simplify_mx(self.Phi().I())
            
    def IC(self, b):
        c = vector(self.T() \ b)
        return self.Phi()*c

    def say_evals(self):
        return r"Since $p=" + latex(self.p()) + ",\,q=" + latex(self.q()) + ",\,\Delta=" + latex(self.disc()) + r"$, the eigenvalues are $r_{1,2}=\frac{p\pm\sqrt\Delta}{2}=" + reduce(lambda x,y: latex(x) + ",\," + latex(y), self.evals()) + "$ <br>"
    
    def say_ivect(self, i=0, cx_exp=False, ui=True, usym=r"\mathbf u"):
        if ui:
            us = r"\mathbf u^{(" + latex(i+1) + r")}"
        else:
            us = usym
        return r"Since $A-(" + latex(self.evals()[i]) +")I=" + latex(self.A - self.evals()[i]) + r"$, we can choose ${}=".format(us) + latex(column_matrix([self.ivects(cx_exp=cx_exp)[i]])) + r"$"
    
    def say_evects(self):
        if self.is_defective():
            aug_mx = column_matrix((self.A - self.evals()[0]).columns() + [self.ivects()[0]])
            aug_mx_latex = r"\left(\begin{array}{cc|c}" + reduce(lambda x,y: reduce(lambda z,w: latex(z) + " & " + latex(w), x) + r" \\ " + reduce(lambda z,w: latex(z) + " & " + latex(w), y), aug_mx) + r"\end{array}\right)"
            return self.say_ivect(ui=False) + r"<br> and as $\left(A-\left({2}\right)I\mid\mathbf u\right)={0}$, we can choose $\mathbf v={1}$<br>".format(aug_mx_latex, latex(column_matrix([self.ivects()[1]])),latex(self.evals()[0]))
        elif self.disc() < 0 and not self.cx_exp:
            return self.say_ivect(i=0,cx_exp=True, ui=False) + r" $=\mathbf a+i\mathbf b=" + reduce(lambda x,y: latex(column_matrix([x]))+"+i"+latex(column_matrix([y])), self.ivects()) + r"$ <br>"
        else:
            return r"{}, <br> {} <br>".format(self.say_ivect(i=0,cx_exp=True),self.say_ivect(i=1,cx_exp=True))
        
    def evsymbols(self):
        if self.is_defective():
            return [r"\mathbf u",r"\mathbf v"]
        elif self.disc() < 0 and not self.cx_exp:
            return [r"\mathbf a", r"\mathbf b"]
        else:
            return [r"\mathbf u^{(1)}", r"\mathbf u^{(2)}"]
        
    def say_T(self):
        return self.say_evects() + r"so we have $T=\left({}\quad {}\right)={}$ <br>".format(self.evsymbols()[0],self.evsymbols()[1],latex(self.T()))
    
    def say_J(self):
        adj = [r"$\Delta\ne0$",r"$A$ is defective",r"$A$ is a scalar matrix"][(SR(bool(self.disc() == 0)))*(2-SR(self.is_defective()))]
        return r"since {}, we have $J={}$".format(adj,latex(self.J()))
    
    def say_Phi(self):
        return self.say_fsols() + r"so we have $\Phi({0})=\left(\mathbf x^{{(1)}}({0})\quad\mathbf x^{{(2)}}({0})\right)={1}$ <br>".format(latex(self.t),latex(self.Phi()))
    
    def say_det_Phi(self):
        if self.disc() < 0 and not self.cx_exp:
            et=[latex(e^(2*self.evals()[0].real()*self.t)),""][SR(bool(self.evals()[0].real() == 0))]
            
            return r"We have $\det\Phi({t})=e^{{2\lambda {t}}}\det T={et}{det_T}={det_Phi}$".format(t=latex(self.t), et=et, det_T=latex_vmatrix(self.T()), det_Phi=latex(self.Phi().det().full_simplify()))
        else:
            return r"We have $\det\Phi({t})={vPhi}={detPhi}$".format(t=latex(self.t), vPhi=latex_vmatrix(self.Phi()), detPhi=latex(self.Phi().det().full_simplify()))
    
    def say_PhiI(self):
        return self.say_det_Phi() + r"<br> so $\Phi({t})^{{-1}}=\frac{{1}}{{\det\Phi({t})}}{Phia}={PhiI}$ <br>".format(t=latex(self.t), Phia=latex(simplify_mx(self.Phi().adjoint())), PhiI=latex(self.PhiI()))
    
    def say_esol(self, i=0):
        return self.say_ivect(i,cx_exp=self.cx_exp) + r", and thus $\mathbf x^{(" + latex(i+1) + r")}(" + latex(self.t) + ")=e^{r_" + latex(i+1) + latex(self.t) + r"}\mathbf u^{(" + latex(i+1) + r")}=" + latex(column_matrix([self.fsols()[i]])) + r"$. <br>"
    
    def say_fsols(self):
        if self.disc() > 0:
            return reduce(lambda x,y: x + y, [self.say_esol(i) for i in range(2)])
        elif self.disc() < 0:
            if self.cx_exp:
                return reduce(lambda x,y: x + y, [self.say_esol(i) for i in range(2)])
            else:
                say = self.say_ivect(0,cx_exp=True) + r" $=\mathbf a+i\mathbf b=" + reduce(lambda x,y: latex(column_matrix([x]))+"+i"+latex(column_matrix([y])), self.ivects()) + r"$, so that <br>"
                say += r"$\mathbf x^{{(1)}}({0})=e^{{\lambda {0}}}\left(\cos(\mu {0})\mathbf a-\sin(\mu {0})\mathbf b\right)= {1}$ and <br>".format(latex(self.t),latex(column_matrix([self.fsols()[0]])))
                say += r"$\mathbf x^{{(2)}}({0})=e^{{\lambda {0}}}\left(\cos(\mu {0})\mathbf b+\sin(\mu {0})\mathbf a\right)={1}$. <br>".format(latex(self.t), latex(column_matrix([self.fsols()[1]])))
                return say
        else:
            if self.A-self.p()/2 == 0:
                return r"Since $A-({0})I=0$, we can choose $\mathbf x^{{(1)}}({1})=e^{{r_1{1}}}{2}={4}$ and $\mathbf x^{{(2)}}({1})=e^{{r_1{1}}}{3}={5}$ <br>".format(latex(self.evals()[0]),latex(self.t),latex(column_matrix([self.ivects()[0]])),latex(column_matrix([self.ivects()[1]])),latex(column_matrix([self.fsols()[0]])),latex(column_matrix([self.fsols()[1]])))
            else:
                aug_mx = column_matrix((self.A - self.evals()[0]).columns() + [self.ivects()[0]])
                aug_mx_latex = r"\left(\begin{array}{cc|c}" + reduce(lambda x,y: reduce(lambda z,w: latex(z) + " & " + latex(w), x) + r" \\ " + reduce(lambda z,w: latex(z) + " & " + latex(w), y), aug_mx) + r"\end{array}\right)"
                say = self.say_esol(0) + ", <br>"
                say += r"and as $\left(A-\left({2}\right)I\mid\mathbf u^{{(1)}}\right)={0}$, we can choose $\mathbf v={1}$, thus <br>".format(aug_mx_latex, latex(column_matrix([self.ivects()[1]])),latex(self.evals()[0]))
                say += r"$\mathbf x^{{(2)}}({0})=e^{{r_1 {0}}}\mathbf v+{0} e^{{r_1 {0}}}\mathbf u^{{(1)}}={1}$ <br>".format(latex(self.t), latex(column_matrix(self.fsols()[1])))
                return say
    
    def say_IC(self, b):
        c = vector(self.T() \ b)
        return r"To solve $\mathbf x({0}=0)=c_1\mathbf x^{{(1)}}({0}=0)+c_2\mathbf x^{{(2)}}({0}=0)=c_1{1}+c_2{2}={3}$, we need $c_1={4},\,c_2={5}$. <br>".format(latex(self.t),latex(column_matrix([self.ivects()[0]])),latex(column_matrix([self.ivects()[1]])),latex(column_matrix([b])),latex(c[0]),latex(c[1]))
    
    def draw_plots(self,b):
        t = self.t
        xmin, xmax = minmax([v[0] for v in self.ivects() + [b,vector([0,0])]])
        ymin, ymax = minmax([v[1] for v in self.ivects() + [b,vector([0,0])]])
        if max(self.evals()) > 0:
            xmin,xmax,ymin,ymax = xmin*2, xmax*2, ymin*2, ymax*2
        p = plot_vector_field(self.A*vector([x,y]), (x,xmin,xmax), (y,ymin,ymax)) # Draw direction field for x'=Ax
        p = p + parametric_plot(self.fsols()[0],(t,0,10),color="red",legend_label=r"$x^{(1)}(t)$", thickness=2) # Draw first fundamental solution
        p = p + parametric_plot(self.fsols()[1],(t,0,10),color="green",legend_label=r"$x^{(2)}(t)$", thickness=2) # Draw second fundamental solution
        p = p + parametric_plot(self.IC(b),(t,0,10),color="blue",legend_label=r"$x(t)$", thickness=2) # Draw solution of IVP
        show(p, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

    def say_system(self,b):
        return r"$$\begin{{aligned}} x_1'(t)&={0},\quad x_1(0)={1},\\ x_2'(t)&={2},\quad x_2(0)={3} \end{{aligned}}$$".format(latex((self.A*vector([x1,x2]))[0]),latex(b[0]),latex((self.A*vector([x1,x2]))[1]),latex(b[1]))
    
    def generate_g(self):
        c = randint(1,4)*(-1)^(randint(0,1))
        p = randint(1,4)*(-1)^(randint(0,1))
        t = self.t
        fl = [1,exp(p*t)]
        if self.disc() < 0:
            fl += [cos(self.mu()*t),sin(self.mu()*t)]
        i = randint(0,1)
        v=[0,0]
        v[i] = c*choice(fl)
        v[1-i] = 0
        return vector(v)


def poly_to_opens(p):
    p=p-a+a
    lc=p.leading_coefficient()
    lcs=lc.sign()
    if p.degree() < 1:
        return[[-infinity,infinity,lcs]]
    if p.degree() == 1:
        return [[-infinity,p.any_root(),-lcs],[p.any_root(),infinity,lcs]]
    if p.degree() == 2:
        if p.discriminant() > 0:
            return [[-infinity,-p.coeffs()[1]/2-sqrt(p.discriminant())/2,lcs],[-p.coeffs()[1]/2-sqrt(p.discriminant())/2,-p.coeffs()[1]/2+sqrt(p.discriminant())/2,-lcs],[-p.coeffs()[1]/2+sqrt(p.discriminant())/2,infinity,lcs]]
        if p.discriminant() == 0:
            return [[-infinity,-p.coeffs()[1]/2,lcs],[-p.coeffs()[1]/2,infinity,lcs]]
        if p.discriminant() < 0:
            return [[-infinity,infinity,lcs]]

def to_arrows(f,fname,fvar,y, armin = -5, armax=5, fontsize=20):
    arlen = armax - armin
    tickd = arlen / 100
    textvd = arlen * ( 3 / 100 )
    texthd = arlen / 50
    pints = arrow((armin,y),(armax,y))
    for oi in poly_to_opens(f):
        if oi[1] != infinity:
            pints += line([(oi[1],y-tickd),(oi[1],y+tickd)])
            pints += text(r"$"+latex(oi[1])+r"$",(oi[1],y+textvd),fontsize=fontsize)
            ub = oi[1]
        else:
            ub = armax
        if oi[0] != -infinity:
            lb = oi[0]
        else:
            lb = armin
        ineq = r"$" + fname
        if oi[2] == 1:
            ineq += r">0$"
        elif oi[2] == -1:
            ineq += r"<0$"
        elif oi[2] == 0:
            ineq += r"=0$"
        pints += text(ineq,((lb+ub)/2,y+textvd),fontsize=fontsize)
        pints += text("$" + fname + "=" + latex(f.substitute(a=fvar).expand()) + "$",(armax+texthd,y),fontsize=fontsize,horizontal_alignment="left")
    return pints
    
def int_cap(i1,i2):
    if i1[0] == -infinity:
        lb = i2[0]
    elif i2[0] == -infinity:
        lb = i1[0]
    else:
        lb = max ( i1[0], i2[0] )
    if i1[1] == infinity:
        ub = i2[1]
    elif i2[1] == infinity:
        ub = i1[1]
    else:
        ub = min ( i1[1], i2[1] )
    return (lb == -infinity or ub == infinity or lb < ub), lb, ub    
    
def minmax(cv):
        if len(cv) == 1:
            return cv[0]-5, cv[0]+5
        elif len(cv) > 1:
            pad = (max(cv)-min(cv))/2
            return min(cv)-pad, max(cv)+pad
        else:
            return -5, 5
    
class pHDE2d:
    def __init__(self, A):
        self.A = A

    def p(self):
        return (self.A).trace()
    
    def pints(self):
        return poly_to_opens(self.p())

    def q(self):
        return (self.A).det()
    
    def qints(self):
        return poly_to_opens(self.q())

    def disc(self):
        return (self.p())^2-4*self.q()
    
    def dints(self):
        return poly_to_opens(self.disc())
    
    def is_scalar(self):
        return self.A == self.p()/2
    
    def open_ints(self): # [q,p,disc]
        l = []
        if self.q()*self.disc() == 0:
            return [[-infinity,infinity,self.q(),self.p(),self.disc()]]
        for qi in poly_to_opens(self.q()):
            if qi[2]==-1:
                l.append(qi[:2]+[-1])
            else:
                lb = qi[0]
                ub = qi[1]
                for pi in poly_to_opens(self.p()):
                    nonempty, lb, ub = int_cap(qi,pi)
                    if nonempty:
                        for di in poly_to_opens(self.disc()):
                            nonempty2, lb2, ub2 = int_cap([lb,ub],di)
                            if nonempty2:
                                l.append([lb2,ub2]+[1]+pi[2:]+di[2:])
        return l
    
    def crit_vals(self):
        l = []
        oi = iter(self.open_ints())
        i = oi.next()
        for j in oi:
            if i[2:] != j[2:]:
                l.append(i[1])
            i = j
        return l
    
    def cv(self):
        return list(set([self.pints()[i][1] for i in range(len(self.pints())-1)]+[self.qints()[i][1] for i in range(len(self.qints())-1)]+[self.dints()[i][1] for i in range(len(self.dints())-1)]))
    
    
    
    def draw_arrows(self):
        armin,armax = minmax(self.cv())
        arrowd = ( armax - armin ) / 10
        pints = to_arrows(self.p(),"p",alpha,2*arrowd, armin=armin, armax=armax) + to_arrows(self.q(),"q",alpha,arrowd, armin=armin, armax=armax) + to_arrows(self.disc(),"\Delta",alpha,0, armin=armin, armax=armax)
        show(pints,aspect_ratio=1,axes=False)
        
    def say_open_ints(self):
        if self.q() == 0:
            return r"For every $\alpha$, the critical point at the origin is non-isolated. <br>"
        elif self.disc() == 0:
            if self.is_scalar():
                return r"For every $\alpha$, the critical point at the origin is a star point. <br>"
            else:
                return r"For every $\alpha$, the critical point at the origin is an improper node, <br>"
        else:
            say = r"The critical point at the origin is: <br>"
            for i in self.open_ints():
                if i[2] < 0:
                    say += r"a saddle point if ${0}<\alpha<{1}$ <br>".format(latex(i[0]),latex(i[1]))
                elif i[4] > 0:
                    if i[3] < 0:
                        say += r"a stable node if ${0}<\alpha<{1}$ <br>".format(latex(i[0]),latex(i[1]))
                    elif i[3] > 0:
                        say += r"an unstable node if ${0}<\alpha<{1}$ <br>".format(latex(i[0]),latex(i[1]))
                elif i[4] < 0:
                    if i[3] < 0:
                        say += r"a stable spiral point if ${0}<\alpha<{1}$ <br>".format(latex(i[0]),latex(i[1]))
                    elif i[3] > 0:
                        say += r"an unstable spiral point if ${0}<\alpha<{1}$ <br>".format(latex(i[0]),latex(i[1]))
                    else:
                        say += r"a center point if ${0}<\alpha<{1}$ <br>".format(latex(i[0]),latex(i[1]))
        return say
    
    def say_crit_vals(self):
        if self.crit_vals() == []:
                    return "There are no critical values"
        say = "The critical point at the origin at the critical value <br>"
        for v in self.crit_vals():
            if self.q().substitute(a=v) == 0:
                say += r"$\alpha={0}$ is non-isolated <br>".format(latex(v))
            elif self.disc().substitute(a=v) == 0:
                if matrix([[f.substitute(a=v) for f in r] for r in self.A.rows()]) == self.p()/2:
                    say += r"$\alpha={0}$ is a star point <br>".format(latex(v))
                else:
                    say += r"$\alpha={0}$ is an improper node <br>".format(latex(v))
            else:
                say += r"$\alpha={0}$ is a center point <br>".format(latex(v))
        return say
    
    def draw_plots(self):
        armin,armax = minmax(self.cv())
        xmin,xmax,ymin,ymax = -5,5,-5,5
        @interact
        def _f(alpha=(0,(armin,armax)),x01=(2,(xmin,xmax)),x02=(1,(ymin,ymax))):
            A = matrix([[f.substitute(a=alpha) for f in r] for r in self.A.rows()])
            vec = vector([x1f,x2f])
            des = [diff(x1f,t) == (A*vec)[0], diff(x2f,t) == (A*vec)[1]]
            sol = map(lambda eq: eq.right_hand_side(), desolve_system(des, [x1f,x2f], ics = [0, 1, -2]))
            g = plot_vector_field(A*vector([x1,x2]),(x1,xmin,xmax),(x2,ymin,ymax))
            g += parametric_plot(sol, (t,0,10), thickness=2, legend_label=r"$x(t)$")
            show(g, xmin=xmin, xmax =xmax, ymin = ymin, ymax = ymax)



class IHDE2d(HDE2d):
    def __init__(self, A, g, cx_exp=False, vari=t, t0=0):
        self.A = A
        self.cx_exp = cx_exp
        self.t = vari
        self.g = g
        assume(vari > 0)
        self.t0=t0
    
    def du(self):
        return vector([SR(c) for c in self.PhiI()*self.g])
    
    def u(self):
        t = self.t
        return vector([SR(integrate(SR(c)(t=s),(s,self.t0,t))) for c in (self.PhiI()*self.g)])
    
    def h(self):
        return self.T().I*self.g
    
    def y(self):
        t = self.t
        if self.disc() != 0:
            return vector([exp(self.evals()[k]*t)*SR(integrate(exp(-self.evals()[k]*s)*SR(self.h()[k])(t=s), (s,self.t0,t))).full_simplify() for k in range(2)])
        else:
            y2 = exp(self.evals()[1]*t)*SR(integrate(exp(-self.evals()[1]*s)*SR(self.h()[1])(t=s), (s,self.t0,t))).full_simplify()
            y1 = exp(self.evals()[0]*t)*SR(integrate(exp(-self.evals()[0]*s)*SR(self.h()[0]+y2)(t=s), (s,self.t0,t))).full_simplify()
            return vector([y1,y2])
        
    def say_h(self):
        return r"$\mathbf h({0})=T^{{-1}}\mathbf g({0})={1}{2}={3}$".format(latex(self.t),latex(self.T().I),latex(self.g.column()),latex(self.h().column()))
    
    def say_yi(self, i=0):
        hj = [self.h()[i],self.h()[i]+self.y()[1-i]][SR(bool(i == 0 and self.is_defective()))]
        return r"$$y_{{{j}}}({t})=e^{{r_{{{j}}}{t}}}\int_{{{t0}}}^{{{t}}}e^{{-r_{{{j}}}s}}\left({hj}(s)\right)\,\mathrm ds={yj}$$".format(j=latex(i+1), t=latex(self.t), t0=latex(self.t0), hj=latex(hj), yj=latex(self.y()[i]))
    
    def say_y(self):
        say = self.say_evals() + self.say_T() + r"Letting $\mathbf x({0})=T\mathbf y({0})$, the system gets transformed to $\mathbf y'({0})=J\mathbf y({0})+\mathbf h({0})$, where $J$ is the Jordan form: <br> {1} and {2} <br> Now we can get a particular solution with transformed coordinates".format(latex(self.t),self.say_J(),self.say_h())
        i = SR(self.is_defective())
        return say + self.say_yi(i=i) + self.say_yi(i=1-i)
    
    def say_u(self):
        t = self.t
        du = vector([SR(c)(t=s) for c in self.du()])
        say = r"A particular solution will be $\mathbf u({t})=\int_{{{t0}}}^{{{t}}}\Phi(s)^{{-1}}\mathbf g(s)$ <br>".format(t=latex(self.t), t0=latex(self.t0))
        say += self.say_Phi() + self.say_PhiI()
        say += r"Therefore $\Phi(s)^{{-1}}\mathbf g(s)={du}$, and so <br>".format(du=latex(du.column()))
        say += reduce(lambda x,y: x + y, [r"$$u_{{{i}}}({t})=\int_{{{t0}}}^{{{t}}}{dui}\,\mathrm ds={ui}$$".format(i=latex(i+1), t=latex(self.t), t0=latex(self.t0), dui=latex(du[i]), ui=latex(self.u()[i])) for i in range(2)])
        return say