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

def generate_evals(max_abs = 8, evtype = 0):
    if evtype == 0:
        evtype = randint(1,3)

    if evtype == 1:
        evlist = range(-max_abs,max_abs+1)
        return [evlist.pop(choice(range(len(evlist)))),evlist.pop(choice(range(len(evlist))))]

    if evtype == 2:
        ip = randint(1,max_abs)
        mabs = floor(sqrt(max_abs^2-ip^2))
        rp = randint(-mabs,mabs)
        return [rp+(-1)^i*I*ip for i in range(2)]

    if evtype == 3:
        ev = randint(-max_abs,max_abs)
        return [ev,ev]

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
        
    return matrix(A)
    
var('t,x,y')

class HDE2:
    def __init__(self, A, t=t):
        self.t = t
        self.A = A
        
    def evals(self):
        p = self.A.trace()
        q = self.A.det()
        return [p/2+sqrt(self.Delta())/2,p/2-sqrt(self.Delta())/2]
    
    def nontrivsoln(self,B, simplify = True):
        r = filter(lambda v: v != vector([0,0]),B.rows())[0]
        if 0 != r[1]:
            ev = vector([-r[1],r[0]])
        else:
            ev = vector([r[1],-r[0]])
        if ( ev[0].abs() > ev[1].abs() and ev[1] < 0 ) or ( ev[1].abs() > ev[0].abs() and ev[0] < 0 ) or max(ev) <= 0:
            ev = -ev
        if simplify:
            m = LCM([QQ(a).denom() for a in [a for b in [[a.real(), a.imag()] for a in ev] for a in b]])
            ev = vector(map(lambda a: a/m, ev))
            m = GCD([QQ(a) for a in [a for b in [[a.real(), a.imag()] for a in ev] for a in b]])
            ev = vector(map(lambda a: a/m, ev))
        return ev
    
    def Delta(self):
        return self.A.trace()^2-4*self.A.determinant()
    
    def is_defective(self):
        return ( self.A - self.evals()[0] ) != zero_matrix(2,2)
    
    def T(self):
        if self.Delta() > 0:
            T = column_matrix([self.nontrivsoln(self.A-r) for r in self.evals()])
        elif self.Delta() < 0:
            u = self.nontrivsoln(self.A-self.evals()[0])
            T = column_matrix([map(lambda a: QQ(a.real()), u),map(lambda a: QQ(a.imag()), u)])
        elif self.Delta() == 0:
            B = (self.A - self.evals()[0])
            if self.is_defective():
                nzrow = filter(lambda r: r != vector([0,0]),B.rows())
                u = self.nontrivsoln(B, simplify = False)
                if u[0] != 0:
                    v = vector([0, u[0]/B[0][1]])
                else:
                    v = vector([u[1]/B[0][1], 0])
            else:
                u = vector([1,0])
                v = vector([0,1])
            T = column_matrix([u,v])
            m = GCD(QQ(a) for a in T.list())
            T = T/m
        return T
    def Psi(self):
        if self.Delta() > 0 or ( self.Delta() == 0 and not self.is_defective() ):
            Psi = column_matrix([exp(self.evals()[i]*self.t)*self.T().columns()[i] for i in range(2)])
        elif self.Delta() < 0:
            la = QQ(self.evals()[0].real())
            mu = QQ(self.evals()[0].imag())
            a = self.T().columns()[0]
            b = self.T().columns()[1]
            Psi = column_matrix([exp(la*self.t)*(a*cos(mu*self.t)-b*sin(mu*self.t)),exp(la*self.t)*(b*cos(mu*self.t)+a*sin(mu*self.t))])
        else:
            r = self.evals()[0]
            u = self.T().columns()[0]
            v = self.T().columns()[1]
            Psi = column_matrix([exp(r*self.t)*u,exp(r*self.t)*(v+self.t*u)])
        return Psi
    def c(self,x0):
        return self.T() \ x0
    def IC(self,x0):
        return self.Psi() * self.c(x0)
    def minmax(self,cv):
        if len(cv) == 1:
            return cv[0]-5, cv[0]+5
        elif len(cv) > 1:
            pad = (max(cv)-min(cv))/2
            return min(cv)-pad, max(cv)+pad
        else:
            return -5, 5
    def draw_plots(self,b = vector([6,2])):
        t = self.t
        xmin, xmax = self.minmax([v[0] for v in self.T().columns() + [b,vector([0,0])]])
        ymin, ymax = self.minmax([v[1] for v in self.T().columns() + [b,vector([0,0])]])
        if max(self.evals()) > 0:
            xmin,xmax,ymin,ymax = xmin*2, xmax*2, ymin*2, ymax*2
        p = plot_vector_field(self.A*vector([x,y]), (x,xmin,xmax), (y,ymin,ymax)) # Draw direction field for x'=Ax
        p = p + parametric_plot(self.Psi().columns()[0],(t,0,10),color="red",legend_label=r"$\psi_1(t)$", thickness=2) # Draw first fundamental solution
        p = p + parametric_plot(self.Psi().columns()[1],(t,0,10),color="green",legend_label=r"$\psi_2(t)$", thickness=2) # Draw second fundamental solution
        p = p + parametric_plot(self.IC(b),(t,0,10),color="blue",legend_label=r"$x(t)$", thickness=2) # Draw solution of IVP
        show(p, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
        
    def say_evals(self):
        if self.Delta() >= 0:
            evals=",".join([latex(r) for r in self.evals()])
        else:
            evals=r"\pm".join([latex(QQ(self.evals()[0].real())),latex(QQ(self.evals()[0].imag()))]) + "i"
        return r"Compute eigenvalues: $$p={0},\,q={1},\Delta={2},\,r_{{1,2}}={3}$$".format(latex(self.A.trace()),latex(self.A.det()),latex(self.Delta()),evals)
    
    def say_number(self,a):
        if a.imag() != 0 or a < 0:
            return latex(a).join(["(",")"])
        else:
            return latex(a)
        
    def say_vector(self,v):
        return latex(column_matrix([v]))
    
    def say_u(self, i=0):
        return self.say_vector(self.T().column(i))
        
    def say_T(self):
        if self.Delta() > 0:
            say = r"Contemplating the matrix $$A-{0}I={1},$$ we see that we can choose a ${0}$-eigenvector $\mathbf u_1={2}$, and contemplating the matrix $$A-{3}I={4}$$, we see that we can choose a ${3}$-eigenvector $\mathbf u_2={5}$.".format(self.say_number(self.evals()[0]),latex(self.A-self.evals()[0]),self.say_u(0),self.say_number(self.evals()[1]),latex(self.A-self.evals()[1]),self.say_u(1))
        elif self.Delta() < 0:
            say = r"Contemplating the matrix $$A-{0}I={1},$$ we see that we can choose a ${0}$-eigenvector $$\mathbf u_1={2}={3}+{4}=\mathbf a+i\mathbf b.$$".format(self.say_number(self.evals()[0]),latex(self.A-self.evals()[0]),self.say_vector(self.T().column(0)+I*self.T().column(1)),self.say_u(0),self.say_u(1))
        elif self.is_defective():
            say = r"Contemplating the matrix $$A-{0}I={1},$$ we see that we can choose a ${0}$-eigenvector $\mathbf u_1={2}$. Since $A$ is defective, contemplating the augmented matrix $$(A-{0}I|\mathbf u_1)={3},$$ we see that we can choose a generalized ${0}$-eigenvector $\mathbf v_1={4}$.".format(self.say_number(self.evals()[0]),latex(self.A-self.evals()[0]),self.say_u(0),latex((self.A-self.evals()[0]).augment(self.T().column(0),subdivide=True)),self.say_u(1))
        else:
            say = r"Since $A-{0}I=0$, we can choose the independent ${0}$-eigenvectors $\mathbf u_1={1},\,\mathbf u_2={2}$.".format(self.say_number(self.evals()[0]),self.say_u(0),self.say_u(1))
        return say
    
    def say_Psi(self):
        say = r"Then the following constitute a set of fundamental solutions."
        if self.Delta() > 0 or ( self.Delta() == 0 and not self.is_defective() ):
            say += r",$$ $$".join([r"\boldsymbol\psi_{0}({1})=e^{{{2}{1}}}\mathbf u_{0}={3}".format(i+1,latex(self.t),latex(self.evals()[i]),self.say_vector(self.Psi().column(i))) for i in range(2)]).join([r"$$",r".$$"])
        elif self.Delta() < 0:
            say += r"$$\boldsymbol\psi_1({0})=e^{{\lambda {0}}}(\mathbf a\cos(\mu {0})-\mathbf b\sin(\mu {0}))={1},$$ $$\boldsymbol\psi_2({0})=e^{{\lambda {0}}}(\mathbf b\cos(\mu {0})+\mathbf a\sin(\mu {0}))={2}.$$".format(latex(self.t),self.say_vector(self.Psi().column(0)),self.say_vector(self.Psi().column(1)))
        else:
            say += r"$$\boldsymbol\psi_1({0})=e^{{{1}{0}}}\mathbf u_1={2},$$ $$\boldsymbol\psi_2({0})=e^{{{1}{0}}}(\mathbf v_1+{0}\mathbf u_1)={3}.$$".format(latex(self.t),latex(self.evals()[0]),self.say_vector(self.Psi().column(0)),self.say_vector(self.Psi().column(1)))
        say += r"Therefore the general solution is $$x({0})=c_1\boldsymbol\psi_1({0})+c_2\boldsymbol\psi_2({0}).$$".format(latex(self.t))
        return say
