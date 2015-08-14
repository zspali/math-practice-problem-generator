var('t')

class HHE1d:
    
    def __init__(self, f, u="u", alpha2=1, bc=[0,0]):
        self.u = u
        self.f = f
        self.alpha2=alpha2
        self.bc = bc
        
    def Tn(self):
        return exp(-self.alpha2 * n^2 * pi^2 * t / self.f.L()^2)
    
    def Xn(self):
        return sin(n*pi*x/self.f.L())
    
    def say_eqs(self):
        LHS = latex(self.alpha2 * function(self.u + "_{xx}",t,x))
        say = r"$${LHS}={u}_{{t}}(t,x)$$ ".format(u=self.u,LHS=LHS)
        say += r"$${u}(t,x=0)={T0},\,{u}(t,x={L})={T1}$$ ".format(u = self.u, T0=latex(self.bc[0]), T1=latex(self.bc[1]), L=latex(self.f.L()))
        say += r"$${u}(t=0,x)={f}$$ ".format(u = self.u, f = latex(self.f))
        return say
    
    def say_fseries(self):
        say = r"We can write "
        say += r"$${u}(t,x)=\sum_{{n=1}}^\infty b_n{Tn}{Xn},$$ ".format(u = self.u, Tn = latex(self.Tn()), Xn = latex(self.Xn()))
        say += r"where the $b_n$ are the Fourier sine coefficients of ${f}$.<br>".format(f=latex(self.f))
        return say + self.f.say_sin_coeff()
    
    def psum(self,m):
        bn = SR(self.f.sin_coeff())
        return sum([bn(n=k)*self.Tn()(n=k)*self.Xn()(n=k) for k in range(1, m+1)])
    
    def say_psum(self, m):
        say = r"We have "
        say += r"$$s_{{{m}}}(t,x)=\sum_{{n=1}}^{{{m}}} b_n{Tn}{Xn}={sm}$$ ".format(m = latex(m), Tn = latex(self.Tn()), Xn = latex(self.Xn()), sm = latex(self.psum(m)))
        return say
    
    def plot_psum(self, m, tmax = 1):
        return contour_plot(self.psum(m), (x,0,self.f.L()), (t,0,tmax), contours = 50)
    
class NHHE1d(HHE1d):
    def v(self):
        return self.bc[0] + (self.bc[1] - self.bc[0])/self.f.L() * x
    
    def w(self):
        fmvl = [[c[0],c[1]-self.v()] for c in self.f.flist]
        fmv = pc(fmvl, fsymbol = function("f",x)-function("v",x))
        return HHE1d(f = fmv, u = "w", alpha2 = self.alpha2)
    
    def say_fseries(self):
        say = r"We have $$u(t,x)=v(x)+w(t,x),$$ where $$v(x)={v}$$, and $w(t,x)$ is a formal solution to ".format(v = latex(self.v()))
        return say + self.w().say_eqs() + self.w().say_fseries()
    
    def psum(self,m):
        return self.v() + self.w().psum(m)
    
    def say_psum(self, m):
        say = r"We have "
        say += r"$$s_{{{m}}}(t,x)=v(x)+\sum_{{n=1}}^{{{m}}} b_n{Tn}{Xn}={sm}$$ ".format(m = latex(m), Tn = latex(self.Tn()), Xn = latex(self.Xn()), sm = latex(self.psum(m)))
        return say
    
class IHE1d(HHE1d):
    
    def Xn(self):
        return cos(n*pi*x/self.f.L())
    
    def say_eqs(self):
        LHS = latex(self.alpha2 * function(self.u + "_{xx}",t,x))
        say = r"$${LHS}={u}_{{t}}(t,x)$$ ".format(u=self.u,LHS=LHS)
        say += r"$${u}_x(t,x=0)={T0},\,{u}_x(t,x={L})={T1}$$ ".format(u = self.u, T0=latex(self.bc[0]), T1=latex(self.bc[1]), L=latex(self.f.L()))
        say += r"$${u}(t=0,x)={f}$$ ".format(u = self.u, f = latex(self.f))
        return say
    
    def say_fseries(self):
        say = r"We can write "
        say += r"$${u}(t,x)=\frac12a_0+\sum_{{n=1}}^\infty a_n{Tn}{Xn},$$ ".format(u = self.u, Tn = latex(self.Tn()), Xn = latex(self.Xn()))
        say += r"where the $a_n$ are the Fourier cosine coefficients of ${f}$.<br>".format(f=latex(self.f))
        return say + self.f.say_cos_coeff()
    
    def psum(self,m):
        an = SR(self.f.cos_coeff())
        return self.f.cos_coeff(m=0)/2+sum([an(n=k)*self.Tn()(n=k)*self.Xn()(n=k) for k in range(1, m+1)])
    
    def say_psum(self, m):
        say = r"We have "
        say += r"$$s_{{{m}}}(t,x)=\frac12a_0+\sum_{{n=1}}^{{{m}}} a_n{Tn}{Xn}={sm}$$ ".format(m = latex(m), Tn = latex(self.Tn()), Xn = latex(self.Xn()), sm = latex(self.psum(m)))
        return say
    
class ZVWE1d(HHE1d):
    
    def __init__(self, f, u="u", a=1):
        self.u = u
        self.f = f
        self.a=a
    
    def Tn(self):
        return cos(self.a*n*pi*t/self.f.L())
    
    def Xn(self):
        return sin(n*pi*x/self.f.L())
    
    def say_eqs(self):
        LHS = latex(self.a^2 * function(self.u + "_{xx}",t,x))
        say = r"$${LHS}={u}_{{tt}}(t,x)$$ ".format(u=self.u,LHS=LHS)
        say += r"$${u}(t,x=0)=0,\,{u}(t,x={L})=0$$ ".format(u = self.u, L=latex(self.f.L()))
        say += r"$${u}(t=0,x)={f}$$ ".format(u = self.u, f = latex(self.f))
        say += r"$${u}_t(t=0,x)=0$$ ".format(u = self.u)
        return say
    
    def say_fseries(self):
        say = r"We can write "
        say += r"$${u}(t,x)=\sum_{{n=1}}^\infty b_n{Tn}{Xn},$$ ".format(u = self.u, Tn = latex(self.Tn()), Xn = latex(self.Xn()))
        say += r"where the $b_n$ are the Fourier sine coefficients of ${f}$.<br>".format(f=latex(self.f))
        return say + self.f.say_sin_coeff()
    
class ZDWE1d(HHE1d):
    
    def __init__(self, f, u="u", a=1):
        self.u = u
        self.f = f
        self.a=a
    
    def Tn(self):
        return sin(self.a*n*pi*t/self.f.L())
    
    def Xn(self):
        return sin(n*pi*x/self.f.L())
    
    def say_eqs(self):
        LHS = latex(self.a^2 * function(self.u + "_{xx}",t,x))
        say = r"$${LHS}={u}_{{tt}}(t,x)$$ ".format(u=self.u,LHS=LHS)
        say += r"$${u}(t,x=0)=0,\,{u}(t,x={L})=0$$ ".format(u = self.u, L=latex(self.f.L()))
        say += r"$${u}(t=0,x)=0$$ ".format(u = self.u)
        say += r"$${u}_t(t=0,x)={f}$$ ".format(u = self.u, f = latex(self.f))
        return say
    
    def say_fseries(self):
        say = r"We can write "
        say += r"$${u}(t,x)=\sum_{{n=1}}^\infty k_n{Tn}{Xn},$$ ".format(u = self.u, Tn = latex(self.Tn()), Xn = latex(self.Xn()))
        say += r"where the $k_n{anpipL}=b_n$ are the Fourier sine coefficients of ${f}$.<br>".format(f=latex(self.f),anpipL=latex(self.a*n*pi/self.f.L()),L=self.f.L())
        return say + self.f.say_sin_coeff()
    
    def kn(self):
        return self.f.L()/self.a/n/pi*self.f.sin_coeff()
    
    def psum(self,m):
        return sum([self.kn()(n=k)*self.Tn()(n=k)*self.Xn()(n=k) for k in range(1, m+1)])
    
    def say_psum(self, m):
        say = r"We have "
        say += r"$$s_{{{m}}}(t,x)=\sum_{{n=1}}^{{{m}}} k_n{Tn}{Xn}={sm}$$ ".format(m = latex(m), Tn = latex(self.Tn()), Xn = latex(self.Xn()), sm = latex(self.psum(m)))
        return say

var('x,y')

def insert_after_first(st, ins):
    Xnl = list(st)
    Xnl.insert(1,ins)
    return "".join(Xnl)
    
class DP2d:
    
    def __init__(self, bc, ab): # bc = [ [x=0, x=a], [y=0, y=b] ]
        self.bc = bc
        self.ab = ab
        
    def XisF(self):
        return self.bc[0] == [0,0]
    
    def x(self):
        return [y,x][self.XisF()]
    
    def a(self):
        return self.ab[1-self.XisF()]
    
    def b(self):
        return self.ab[self.XisF()]
    
    def y(self):
        return [x,y][self.XisF()]
    
    def v_is_y(self):
        return bool(self.bc[self.XisF()][0] == 0)
    
    def v(self):
        return self.y() - (1 - self.v_is_y()) * self.b()
    
    def Fn(self):
        return sin ( n * pi * self.x() / self.a() )
    
    def Kn(self):
        return sinh ( n * pi * self.v() / self. a() )
    
    def Knev(self):
        vari = self.y()
        valu = self.y_subs()
        return self.Kn().subs(vari==valu)
    
    def f(self):
        return self.bc[self.XisF()][self.v_is_y()]
    
    def y_subs(self):
        return [0,self.b()][self.v_is_y()]
    
    def kn(self):
        vari = self.y()
        valu = self.y_subs()
        return self.f().sin_coeff() / self.Knev()
                
    def psum(self, m):
        return sum([self.kn()(n=k)*self.Kn()(n=k)*self.Fn()(n=k) for k in range(1, m+1)])
    
    def plot_psum(self, m):
        return contour_plot(self.psum(m), (x,0,self.ab[0]), (y,0,self.ab[1]), contours = 50)
    
    
    def say_eqs(self):
        say = r"$$ v_{xx}(x,y)+v_{yy}(x,y)=0 $$ "
        say += r"$$ v(x=0,y)={k},\,v(x={a},y)={f} $$ ".format(k=latex(self.bc[0][0]), a = latex(self.ab[0]), f = latex(self.bc[0][1]))
        say += r"$$ v(x,y=0)={h},\,v(x,y={b})={g} $$ ".format(h=latex(self.bc[1][0]), b = latex(self.ab[1]), g = latex(self.bc[1][1]))
        return say
    
    def say_fseries(self):
        say = r"Separating the variables: $v(x,y)=X(x)Y(y)$, the PDE gives the parametric ODE's "
        pdelist = [r"X''(x)-\lambda X(x)=0",r"Y''(y)+\lambda Y(y)"]
        say += r"$$ X''(x)-\lambda X(x)=0,\,Y''(y)+\lambda Y(y)$$ "
        say += r"The homogeneous boundary conditions become "
        bclist = [ "X(x={})=0".format(latex(k)) for k in [0,self.ab[0]]] + [ "Y(y={})=0".format(latex(k)) for k in [0,self.ab[1]]]
        bclist.pop(2*self.XisF() + self.v_is_y())
        say += "$${}$$ ".format(reduce( lambda x,y: x + r",\," + y, bclist))
        X = r"{X}({x})".format(X=latex(self.x()).capitalize(), x=latex(self.x()))
        Y = r"{Y}({y})".format(Y=latex(self.y()).capitalize(), y=latex(self.y()))
        say += r"Since there are two HBC's for ${X}$, we need the solution of the eigenfunction problem for ".format(X=X)
        say += r"$$ {PODE},\,{BC1},\,{BC2} $$ ".format(PODE=pdelist[1-self.XisF()], BC1=bclist[1-self.XisF()], BC2=bclist[2-self.XisF()])
        say += r"The eigenvalues and eigenfunctions are "
        Xnl = list(X)
        Xnl.insert(1,"_n")
        Xn = reduce(lambda x,y: x+y, Xnl)
        say += r"$$ \lambda_n={lambdan},\,{Xn}={Fn},\,n=1,2,\dotsc $$ ".format(lambdan = latex((-1)^self.XisF()*n^2*pi^2/self.a()^2), Xn=Xn, Fn=latex(self.Fn()))
        Ydd = insert_after_first(Y, "''")
        say += r"Substituting $\lambda=\lambda_n$ to the other ODE, we need to find a nonzero solution to "
        say += r"$${Ydd}-{npipa}{Y}=0,\,{BC}$$ ".format(Y=Y, Ydd=Ydd, npipa=latex(n^2*pi^2/self.a()^2),BC=bclist[2 * self.XisF()])
        Yn = insert_after_first(Y, "_n")
        say += r"Which is $${Yn}={Kn}$$ ".format(Yn=Yn, Kn=latex(self.Kn()))
        say += r"If we take the formal function $$v(x,y)=\sum_{{n=1}}^\infty k_n{Kn}{Fn}$$ ".format(Kn=latex(self.Kn()), Fn=latex(self.Fn()))
        say += r"Then the nonhomogeneous boundary condition reads "
        say += r"$$ {f}=v({x},{y})=\sum_{{n=1}}^\infty k_n{Knev}{Fn} $$ ".format(x=latex(self.x()), y=latex(self.y() == self.y_subs()), Knev=latex(self.Knev()), Fn=latex(self.Fn()), f=latex(self.f()))
        say += r"That is we need the $k_n{Knev}$ to be the Fourier sine coefficients $b_n$ of ${f}$. <br>".format(Knev = latex(self.Knev()), f=latex(self.f()))
        return say + self.f().say_sin_coeff()
    
    def say_psum(self, m):
        return r"We have $$s_{{{m}}}(x,y)=\sum_{{n=1}}^{{{m}}} k_n{Kn}{Fn}={sm}$$ ".format(Kn=latex(self.Kn()), Fn=latex(self.Fn()), m=latex(m), sm=latex(self.psum(m)))
    
class NP2d(DP2d):
    
    def Fn(self):
        return cos ( n * pi * self.x() / self.a() )
    
    def Kn(self):
        return cosh ( n * pi * self.v() / self. a() )
    
    def Knev(self):
        vari = self.y()
        valu = self.y_subs()
        return diff(self.Kn(), vari).subs(vari==valu)
    
    def kn(self):
        vari = self.y()
        valu = self.y_subs()
        return self.f().cos_coeff() / self.Knev()
    
    def say_eqs(self):
        say = r"$$ v_{xx}(x,y)+v_{yy}(x,y)=0 $$ "
        say += r"$$ v_x(x=0,y)={k},\,v_x(x={a},y)={f} $$ ".format(k=latex(self.bc[0][0]), a = latex(self.ab[0]), f = latex(self.bc[0][1]))
        say += r"$$ v_y(x,y=0)={h},\,v_y(x,y={b})={g} $$ ".format(h=latex(self.bc[1][0]), b = latex(self.ab[1]), g = latex(self.bc[1][1]))
        return say
    
    def say_fseries(self):
        say = r"Separating the variables: $v(x,y)=X(x)Y(y)$, the PDE gives the parametric ODE's "
        pdelist = [r"X''(x)-\lambda X(x)=0",r"Y''(y)+\lambda Y(y)"]
        say += r"$$ X''(x)-\lambda X(x)=0,\,Y''(y)+\lambda Y(y)$$ "
        say += r"The homogeneous boundary conditions become "
        bclist = [ "X(x={})=0".format(latex(k)) for k in [0,self.ab[0]]] + [ "Y(y={})=0".format(latex(k)) for k in [0,self.ab[1]]]
        bclist.pop(2*self.XisF() + self.v_is_y())
        bcplist = map(lambda x: insert_after_first(x, "'"), bclist)
        say += r"$${}$$ ".format(reduce(lambda x,y: x + r",\," + y, bcplist))
        X = r"{X}({x})".format(X=latex(self.x()).capitalize(), x=latex(self.x()))
        Y = r"{Y}({y})".format(Y=latex(self.y()).capitalize(), y=latex(self.y()))
        say += r"Since there are two HBC's for ${X}$, we need the solution of the eigenfunction problem for ".format(X=X)
        say += r"$$ {PODE},\,{BC1},\,{BC2} $$ ".format(PODE=pdelist[1-self.XisF()], BC1=bcplist[1-self.XisF()], BC2=bcplist[2-self.XisF()])
        say += r"The eigenvalues and eigenfunctions are "
        Xnl = list(X)
        Xnl.insert(1,"_n")
        Xn = reduce(lambda x,y: x+y, Xnl)
        say += r"$$ \lambda_n={lambdan},\,{Xn}={Fn},\,n=0,1,2,\dotsc $$ ".format(lambdan = latex((-1)^self.XisF()*n^2*pi^2/self.a()^2), Xn=Xn, Fn=latex(self.Fn()))
        Ydd = insert_after_first(Y, "''")
        say += r"Substituting $\lambda=\lambda_n$ to the other ODE, we need to find a nonzero solution to "
        say += r"$${Ydd}-{npipa}{Y}=0,\,{BC}$$ ".format(Y=Y, Ydd=Ydd, npipa=latex(n^2*pi^2/self.a()^2),BC=insert_after_first(bclist[2 * self.XisF()],"'"))
        Yn = insert_after_first(Y, "_n")
        say += r"Which is $${Yn}={Kn}$$ ".format(Yn=Yn, Kn=latex(self.Kn()))
        say += r"If we take the formal function $$v(x,y)=\sum_{{n=0}}^\infty k_n{Kn}{Fn}$$ ".format(Kn=latex(self.Kn()), Fn=latex(self.Fn()))
        say += r"Then the nonhomogeneous boundary condition reads "
        say += r"$$ {f}=v_{{{dy}}}({x},{y})=\sum_{{n=1}}^\infty k_n{Knev}{Fn} $$ ".format(dy = latex(self.y()), x=latex(self.x()), y=latex(self.y() == self.y_subs()), Knev=latex(self.Knev()), Fn=latex(self.Fn()), f=latex(self.f()))
        say += r"That is we need the $k_n{Knev}$ to be the Fourier cosine coefficients $a_n$ of ${f}$ for $n>0$, and we need to check that $a_0=0$. We can let $k_0=0$. <br>".format(Knev = latex(self.Knev()), f=latex(self.f()))
        return say + self.f().say_cos_coeff()
    
