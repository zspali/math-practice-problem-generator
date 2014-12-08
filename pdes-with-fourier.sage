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
        say += r"$${u}(t=0,x)={f}$$".format(u = self.u, f = latex(self.f))
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
        say += r"$$s_{{{m}}}(t,x)=\sum_{{n=1}}^{{{m}}} b_n{Tn}{Xn}={sm}$$".format(m = latex(m), Tn = latex(self.Tn()), Xn = latex(self.Xn()), sm = latex(self.psum(m)))
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
        say += r"$$s_{{{m}}}(t,x)=v(x)+\sum_{{n=1}}^{{{m}}} b_n{Tn}{Xn}={sm}$$".format(m = latex(m), Tn = latex(self.Tn()), Xn = latex(self.Xn()), sm = latex(self.psum(m)))
        return say
    
class IHE1d(HHE1d):
    
    def Xn(self):
        return cos(n*pi*x/self.f.L())
    
    def say_eqs(self):
        LHS = latex(self.alpha2 * function(self.u + "_{xx}",t,x))
        say = r"$${LHS}={u}_{{t}}(t,x)$$ ".format(u=self.u,LHS=LHS)
        say += r"$${u}_x(t,x=0)={T0},\,{u}_x(t,x={L})={T1}$$ ".format(u = self.u, T0=latex(self.bc[0]), T1=latex(self.bc[1]), L=latex(self.f.L()))
        say += r"$${u}(t=0,x)={f}$$".format(u = self.u, f = latex(self.f))
        return say
    
    def say_fseries(self):
        say = r"We can write "
        say += r"$${u}(t,x)=\frac12a_0+\sum_{{n=1}}^\infty a_n{Tn}{Xn},$$ ".format(u = self.u, Tn = latex(self.Tn()), Xn = latex(self.Xn()))
        say += r"where the $a_n$ are the Fourier sine coefficients of ${f}$.<br>".format(f=latex(self.f))
        return say + self.f.say_cos_coeff()
    
    def psum(self,m):
        an = SR(self.f.cos_coeff())
        return self.f.cos_coeff(m=0)/2+sum([an(n=k)*self.Tn()(n=k)*self.Xn()(n=k) for k in range(1, m+1)])
    
    def say_psum(self, m):
        say = r"We have "
        say += r"$$s_{{{m}}}(t,x)=\frac12a_0+\sum_{{n=1}}^{{{m}}} a_n{Tn}{Xn}={sm}$$".format(m = latex(m), Tn = latex(self.Tn()), Xn = latex(self.Xn()), sm = latex(self.psum(m)))
        return say