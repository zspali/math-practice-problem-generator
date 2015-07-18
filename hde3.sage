var('t,x,y')

class HDE3:
    def __init__(self, A=zero_matrix(3,3), J=zero_matrix(3,3), is_defective=0, t=t, max_abs = 8, max_denom = 2, gen_n_tries = 1000):
        self.t = t
        self.A = A
        self.J = J
        self.is_defective = is_defective
        self.max_abs = max_abs
        self.max_denom = max_denom
        self.gen_n_tries = gen_n_tries
        
        if self.A == zero_matrix(3,3):
            if self.J == zero_matrix(3,3):
                if self.is_defective == 0:
                    self.is_defective = randint(1,2)
                self.is_defective = bool(self.is_defective - 1)
                if not self.is_defective:
                    self.J = diagonal_matrix([randint(-self.max_abs,self.max_abs) for i in range(3)])
                else:
                    ints = range(-self.max_abs, self.max_abs+1)
                    i = randint(0,2*self.max_abs)
                    r1 = ints.pop(i)
                    i = randint(0,2*self.max_abs-1)
                    r2 = ints.pop(i)
                    self.J = matrix(QQ,[[r1,0,0],[0,r2,1],[0,0,r2]])
            self.A = self.generate_mx()
        self.J = self.A.jordan_form()
    def generate_trafo(self):
        T = [[1,0,0],[0,1,0],[0,0,1]]
        d = randint(0,3)
        if d == 0:
            i = randint(0,2)
            T[i][i] = T[i][i] * randint(1,self.max_abs)*(-1)^randint(0,1)
        else:
            ids = [0,1,2]
            ids.pop(randint(0,2))
            i = randint(0,1)
            T[ids[i]][ids[1-i]] = randint(-self.max_abs,self.max_abs)
        return matrix(T)

    def randomly_conjugate(self, A):
        T = self.generate_trafo()
        A = T^(-1)*A*T
        return A

    def OK_mx(self, A):
        return max(map(abs,A.list())) <= self.max_abs and lcm(map(denominator,A.list())) <= self.max_denom

    def great_mx(self, A):
        Av = A - t*matrix.identity(3)
        return max(map(lambda v:v.list().count(0), Av.rows()+Av.columns())) == 1

    def generate_mx(self):
        A = self.J
        for i in range(self.gen_n_tries):
            Ao = A
            A = self.randomly_conjugate(A)
            if not self.OK_mx(A):
                A = Ao
            if self.great_mx(A):
                break
        return matrix(QQ,A)
    
    def sing_mx(self, ev):
        return self.A - ev
    
    def two_nz_entries(self,r, simplify = True):
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
    
    def set_two_nz_entries(self, r, i, simplify = True):
        ii = range(3)
        ii.pop(i)
        [i0,i1]=ii
        rs = [r[j] for j in ii]
        v = range(3)
        vs = self.two_nz_entries(rs, simplify = simplify)
        v[i0]=vs[0]
        v[i1]=vs[1]
        v[i]=0
        return vector(v)
    
    def one_nz(self, B, ii = range(3)):
        l = filter(lambda r: len(filter(lambda x: x != 0, r[1])) == 1, [[i,B.row(i)] for i in ii])
        return [r + [filter(lambda x: x[1] != 0,[[i,r[1][i]] for i in range(3)])[0][0]] for r in l]
    
    def num_nz(self, l):
        return len(filter(lambda x: x != 0, l))
    
    def say_number(self,a):
        if a.imag() != 0 or a < 0:
            return latex(a).join(["(",")"])
        else:
            return latex(a)
        
    def say_vector(self,v):
        return latex(column_matrix([v]))
    
    def latex_evlist(self, evlist, start_index):
        return r"\,".join(["=".join([r"\mathbf u_{0}".format(i+start_index+1),self.say_vector(evlist[i])]) for i in range(len(evlist))])
    
    def zeros_diff_columns(self,B, ii = range(3)):
        return filter(lambda p: ([] not in [k[2] for k in p]) and (p[0][2][0] != p[1][2][0]) and p[0][1] != 0 and p[1][1] != 0,[[[i,B.row(i),filter(lambda k: B.row(i)[k] == 0, range(3))] for i in j] for j in [[j0,j1] for j0 in ii for j1 in filter(lambda x: x != j0, ii)]])
    
    def find_evects(self, ev, ev_mult, say, start_index):
        B = self.A-ev
        r = B.rank()
        if say:
            html(r"Contemplate the matrix $A-{0}\cdot I={1}$.<br>".format(self.say_number(ev),latex(B)))
        if r == 0:
            evlist = identity_matrix(3).columns()
            latex_evlist = r"\,".join(["=".join([r"\mathbf u_{0}".format(i),self.say_vector(evlist[i])]) for i in range(3)])
            if say:
                html(r"Since $A={0}I$, we have an independent set consisting of 3 ${0}$-eigenvectors ${1}$.<br>".format(self.say_number(ev), latex_evlist))
            return [identity_matrix(3).columns()]
        elif r == 1:
            nzcolumns = filter(lambda v: v[1] != 0, [[i, B.column(i)] for i in range(3)])
            nnz = len(nzcolumns)
            if nnz == 1:
                evlist = [identity_matrix(3).column(i) for i in filter(lambda j: j != nzcolumns[0][0], range(3))]
            elif nnz == 2:
                l = [v[0] for v in nzcolumns]
                i = filter(lambda j: j not in l, range(3))[0]
                v1 = identity_matrix(3).column(i)
                l2 = range(3)
                l2[i] = 0
                r = filter(lambda row: row != 0, B.rows())[0]
                l3 = self.two_nz_entries(filter(lambda j: j != 0, r))
                l2[l[0]]  = l3[0]
                l2[l[1]] = l3[1]
                evlist = [v1, vector(l2)]
            else:
                r = filter(lambda v: v != 0, B.rows())[0]
                zero_indices=[j[0] for j in sorted([[i,r[i]] for i in range(3)],key = lambda v: -abs(v[1]))[:2]]
                evlist=[self.set_two_nz_entries(r, zero_indices[j]) for j in range(2)]
            if say:
                html(r"It has rank 1, and we get 2 independent ${0}$-eigenvectors ${2}$.<br>".format(self.say_number(ev),latex(B), self.latex_evlist(evlist, start_index)))
            return evlist
        else:
            reduction_steps=[] # 0:R_0+1*R_2, 1:0*R_1
            return self.reduction_step(B, ev, ev_mult, [], say, start_index)
                
    def reduction_step(self, B, ev, ev_mult, reduction_steps, say, start_index):
        ii = range(3)
        zr = filter(lambda r: r[1] == 0, [[i, B.row(i)] for i in ii])
        if zr != []:
            ii.remove(zr[0][0])
        dep_l = filter(lambda p: matrix([r[1] for r in p]).rank() == 1, [[[j,B.row(j)] for j in i] for i in [[k0,k1] for k0 in ii for k1 in filter(lambda x: x != k0, ii)]])
        if dep_l != []:
            ii.remove(max(dep_l[0],key = lambda x: max(x[1],key = lambda y: abs(y)))[0])
        z_l = [[i,B.row(i),filter(lambda j: B.row(i)[j] == 0,range(3))] for i in ii]
        zdc = filter(lambda x: min([len(set(z_l[x[0]][2]+z_l[x[1]][2])), len(set(z_l[x[0]][2] + [x[2]])), len(set(z_l[x[1]][2] + [x[2]]))]) > 1,[i+[k] for i in [[j0,j1] for j0 in range(len(ii)) for j1 in filter(lambda x: x != j0, range(len(ii)))] for k in range(3)])
        
        if zdc != []:
            dc = zdc[0]
            fi = dc[2]
            el = range(3)
            el[fi] = 1
            sc = filter(lambda x: x != fi, range(3))
            sr = [max(range(2), key = lambda x: abs(z_l[dc[x]][1][j])) for j in sc]
            tzil = filter(lambda i: len(z_l[i][2]) == 2, dc[:2])
            if tzil != []:
                tzi = filter(lambda i: i not in z_l[tzil[0]][2], range(3))[0]
                oi = filter(lambda i: i != tzil[0], dc[:2])[0]
                eve = self.set_two_nz_entries(z_l[oi][1],tzi)
            else:
                el[sc[0]] = -z_l[dc[sr[0]]][1][fi]/z_l[dc[sr[0]]][1][sc[0]]
                el[sc[1]] = -z_l[dc[sr[1]]][1][fi]/z_l[dc[sr[1]]][1][sc[1]]
                eve = vector(el) * lcm([x.denominator() for x in el])
                eve *= 1/gcd(eve.list())
            
            if ev_mult > 1:
                reve = copy(eve)
                for rs in reduction_steps:
                    if rs[0] == 0:
                        reve[rs[1]] += rs[2] * reve[rs[3]]
                    else:
                        reve[rs[2]] *= rs[1]
                BA = B.augment(reve)
                i0 = max(filter(lambda i: BA.row(i)[3] != 0, ii), key = lambda i: len(filter(lambda y: y == 0, BA.row(i))))
                i1 = max(filter(lambda i: i != i0, ii), key = lambda i: len(filter(lambda y: y == 0, BA.row(i))))
                
                
                
                ri = [i0,i1]
                rtd = filter(lambda x: x not in ri, range(3))[0]
                Br = B.delete_rows([rtd])
                for rs in reduction_steps:
                    if rs[0] == 0:
                        reve[rs[1]] += rs[2] * reve[rs[3]]
                    else:
                        reve[rs[2]] *= rs[1]
                    
                ever = vector([reve[i] for i in filter(lambda i: i != rtd, range(3))])
                ctd = max(filter(lambda x: Br.delete_columns([x]).rank() == 2, range(3)), key = lambda x: len(filter(lambda i: i == 0, Br.delete_columns([x]).list() )))
                Brr = Br.delete_columns([ctd])
                geve = vector(SR,range(3))
                geve[ctd] = 0
                gever = Brr \ ever
                vi = filter(lambda i: i != ctd, range(3))
                geve[vi[0]] = gever[0]
                geve[vi[1]] = gever[1]
                m1 = lcm([x.denominator() for x in geve])
                geve *= m1
                eve *= m1
                m2 = gcd(eve.list() + geve.list())
                geve *= m2
                eve *= m2
                if say:
                    html(r"This shows that we can choose a ${0}$-eigenvector $\mathbf u_{1}={2}$.<br>".format(self.say_number(ev),latex(start_index+1),self.say_vector(eve)))
                    text = r"A generalized eigenvector $\mathbf v_{0}$ is a solution to the equation defined by the augmented matrix $(A-{1}I|\mathbf u_{0})$.".format(latex(start_index+1),self.say_number(ev))
                    if reduction_steps != []:
                        text += r" Using the same row reduction, this turns into ${0}$.".format(latex(B.augment(reve*m1*m2, subdivide=True)))
                    text += r" Contemplating this, we see that we can choose $\mathbf v_{0}={1}$.<br>".format(latex(start_index+1),self.say_vector(geve))
                    html(text)
                return [eve,geve]
                

            if say:
                html(r"This shows that we can choose a ${0}$-eigenvector $\mathbf u_{1}={2}$.<br>".format(self.say_number(ev),latex(start_index+1),self.say_vector(eve)))
            
            return [eve]
        
        p = min([[i[0],i[1],k] for i in [[j0,j1] for j0 in range(len(ii)) for j1 in filter(lambda x: x != j0, range(len(ii)))] for k in filter(lambda x: z_l[i[0]][1][x] != 0, range(3))], key = lambda x:[-len(z_l[x[0]][2]),abs(z_l[x[1]][1][x[2]]/z_l[x[0]][1][x[2]])])
        r1 = lcm([x.denominator() for x in z_l[p[0]][1]])
        B[z_l[p[0]][0]] *= r1
        r2 = 1/(gcd(B[z_l[p[0]][0]])/gcd(B[z_l[p[0]][0]].list() + [B[z_l[p[1]][0]][p[2]]]))
        B[z_l[p[0]][0]] *= r2
        if r1*r2 != 1:
            reduction_steps += [[1,r1*r2,z_l[p[0]][0]]]
            if say:
                html(r"The operation ${0}\cdot R_{1}$ gives ${2}$. <br>".format(self.say_number(r1*r2),latex(z_l[p[0]][0]+1),latex(B)))
        r3 = -B[z_l[p[1]][0]][p[2]] / B[z_l[p[0]][0]][p[2]]
        if r3 != 0:
            reduction_steps += [[0,z_l[p[1]][0],r3,z_l[p[0]][0]]]
            B[z_l[p[1]][0]] += r3 * B[z_l[p[0]][0]]
            if say:
                html(r"The operation $R_{3}+{0}\cdot R_{1}$ gives ${2}$. <br>".format(self.say_number(r3),latex(z_l[p[0]][0]+1),latex(B),self.say_number(z_l[p[1]][0]+1)))
        return self.reduction_step(B,ev,ev_mult,reduction_steps,say,start_index)
    
    def say_gen_sol(self):
        ev_list = []
        mult = 1
        ev_list = []
        ev0 = self.J[0][0]
        for i in range(1,3):
            ev1 = self.J[i][i]
            if ev0 == ev1:
                mult += 1
            else:
                ev_list.append([ev0,mult])
                ev0 = ev1
                mult = 1
        ev_list.append([ev0,mult])
        fsl = []
        fin = 0
        for i in range(len(ev_list)):
            evects = self.find_evects(ev_list[i][0],ev_list[i][1],True,i)
            if ev_list[i][1] > 1 and rank(self.A-ev_list[i][0]) == 2:
                fsl += [r"c_{0}e^{{{1}t}}\mathbf u_{2} + c_{3}e^{{{1}t}}(\mathbf v_{2} + t\mathbf u_{2})".format(latex(fin+1),latex(ev_list[i][0]),latex(i+1),latex(fin+2))]
            else:
                fsl += [" + ".join([r"c_{0}e^{{{1}t}}\mathbf u_{2}".format(latex(fin+1+j),latex(ev_list[i][0]),latex(fin+1+j)) for j in range(ev_list[i][1])])]
            fin += ev_list[i][1]
        html(r"So a general solution is $$\mathbf x(t)={0}$$".format(" + ".join(fsl)))