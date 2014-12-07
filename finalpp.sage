min_m = 3
max_m = 5

max_step = 4
max_abs = 4

plist = ["Random", "Heat Equation"]

@interact
def _f(psel = Selector(plist, label = "PDE type:", selector_type = "button"), regen = Button(text="Regenerate Problem", default=True, value=True, label = "")):
    if psel == plist[0]:
        psel = choice(plist[1:])
    if psel == plist[1]:
    
        blist = ["Random", "Temperatures fixed at ends"]
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
            
            if bsel == blist[1]:
            
                bc = [randint(-max_abs,max_abs), randint(-max_abs,max_abs)]
                
                if bc == [0,0]:
                    problem = HHE1d(f = f, bc = bc, alpha2 = alpha2)
                else:
                    problem = NHHE1d(f = f, bc = bc, alpha2 = alpha2)
            
                html(r"Find the {rsel} to the problem".format(rsel = rsel))
                html(problem.say_eqs())
                
                if lsel == llist[1]:
                    html("Where $f(x)$ has graph<br>")
                    show(f.plot_function())
                else:
                    html("Where $f(x)$ has formula" + f.say_function())