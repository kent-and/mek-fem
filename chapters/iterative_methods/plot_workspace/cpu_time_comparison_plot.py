import pylab
from numpy import array

Ns= array([32, 64, 128, 256, 512, 1024])**2
lu_time = [0.00400400161743, 0.0155980587006, 0.0793328285217, 0.482635974884, 3.27159714699, 28.2356829643]
cg_time = [ 0.00286412239075, 0.015839099884, 0.119037151337, 0.85956287384, 9.34541797638, 83.0514290333]
cgilu_time =[ 0.00153207778931,  0.00740194320679, 0.0501880645752, 0.371992111206,3.67462015152,  29.5216519833]
cgamg_time= [0.00507497787476,  0.0127878189087, 0.0421559810638, 0.175911903381, 0.791210174561, 3.56511092186]


pylab.plot(Ns, lu_time)
pylab.plot(Ns, cg_time)
pylab.plot(Ns, cgilu_time)
pylab.plot(Ns, cgamg_time)  
pylab.xlabel('Degrees of freedom')
pylab.ylabel('Time(sec)')
pylab.legend(["lu", "cg", "cg/ilu", "cg/amg"])
pylab.savefig("cpu_time_comparison_plot.pdf")
pylab.show()

#pylab.loglog(Ns, lu_time)
#pylab.loglog(Ns, cg_time)
#pylab.loglog(Ns, cgilu_time)
#pylab.loglog(Ns, cgamg_time)
#pylab.legend(["lu", "cg", "cg/ilu", "cg/amg"])
#pylab.savefig('tmp_cpu.pdf')
#pylab.show()