import numpy
# site.addsitedir('/scratch/NICO/dadi-1.6.3_modif')  # Seems to be the important line!
import dadi
from dadi import Numerics, PhiManip, Integration, Inference
# import Inference
from dadi.Spectrum_mod import Spectrum


###########################################################
#
#         |
#         |
#     --------
#     |  <->  |
#   -----     |
#   |<-->|<-->|
#   |<-->|<-->|
#   |<-->|<-->| 
#  1c   1    2
#   <-------->
############################################################


def SIC_M12contemp(params, ns, pts):
    N1, N2a, N2, N3, T2, T1, m12, m21, m21c, m1c2, m21i, m1i2, m1i1c, m1c1i = params
    xx = dadi.Numerics.default_grid(pts)
    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    ## allow drift to occur along each of these branches
    phi = dadi.Integration.two_pops(phi, xx, T2, nu1=N1, nu2=N2a, m12=m12, m21=m21)

    phi = dadi.PhiManip.phi_2D_to_3D_split_2(xx, phi)  #split between 2 and 3

    phi = dadi.Integration.three_pops(phi, xx, T1, nu1=N1, nu2=N2, nu3=N3, m12=m21i, m13=m21c, m21=m1i2, m23=m1i1c, m31=m1c2, m32=m1c1i)
    ## simulate the fs
    fs = dadi.Spectrum.from_phi(phi, ns, (xx, xx, xx), pop_ids=['2', '1i', '1c'])
    return fs


best, opt = '', -100000.0
for jj in range(0, 10):

    func = SIC_M12contemp
    dd = dadi.Misc.make_data_dict('../SFS_3pops.xls')
    data = dadi.Spectrum.from_data_dict(dd, ['2', '1i', '1c'], [20, 20, 20], polarized=True)
    ns = data.sample_sizes
    print
    ns
    pts_l = [ns[0] + 5, ns[0] + 10, ns[0] + 15]

    # nu1, nu2, nu1fa, nu1fb, nu2f, TS1, TS2
    params =      [1,       1,    1,    1,    1,    1,    1,   1,     0.1,   0.1,  0.1, 0.1, 1, 1]
    upper_bound = [10,     10,   10,   10,   10,   10,   10,  10,    1,  1, 1, 1, 10, 10]
    lower_bound = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.001, 0.001, 0.001, 0.001, 0.01, 0.01]
    func_ex = dadi.Numerics.make_extrap_log_func(
        func)  # Makde the extrapolating version of our demographic model function.
    model = func_ex(params, ns, pts_l)  # Calculate the model AFS.
    ll_model = dadi.Inference.ll_multinom(model, data)  # Likelihood of the data given the model AFS.
    theta = dadi.Inference.optimal_sfs_scaling(model, data)  # The optimal value of theta given the model.
    print
    'Model log-likelihood:', ll_model, 'and theta', theta
    p0 = dadi.Misc.perturb_params(params, fold=1, lower_bound=lower_bound, upper_bound=upper_bound)

    popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l, lower_bound=lower_bound, upper_bound=upper_bound,
                                       maxiter=10000, verbose=1)

    model = func_ex(popt, ns, pts_l)
    ll_opt = dadi.Inference.ll_multinom(model, data)
    theta = dadi.Inference.optimal_sfs_scaling(model, data)
    AIC = 2 * len(params) - 2 * ll_opt
    print
    'Optimized log-likelihood:', ll_opt, 'AIC:', AIC, 'Theta:', theta
    print
    'Parameters:', popt
    fout = open('SIc_AllM.xls', 'a')
    fout.write('1st_anneal_optimization\t3pops_SIC_Mallways\t' + str(ll_opt) + '\t' + str(theta) + '\t' + str(popt).replace('\n',
                                                                                                                   '') + '\t' + str(
        AIC) + '\n')

    fout.close()

    if float(ll_opt) > opt:
        opt = ll_opt
        best = str(popt).replace('\n', '')

print
opt, best
newparm = []
print
best.split(' ')
for x in best.split(' '):
    if x != '' and x != ']' and x != '[':
        newparm.append(float(x.replace('[', '').replace(' ', '').replace(']', '').replace(',', '')))
    up, lw = [], []
    for x, y, z in zip(newparm, upper_bound, lower_bound):
        if x == y:
            up.append(x * 10)
        else:
            up.append(y)
        if x == z:
            lw.append(x / 10)
        else:
            lw.append(z)

print
'second optimization'
print
'new parameters', newparm
print
'new upper bounds', up
print
'new lower bounds', lw

for jj in range(0, 5):
    func = SIC_M12contemp
    dd = dadi.Misc.make_data_dict('../SFS_3pops.xls')
    data = dadi.Spectrum.from_data_dict(dd, ['2', '1i', '1c'], [20, 20, 20], polarized=True)
    ns = data.sample_sizes
    print
    ns
    pts_l = [ns[0] + 5, ns[0] + 10, ns[0] + 15]

    func_ex = dadi.Numerics.make_extrap_log_func(
        func)  # Makde the extrapolating version of our demographic model function.
    model = func_ex(newparm, ns, pts_l)  # Calculate the model AFS.
    ll_model = dadi.Inference.ll_multinom(model, data)  # Likelihood of the data given the model AFS.
    theta = dadi.Inference.optimal_sfs_scaling(model, data)  # The optimal value of theta given the model.
    print
    'Model log-likelihood:', ll_model, 'and theta', theta
    p0 = dadi.Misc.perturb_params(newparm, fold=1, lower_bound=lw, upper_bound=up)
    popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l, lower_bound=lw, upper_bound=up, maxiter=1000,
                                       verbose=1)
    model = func_ex(popt, ns, pts_l)
    ll_opt = dadi.Inference.ll_multinom(model, data)
    theta = dadi.Inference.optimal_sfs_scaling(model, data)

    AIC = 2 * len(params) - 2 * ll_opt
    print
    'Optimized log-likelihood:', ll_opt, 'AIC:', AIC, 'Theta:', theta
    print
    'Parameters:', popt

    fout = open('SIc_AllM.xls', 'a')
    fout.write('2nd_otimization\t3pops_SIC_Mallways\t' + str(ll_opt) + '\t' + str(theta) + '\t' + str(popt).replace('\n',
                                                                                                           '') + '\t' + str(
        AIC) + '\n')
    fout.close()