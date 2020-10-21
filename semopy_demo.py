# structural equation modeling to discover functional connectivity
# data generated with get_betas_for_tetrad.m

# install with:
# pip install semopy==2.0.0a4
# follow http://semopy.com/tutorial.html not https://pypi.org/project/semopy/

# old:
# had to vi ~/anaconda3/lib/python3.7/site-packages/semopy/model_base.py 
# and sys.path.append('/Users/momchil/anaconda3/lib/python3.7/site-packages/semopy')
# and also rename parser.py to e.g. parser_sem.py and import it like that...
from semopy import Model
from semopy import stats
import numpy as np
import os
import pandas as pd
from IPython import embed
import scipy.io
from scipy.stats import wishart


names = ['parallel', 'serial A', 'serial B', 'hybrid']

'''
descs = [""" Put ~ MTG
             VS ~ MTG
         """,
         """ Put ~ VS
             VS ~ MTG 
         """,
         """ Put ~ MTG + VS
         """,
         """ Put ~ MTG + VS
             VS ~ MTG
         """
         ]

descs = [""" Put ~ MTG + RPEpsi
             VS ~ MTG + RPE
             MTG ~ psi
         """,
         """ Put ~ VS + RPEpsi
             VS ~ MTG + RPE
             MTG ~ psi
         """,
         """ Put ~ MTG + VS + RPEpsi
             VS ~ MTG + RPE
             MTG ~ psi
         """,
         """ Put ~ MTG + VS + RPEpsi
             VS ~ RPE
             MTG ~ psi
         """
         ]
'''


'''
descs = [""" put ~ mtg + RPEpsi
             vs ~ mtg + RPE
             mtg ~ psi
             put =~ Put
             mtg =~ MTG
             vs =~ VS
         """,
         """ put ~ vs + RPEpsi
             vs ~ mtg + RPE
             mtg ~ psi
             put =~ Put
             mtg =~ MTG
             vs =~ VS
         """,
         """ put ~ mtg + vs + RPEpsi
             vs ~ mtg + RPE
             mtg ~ psi
             put =~ Put
             mtg =~ MTG
             vs =~ VS
         """,
         """ put ~ mtg + vs + RPEpsi
             vs ~ RPE
             mtg ~ psi
             put =~ Put
             mtg =~ MTG
             vs =~ VS
         """
         ]
'''

'''
descs = [""" Put ~ IFG
             VS ~ IFG
         """,
         """ Put ~ VS
             VS ~ IFG 
         """,
         """ Put ~ IFG + VS
         """,
         """ Put ~ IFG + VS
             VS ~ IFG
         """
         ]
'''

descs = [""" Put ~ RPEpsi + IFG
             VS ~ RPE + IFG
             IFG ~ psi
             Ins ~ psi
         """,
         """ Put ~ RPEpsi + VS
             VS ~ RPE + IFG 
             IFG ~ psi
             Ins ~ psi
         """,
         """ Put ~ RPEpsi + IFG + VS 
             VS ~ RPE
             IFG ~ psi
             Ins ~ psi
         """,
         """ Put ~ RPEpsi + IFG + VS
             VS ~ RPE + IFG
             IFG ~ psi
             Ins ~ psi
         """,
         """ Put ~ RPEpsi + IFG
             VS ~ RPE
             IFG ~ psi
             Ins ~ psi
         """,
         """ Put ~ RPEpsi + IFG
             VS ~ RPE
             IFG ~ psi
             Ins ~ psi
         """,
         """ Put ~ RPEpsi + IFG
             VS ~ RPE
             IFG ~ psi
             Ins ~ psi
         """,
         ]


# IFG and Ins
descs_IFG_Ins = [""" 
             IFG ~ psi
             Ins ~ psi
         """,
         """  
             IFG ~ psi + Ins
             Ins ~ psi
         """,
         """ 
             IFG ~ psi
             Ins ~ psi + IFG
         """,
         """ 
             IFG ~ psi + Ins
             Ins ~ psi + IFG
         """,
         ]


descs_ifg_ins = [""" 
             ifg ~ psi
             ins ~ psi
             ifg =~ IFG
             ins =~ Ins
         """,
         """  
             ifg ~ psi + ins
             ins ~ psi
             ifg =~ IFG
             ins =~ Ins
         """,
         """ 
             ifg ~ psi
             ins ~ psi + ifg
             ifg =~ IFG
             ins =~ Ins
         """,
         """ 
             ifg ~ psi + ins
             ins ~ psi + ifg
             ifg =~ IFG
             ins =~ Ins
         """,
         ]

# IFG and Ins

'''
descs = [""" 
             ifg ~ psi
             ins ~ psi
             ifg =~ IFG
             ins =~ Ins
         """,
         """  
             ifg ~ psi + ins
             ins ~ psi
             ifg =~ IFG
             ins =~ Ins
         """,
         """ 
             ifg ~ psi
             ins ~ psi + ifg
             ifg =~ IFG
             ins =~ Ins
         """,
         """ 
             ifg ~ psi + ins
             ins ~ psi + ifg
             ifg =~ IFG
             ins =~ Ins
         """,
         ]

'''


# VS and Put
descs_VS_Put = [""" Put ~ RPEpsi
             VS ~ RPE
         """,
         """ Put ~ RPEpsi + VS
             VS ~ RPE
         """,
         """ Put ~ RPEpsi + VS
             VS ~ RPE + Put
         """,
         ]

descs_vs_put = [""" put ~ RPEpsi
             vs ~ RPE
             vs =~ VS
             put =~ Put
         """,
         """ put ~ RPEpsi + vs
             vs ~ RPE
             vs =~ VS
             put =~ Put
         """,
         """ put ~ RPEpsi + vs
             vs ~ RPE + put
             vs =~ VS
             put =~ Put
         """,
         ]



'''
# VS and Put
descs = [""" put ~ RPEpsi
             vs ~ RPE
             put =~ Put
             vs =~ VS
         """,
         """ put ~ RPEpsi + vs
             vs ~ RPE
             put =~ Put
             vs =~ VS
         """,
         """ put ~ RPEpsi
             vs ~ RPE + put
             put =~ Put
             vs =~ VS
         """,
         """ put ~ RPEpsi + vs
             vs ~ RPE + put
             put =~ Put
             vs =~ VS
         """,
         ]

'''

'''
descs = [""" Put ~ RPEpsi
             VS ~ RPE
             IFG ~ psi
             Ins ~ psi
         """,
         """ Put ~ VS + IFG
             VS ~ RPE
             IFG ~ psi
             Ins ~ psi
         """,
         """ Put ~ VS + Ins
             VS ~ RPE
             IFG ~ psi
             Ins ~ psi
         """,
         """ Put ~ VS + Ins
             VS ~ RPE + IFG
             IFG ~ psi
             Ins ~ psi
         """,
         """ Put ~ RPEpsi
             VS ~ RPE
             IFG ~ Ins
             Ins ~ psi
         """,
         """ Put ~ VS + IFG
             VS ~ RPE
             IFG ~ Ins
             Ins ~ psi
         """,
         """ Put ~ VS + Ins
             VS ~ RPE
             IFG ~ Ins
             Ins ~ psi
         """,
         """ Put ~ VS + Ins
             VS ~ RPE + IFG
             IFG ~ Ins
             Ins ~ psi
         """,

         ]
'''

descs_all = [
         """ Put ~ RPEpsi + VS
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
         """ Put ~ RPEpsi + VS + IFG
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
         """ Put ~ RPEpsi + VS + Ins
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
         """ Put ~ RPEpsi + VS + IFG
             VS ~ RPE + Ins
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
         """ Put ~ RPEpsi + VS + Ins
             VS ~ RPE + IFG
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
         """ Put ~ RPEpsi + VS + Ins + IFG
             VS ~ RPE + IFG + Ins
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
        ]


descs_Cx_Str = [
         """ Put ~ RPEpsi + VS
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
         """ Put ~ RPEpsi + VS + IFG
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
         """ Put ~ RPEpsi + VS + Ins
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
         """ Put ~ RPEpsi + VS + Ins + IFG
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
         """ Put ~ RPEpsi + Ins
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
         """ Put ~ RPEpsi + VS + Ins
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi 
         """,
         """ Put ~ RPEpsi + Ins
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi 
         """,
        ]



descs_Sam = [
         """ Put ~ RPEpsi
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi
         """,
         """ Put ~ RPEpsi + Ins
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi
         """,
         """ Put ~ RPEpsi + IFG
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi
         """,
         """ Put ~ RPEpsi + Ins + IFG
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi
         """,
         """ Put ~ RPEpsi + VS
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi
         """,
         """ Put ~ RPEpsi + Ins + VS
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi
         """,
         """ Put ~ RPEpsi + IFG + VS
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi
         """,
         """ Put ~ RPEpsi + Ins + IFG + VS
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi
         """,

         """ Put ~ RPEpsi
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
         """ Put ~ RPEpsi + Ins
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
         """ Put ~ RPEpsi + IFG
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
         """ Put ~ RPEpsi + Ins + IFG
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
         """ Put ~ RPEpsi + VS
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
         """ Put ~ RPEpsi + Ins + VS
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
         """ Put ~ RPEpsi + IFG + VS
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
         """ Put ~ RPEpsi + Ins + IFG + VS
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
        ]


descs_confirm_0 = [
         """ Put ~ RPEpsi + Ins
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
         """ Put ~ RPEpsi + Ins
             VS ~ RPE
             IFG ~ psi + Ins
             Ins ~ psi
         """,
         ]


descs_confirm_1 = [
         """ Put ~ RPEpsi + Ins
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
         """ Put ~ RPEpsi + Ins + VS
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
         """ Put ~ RPEpsi
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
         """ Put ~ RPEpsi + IFG
             VS ~ RPE
             IFG ~ psi 
             Ins ~ psi + IFG
         """,
         ]

descs_confirm_2 = [
         """ Put ~ RPEpsi + Ins
             VS ~ RPE
             IFG ~ psi + Ins
             Ins ~ psi
         """,
         """ Put ~ RPEpsi + Ins + VS
             VS ~ RPE
             IFG ~ psi + Ins
             Ins ~ psi
         """,
         """ Put ~ RPEpsi
             VS ~ RPE
             IFG ~ psi + Ins
             Ins ~ psi
         """,
         """ Put ~ RPEpsi + IFG
             VS ~ RPE
             IFG ~ psi + Ins
             Ins ~ psi
         """,
         ]



def gen_descs_helper(cur, cand, ix, descs):

    if ix == len(cand):
        desc = ''
        for k, v in cur.items():
            desc += k + ' ~ ' + ' + '.join(v) + '   \n'

        #if 'Ins' in cur['IFG'] and 'IFG' in cur['Ins']:
        #    return # TODO hack ensure DAG
        descs.append(desc)
        return

    # without A -> B
    gen_descs_helper(cur, cand, ix + 1, descs)

    A = cand[ix][0]
    B = cand[ix][1]
    cur[B].append(A)
    # with A -> B
    gen_descs_helper(cur, cand, ix + 1, descs)
    cur[B].remove(A)

def gen_descs(base, cand):
    descs = []
    gen_descs_helper(base, cand, 0, descs)
    return descs

base = {'Put': ['RPEpsi'], 'VS': ['RPE'], 'IFG': ['psi'], 'Ins': ['psi']}
cand = [('Ins', 'Put'), ('IFG', 'Put'), ('VS', 'Put'), ('IFG', 'Ins'), ('Ins', 'IFG')]

# auto generate models
#
#descs = gen_descs(base, cand)

#print(descs)
#embed()


#descs = descs_VS_Put
#descs = descs_IFG_Ins
#descs = descs_all
#descs = descs_Sam
#descs = descs_Cx_Str

descs = descs_confirm_2



for dirname in ['S1_feedback_onset', 'S1_trial_onset']:
#for dirname in ['feedback_onset', 'trial_onset', 'residuals']:
#for dirname in ['residuals']:

    print(dirname)

    files = os.listdir(dirname)
    files = [f for f in files if f.endswith('.txt') ]

    n_subjects = len(files) # num subjects
    n_models = len(descs) # num models

    logliks = np.zeros((n_subjects, n_models)) # likelihoods
    lmes = np.zeros((n_subjects, n_models)) # LMEs = -0.5 * BICs
    bics = np.zeros((n_subjects, n_models)) # BICs
    ks = np.zeros((n_subjects, n_models)) # # params
    ns = np.zeros((n_subjects, n_models)) # # data points

    for i in range(n_subjects):

        for j in range(n_models):

            model = Model(descs[j])

            filepath = os.path.join(dirname, files[i])
            data = pd.read_csv(filepath, sep='\t')

            opt_res = model.fit(data)

            # from /Users/momchil/anaconda3/lib/python3.7/site-packages/semopy/stats.py: calc_bic()
            # WARNING: all of these are up to a proportionality constant that DIFFERS ACROSS SUBJECTS => do not use for BMS
            #
            '''
            logliks[i,j] = stats.calc_likelihood(model) # note up to proportionality constant, b/c w.r.t. saturated model, but that's the same for all models so it's fine
            ks[i,j], ns[i,j] = len(model.param_vals), model.mx_data.shape[0]
            bic = stats.calc_bic(model)
            lmes[i,j] = -0.5 * bic
            '''

            # calculate likelihood manually, following https://en.wikipedia.org/wiki/Wishart_distribution
            # and https://math.stackexchange.com/questions/2803164/degrees-of-freedom-in-a-wishart-distribution
            # and /Users/momchil/anaconda3/lib/python3.7/site-packages/semopy/stats.py: obj_mlw()
            # and /Users/momchil/anaconda3/lib/python3.7/site-packages/semopy/stats.py: calc_bic()
            # and /Users/momchil/anaconda3/lib/python3.7/site-packages/semopy/stats.py: calc_likelihood()
            model.update_matrices(model.param_vals)
            S = model.mx_cov # empirical covariance matrix
            sigma, _ = model.calc_sigma() # model covariance matrix
            p = sigma.shape[0] # # of variables
            n = model.mx_data.shape[0] # dof = # of data points (b/c Wishart is generalization of Chi2 to MVN => sum of n squared i.i.d. MVNs)
            k = len(model.param_vals) # # of free params

            loglik = wishart.logpdf(S, df=n, scale=sigma)
            bic = k*np.log(n) - 2*loglik
            lme = -0.5 * bic 

            logliks[i,j] = loglik
            ks[i,j] = k 
            ns[i,j] = n
            bics[i,j] = bic 
            lmes[i,j] = lme

            print('subj=', i, ' mod=', j, filepath, ' k=', k, ' n=', n, ' loglik=', logliks[i,j])
            print(model.mx_cov.shape)

    print(lmes.shape)

    d = {'lmes': lmes,
         'bics': bics,
         'logliks': logliks,
         'ks': ks,
         'ns': ns,
         'files': files,
         'descs': descs}

    scipy.io.savemat('semopy_Sam_%s_lmes.mat' % dirname, d)



