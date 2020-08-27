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


names = ['parallel', 'serial', 'hybrid']

descs = [""" Put ~ MTG
             VS ~ MTG
         """,
         """ Put ~ VS
             VS ~ MTG 
         """,
         """ Put ~ MTG + VS
             VS ~ MTG
         """]

#dirname = 'feedback_onset';
#dirname = 'trial_onset';
dirname = 'residuals';
files = os.listdir(dirname)

n = len(files) # num subjects
m = len(descs) # num models

lme = np.zeros((n, m))

for i in range(n):

    for j in range(m):

        model = Model(descs[j])

        filepath = os.path.join(dirname, files[i])
        data = pd.read_csv(filepath, sep='\t')

        opt_res = model.fit(data)

        bic = stats.calc_bic(model)
        lme[i,j] = -0.5 * bic


        print('subj', i, ' mod ', j, filepath)


d = {'lme': lme,
     'files': files,
     'descs': descs}

scipy.io.savemat('semopy_%s_lmes.mat' % dirname, d)



