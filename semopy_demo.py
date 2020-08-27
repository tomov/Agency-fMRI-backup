# install with:
# pip install semopy==2.0.0a4
# follow http://semopy.com/tutorial.html not https://pypi.org/project/semopy/

# old:
# had to vi ~/anaconda3/lib/python3.7/site-packages/semopy/model_base.py 
# and sys.path.append('/Users/momchil/anaconda3/lib/python3.7/site-packages/semopy')
# and also rename parser.py to e.g. parser_sem.py and import it like that...
from semopy import Model

from IPython import embed


mod = """ Put ~ MTG + VS
          VS ~ MTG 
          MTG ~ VS
      """
model = Model(mod)

import pandas as pd
data = pd.read_csv('feedback_onset/feedback_onset_SPMsubj1.txt', sep='\t')
#model.load_data(data)

from semopy import stats

opt_res = model.fit(data)
estimates = model.inspect()

print(stats.calc_bic(model))

embed()

