from tdc import Oracle
from scscore_standalone_model_numpy import SCScorer
import os
from RAscore import RAscore_XGB
import sys

project_root = '/c7/home/margonza/scscore'

oracle = Oracle(name = 'SA')

input_file = sys.argv[1]
generated_comps = open(input_file, 'r')
name = input_file.split('.')[0]
generated_comps_lines = generated_comps.readlines()

scores = open('scores_%s.txt' % name,'w')
scores.write("smiles\tsascore\tscscore\tRAscore\n")
model =SCScorer()
model.restore(os.path.join(project_root, 'models', 'full_reaxys_model_1024uint8', 'model.ckpt-10654.as_numpy.json.gz'))

xgb_scorer = RAscore_XGB.RAScorerXGB()
for line in generated_comps_lines:
    try:
        smi = line.rstrip()
        print(smi)
        sa = oracle(smi)
        smi2, sco = model.get_score_from_smi(smi)
        RA = xgb_scorer.predict(smi)
        scores.write("%s\t%0.3f\t%0.3f\t%0.3f\n" % (smi, sa, sco, RA))
    except:
        RA= 'NA'
        scores.write("%s\t%0.3f\t%0.3f\t%s\n" % (smi, sa, sco, RA))
        print('error on %s' %smi)
                                                  