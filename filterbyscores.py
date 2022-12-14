import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from os.path import join
from rdkit.Chem import PandasTools, Descriptors
from janitor import chemistry
from rdkit.Chem import Descriptors
import sys

name = sys.argv[1]
_file = sys.argv[2]
#_file = 'decoys_%s_scores.txt' % name
scores = pd.read_csv(_file, sep = '\t')
print(scores.shape[0])
scores = scores.drop_duplicates(subset='smiles', keep="last")

# plot distribution
plt.figure()
sns.scatterplot(data=scores, x="sascore", y="RAscore")
plt.xlabel('Sa score')
plt.ylabel('RA score')
plt.savefig('%s_saRA.png' % name)
plt.show()

plt.figure()
sns.scatterplot(data=scores, x="scscore", y="RAscore")
plt.xlabel('SC score')
plt.ylabel('RA score')
plt.savefig('%s_scRA.png' % name)
plt.show()

plt.figure()
sns.scatterplot(data=scores, x="sascore", y="scscore")
plt.xlabel('Sa score')
plt.ylabel('Sc score')
plt.savefig('%sa_sc.png' % name)
plt.show()

# get comps with RAscore higher than 0.7
scores_RA = scores.loc[(scores['RAscore'] >= 0.7)]
print(scores_RA.shape[0])
# get comps with sascore lower than 3.5
scores_sascore = scores.loc[(scores['sascore'] <= 3.5)]
print(scores_sascore.shape[0])
# get comps with scscore lower than 3.5
scores_scscore = scores.loc[(scores['scscore'] <= 3.5)]
print(scores_scscore.shape[0])


sgcs = pd.concat([scores_RA,scores_sascore,scores_scscore])
sgcs = sgcs.drop_duplicates(subset='smiles', keep="last")

print(sgcs.shape[0])
print(sgcs.head(3))
PandasTools.AddMoleculeColumnToFrame(sgcs, smilesCol='smiles')
sgcs = sgcs.mask(sgcs.astype(object).eq('None')).dropna()

# remove molecules with MW < 180
sgcs = sgcs.add_column('MolWt', [Descriptors.MolWt(mol) for mol in sgcs.ROMol])
selected_generated_compounds = sgcs.loc[(sgcs['MolWt'] >= 180)]
print(selected_generated_compounds.shape[0])

selected_generated_compounds['smiles'].to_csv('%s_filtered_decoys.smi'% (name), index=False)


