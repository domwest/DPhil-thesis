import tensorflow as tf
# Make sure the GPU is enabled
#assert tf.config.list_physical_devices('GPU'), 'Start the colab kernel with GPU: Runtime -> Change runtime type -> GPU'

# !pip install kipoiseq==0.5.2 --quiet > /dev/null
# You can ignore the pyYAML error

"""### Imports"""

# Commented out IPython magic to ensure Python compatibility.
import tensorflow_hub as hub
import joblib ##
import gzip
import kipoiseq
from kipoiseq import Interval
import pyfaidx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt ##
import matplotlib as mpl
import seaborn as sns ##
import sys
print(sys.argv[0]) # prints step1_4.py
print(sys.argv[1]) # prints var1
print(sys.argv[2]) # prints var2
print(sys.argv[3]) # prints var3
print(sys.argv[4]) # prints var4
print(sys.argv[5]) # prints var5

# %matplotlib inline
# %config InlineBackend.figure_format = 'retina'

transform_path = 'gs://dm-enformer/models/enformer.finetuned.SAD.robustscaler-PCA500-robustscaler.transform.pkl'
model_path = 'https://tfhub.dev/deepmind/enformer/1'
fasta_file = '/well/PROCARDIS/domwest/deepmind/data/genome.fa'
clinvar_vcf = '/well/PROCARDIS/domwest/deepmind/data/clinvar.vcf.gz'

# Download targets from Basenji2 dataset
# Cite: Kelley et al Cross-species regulatory sequence activity prediction. PLoS Comput. Biol. 16, e1008050 (2020).
targets_txt = '/well/PROCARDIS/domwest/deepmind/data/targets_human.txt'
df_targets = pd.read_csv(targets_txt, sep='\t')
df_targets.head(3)

df_targets_filt = df_targets[(df_targets['description'].str.contains('DNASE:')) | (df_targets['description'].str.contains('ATAC:'))]

# df_targets['index'][(df_targets['description'].str.contains('DNASE:'))&((df_targets['description'].str.contains('cardiac muscle cell'))|(df_targets['description'].str.contains('Cardiac Myocyte'))|(df_targets['description'].str.contains('heart left ventricle')))].unique()
df_targets_filt = df_targets_filt[~(df_targets_filt.index.isin([85, 198, 415, 586, 642]))]

pyfaidx.Faidx(fasta_file)

"""### Code (double click on the title to show the code)"""

# @title `Enformer`, `EnformerScoreVariantsNormalized`, `EnformerScoreVariantsPCANormalized`,
SEQUENCE_LENGTH = 393216

class Enformer:
  def __init__(self, tfhub_url):
    self._model = hub.load(tfhub_url).model
  def predict_on_batch(self, inputs):
    predictions = self._model.predict_on_batch(inputs)
    return {k: v.numpy() for k, v in predictions.items()}
  @tf.function
  def contribution_input_grad(self, input_sequence,
                              target_mask, output_head='human'):
    input_sequence = input_sequence[tf.newaxis]
    target_mask_mass = tf.reduce_sum(target_mask)
    with tf.GradientTape() as tape:
      tape.watch(input_sequence)
      prediction = tf.reduce_sum(
          target_mask[tf.newaxis] *
          self._model.predict_on_batch(input_sequence)[output_head]) / target_mask_mass
    input_grad = tape.gradient(prediction, input_sequence) * input_sequence
    input_grad = tf.squeeze(input_grad, axis=0)
    return tf.reduce_sum(input_grad, axis=-1)

class EnformerScoreVariantsRaw:
  def __init__(self, tfhub_url, organism='human'):
    self._model = Enformer(tfhub_url)
    self._organism = organism
  def predict_on_batch(self, inputs):
    ref_prediction = self._model.predict_on_batch(inputs['ref'])[self._organism]
    alt_prediction = self._model.predict_on_batch(inputs['alt'])[self._organism]
    return alt_prediction.mean(axis=1) - ref_prediction.mean(axis=1)

class EnformerScoreVariantsNormalized:
  def __init__(self, tfhub_url, transform_pkl_path,
               organism='human'):
    assert organism == 'human', 'Transforms only compatible with organism=human'
    self._model = EnformerScoreVariantsRaw(tfhub_url, organism)
    with tf.io.gfile.GFile(transform_pkl_path, 'rb') as f:
      transform_pipeline = joblib.load(f)
    self._transform = transform_pipeline.steps[0][1]  # StandardScaler.    
  def predict_on_batch(self, inputs):
    scores = self._model.predict_on_batch(inputs)
    return self._transform.transform(scores)

class EnformerScoreVariantsPCANormalized:
  def __init__(self, tfhub_url, transform_pkl_path,
               organism='human', num_top_features=500):
    self._model = EnformerScoreVariantsRaw(tfhub_url, organism)
    with tf.io.gfile.GFile(transform_pkl_path, 'rb') as f:
      self._transform = joblib.load(f)
    self._num_top_features = num_top_features
  def predict_on_batch(self, inputs):
    scores = self._model.predict_on_batch(inputs)
    return self._transform.transform(scores)[:, :self._num_top_features]


# TODO(avsec): Add feature description: Either PCX, or full names.

# @title `variant_centered_sequences`

class FastaStringExtractor:
    def __init__(self, fasta_file):
        self.fasta = pyfaidx.Fasta(fasta_file)
        self._chromosome_sizes = {k: len(v) for k, v in self.fasta.items()}
    def extract(self, interval: Interval, **kwargs) -> str:
        # Truncate interval if it extends beyond the chromosome lengths.
        chromosome_length = self._chromosome_sizes[interval.chrom]
        trimmed_interval = Interval(interval.chrom,
                                    max(interval.start, 0),
                                    min(interval.end, chromosome_length),
                                    )
        # pyfaidx wants a 1-based interval
        sequence = str(self.fasta.get_seq(trimmed_interval.chrom,
                                          trimmed_interval.start + 1,
                                          trimmed_interval.stop).seq).upper()
        # Fill truncated values with N's.
        pad_upstream = 'N' * max(-interval.start, 0)
        pad_downstream = 'N' * max(interval.end - chromosome_length, 0)
        return pad_upstream + sequence + pad_downstream
    def close(self):
        return self.fasta.close()

def variant_generator(vcf_file, gzipped=False):
  """Yields a kipoiseq.dataclasses.Variant for each row in VCF file."""
  def _open(file):
    return gzip.open(vcf_file, 'rt') if gzipped else open(vcf_file)
  with _open(vcf_file) as f:
    for line in f:
      if line.startswith('#'):
        continue
      chrom, pos, id, ref, alt_list = line.split('\t')[:5]
      # Split ALT alleles and return individual variants as output.
      for alt in alt_list.split(','):
        yield kipoiseq.dataclasses.Variant(chrom=chrom, pos=pos,
                                           ref=ref, alt=alt, id=id)

def one_hot_encode(sequence):
  return kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.float32)

def variant_centered_sequences(vcf_file, sequence_length, gzipped=False,
                               chr_prefix=''):
  seq_extractor = kipoiseq.extractors.VariantSeqExtractor(
    reference_sequence=FastaStringExtractor(fasta_file))
  for variant in variant_generator(vcf_file, gzipped=gzipped):
    interval = Interval(chr_prefix + variant.chrom,
                        variant.pos, variant.pos)
    interval = interval.resize(sequence_length)
    center = interval.center() - interval.start
    reference = seq_extractor.extract(interval, [], anchor=center)
    alternate = seq_extractor.extract(interval, [variant], anchor=center)
    yield {'inputs': {'ref': one_hot_encode(reference),
                      'alt': one_hot_encode(alternate)},
           'metadata': {'chrom': chr_prefix + variant.chrom,
                        'pos': variant.pos,
                        'id': variant.id,
                        'ref': variant.ref,
                        'alt': variant.alt}}

# @title `plot_tracks`

def plot_tracks(tracks, interval, classif, chrom, pos, ref, alt, var_type, height=1.5):
  fig, axes = plt.subplots(len(tracks), 1, figsize=(20, height * len(tracks)), sharex=True)
  for ax, (title, y) in zip(axes, tracks.items()):
    ax.fill_between(np.linspace(interval.start, interval.end, num=len(y)), y)
    ax.set_title(title)
    sns.despine(top=True, right=True, bottom=True)
  ax.set_xlabel(str(interval))
  plt.tight_layout()
  plt.savefig('/well/PROCARDIS/domwest/deepmind/plot_tracks_extreme/' + classif + '_' + chrom + '_' + str(pos) + '_' + ref + '_' + alt + '_' + var_type + '_tracks.png')

"""## Make predictions for a genetic sequenece"""

model = Enformer(model_path)

fasta_extractor = FastaStringExtractor(fasta_file)

# @title Score the variant
def score_the_variant(chrom, pos, ref, alt, rsid, classifier, index, var_type):
  variant = kipoiseq.Variant(chrom, pos, ref, alt, id=rsid)
  #
  print(chrom + '_' + str(pos) + '_' + ref + '_' + alt)
  classif = classifier.replace(' ','_')
  classif = classif.replace(':','_')
  #
  # Center the interval at the variant
  interval = kipoiseq.Interval(variant.chrom, variant.start, variant.start).resize(SEQUENCE_LENGTH)
  seq_extractor = kipoiseq.extractors.VariantSeqExtractor(reference_sequence=fasta_extractor)
  center = interval.center() - interval.start
  #
  reference = seq_extractor.extract(interval, [], anchor=center)
  alternate = seq_extractor.extract(interval, [variant], anchor=center)
  #
  # Make predictions for the refernece and alternate allele
  reference_prediction = model.predict_on_batch(one_hot_encode(reference)[np.newaxis])['human'][0]
  alternate_prediction = model.predict_on_batch(one_hot_encode(alternate)[np.newaxis])['human'][0]
  #print(len(reference_prediction))
  #print(len(alternate_prediction))
  #
  middle_ref = reference_prediction[:, index][len(reference_prediction[:, index])//2-1:len(reference_prediction[:, index])//2+1]
  # print(middle_ref)
  # print(np.mean((middle_ref)))
  middle_alt = alternate_prediction[:, index][len(alternate_prediction[:, index])//2-1:len(alternate_prediction[:, index])//2+1]
  # print(middle_alt)
  # print(np.mean((middle_alt)))
  ref_minus_alt = alternate_prediction[:, index] - reference_prediction[:, index]
  middle_ds = ref_minus_alt[len(ref_minus_alt)//2-1:len(ref_minus_alt)//2+1]
  # print(middle_ds)
  # print(np.mean((middle_ds)))
  #
  # reference_prediction_df = pd.DataFrame()
  # reference_prediction_df[classifier] = ''
  # reference_prediction_df[classifier] = reference_prediction[:, index]
  # reference_prediction_df.to_csv('/content/ref_prediction_' + classif + '_' + chrom + '_' + str(pos) + '_' + ref + '_' + alt + '.csv', index=False)
  # alternate_prediction_df = pd.DataFrame()
  # alternate_prediction_df[classifier] = ''
  # alternate_prediction_df[classifier] = alternate_prediction[:, index]
  # alternate_prediction_df.to_csv('/content/alt_prediction_' + classif + '_' + chrom + '_' + str(pos) + '_' + ref + '_' + alt + '.csv', index=False)
  # ref_minus_alt = alternate_prediction[:, index] - reference_prediction[:, index]
  # ref_minus_alt_df = pd.DataFrame()
  # ref_minus_alt_df[classifier] = ''
  # ref_minus_alt_df[classifier] = ref_minus_alt
  # ref_minus_alt_df.to_csv('/content/ds_prediction_' + classif + '_' + chrom + '_' + str(pos) + '_' + ref + '_' + alt + '.csv', index=False)
  #
  # @title Visualize some tracks
  variant_track = np.zeros_like(reference_prediction[:, 0], dtype=bool)
  variant_track[variant_track.shape[0] // 2] = True
  tracks = {'variant': variant_track,
            classifier + ' ref': reference_prediction[:, index],
            classifier + ' alt': alternate_prediction[:, index],
            classifier + ' alt-ref': alternate_prediction[:, index] - reference_prediction[:, index],
            }
  #
  classifier = classifier.replace(' ','_')
  plot_tracks(tracks, interval.resize(reference_prediction.shape[0] * 128), classif, chrom, pos, ref, alt, var_type, height=1)
  #
  pred_df = pd.DataFrame()
  pred_df['classifier'] = ''
  pred_df['chr'] = ''
  pred_df['pos'] = ''
  pred_df['ref'] = ''
  pred_df['alt'] = ''
  pred_df['ref_score'] = ''
  pred_df['alt_score'] = ''
  pred_df['damage_score'] = ''
  pred_df['var_type'] = ''
  pred_df.loc[len(pred_df)] = [classifier, chrom, pos, ref, alt, np.mean((middle_ref)), np.mean((middle_alt)), np.mean((middle_ds)), var_type]
  pred_df['snp_test_id'] = pred_df['chr'].str.replace('chr','') + '_' + pred_df['pos'].astype(str) + '_' + pred_df['ref'] + '_' + pred_df['alt']
  #
  pred_df.to_csv('/well/PROCARDIS/domwest/deepmind/bi_allelic_extreme_prediction_dfs/' + classifier + '_' + chrom + '_' + str(pos) + '_' + ref + '_' + alt + '.csv')

# for index, row in df_targets_filt.iterrows():
#   # print(index)
#   # print(row['description'].replace(' ','_'))
#   score_the_variant(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[1], row['description'], index, 'bi_allelic_extreme')

score_the_variant('chr1',2212668,'A','G','rs2503715', sys.argv[2], sys.argv[1], 'bi_allelic_extreme') 
score_the_variant('chr5',57700364,'A','T','rs4699910', sys.argv[2], sys.argv[1], 'bi_allelic_extreme') 
score_the_variant('chr9',131640601,'ACT','A','rs66925847', sys.argv[2], sys.argv[1], 'bi_allelic_extreme') 
score_the_variant('chr5',57700364,'A','T','rs4699910', sys.argv[2], sys.argv[1], 'bi_allelic_extreme') 
score_the_variant('chr13',110388334,'G','A','rs11838776', sys.argv[2], sys.argv[1], 'bi_allelic_extreme') 
score_the_variant('chr11',14288736,'G','A','rs11023178', sys.argv[2], sys.argv[1], 'bi_allelic_extreme') 
score_the_variant('chr1',2212668,'A','G','rs2503715', sys.argv[2], sys.argv[1], 'bi_allelic_extreme') 
score_the_variant('chr6',149766491,'C','A','rs1112729', sys.argv[2], sys.argv[1], 'bi_allelic_extreme') 
score_the_variant('chr12',120230731,'G','A','rs116904997', sys.argv[2], sys.argv[1], 'bi_allelic_extreme') 
score_the_variant('chr14',71693604,'T','A','rs7148679', sys.argv[2], sys.argv[1], 'bi_allelic_extreme') 
score_the_variant('chr14',71693604,'T','C','rs7148679', sys.argv[2], sys.argv[1], 'bi_allelic_extreme') 
score_the_variant('chr19',45810500,'G','C','rs10402263', sys.argv[2], sys.argv[1], 'bi_allelic_extreme') 
score_the_variant('chr12',120230731,'G','A','rs116904997', sys.argv[2], sys.argv[1], 'bi_allelic_extreme') 
score_the_variant('chr6',118370989,'G','C','rs10457327', sys.argv[2], sys.argv[1], 'bi_allelic_extreme') 
score_the_variant('chr12',120230731,'G','A','rs116904997', sys.argv[2], sys.argv[1], 'bi_allelic_extreme') 
score_the_variant('chr18',35090057,'T','C','rs476348', sys.argv[2], sys.argv[1], 'bi_allelic_extreme') 
score_the_variant('chr12',120230731,'G','A','rs116904997', sys.argv[2], sys.argv[1], 'bi_allelic_extreme') 
score_the_variant('chr6',118370989,'G','C','rs10457327', sys.argv[2], sys.argv[1], 'bi_allelic_extreme') 
score_the_variant('chr18',35090057,'T','C','rs476348', sys.argv[2], sys.argv[1], 'bi_allelic_extreme') 

# #REgulamentary info

# """## Score variants in a VCF file

# ### Report top 20 PCs
# """

# enformer_score_variants = EnformerScoreVariantsPCANormalized(model_path, transform_path, num_top_features=20)

# # Score the first 5 variants from ClinVar
# # Lower-dimensional scores (20 PCs)
# it = variant_centered_sequences(clinvar_vcf, sequence_length=SEQUENCE_LENGTH,
#                                 gzipped=True, chr_prefix='chr')
# example_list = []
# for i, example in enumerate(it):
#   if i >= 5:
#     break
#   variant_scores = enformer_score_variants.predict_on_batch(
#       {k: v[tf.newaxis] for k,v in example['inputs'].items()})[0]
#   variant_scores = {f'PC{i}': score for i, score in enumerate(variant_scores)}
#   example_list.append({**example['metadata'],
#                        **variant_scores})
#   if i % 2 == 0:
#     print(f'Done {i}')
# df = pd.DataFrame(example_list)
# df

# """### Report all 5,313 features (z-score normalized)"""

# enformer_score_variants_all = EnformerScoreVariantsNormalized(model_path, transform_path)

# # Score the first 5 variants from ClinVar
# # All Scores
# it = variant_centered_sequences(clinvar_vcf, sequence_length=SEQUENCE_LENGTH,
#                                 gzipped=True, chr_prefix='chr')
# example_list = []
# for i, example in enumerate(it):
#   if i >= 5:
#     break
#   variant_scores = enformer_score_variants_all.predict_on_batch(
#       {k: v[tf.newaxis] for k,v in example['inputs'].items()})[0]
#   variant_scores = {f'{i}_{name[:20]}': score for i, (name, score) in enumerate(zip(df_targets.description, variant_scores))}
#   example_list.append({**example['metadata'],
#                        **variant_scores})
#   if i % 2 == 0:
#     print(f'Done {i}')
# df = pd.DataFrame(example_list)
# df
