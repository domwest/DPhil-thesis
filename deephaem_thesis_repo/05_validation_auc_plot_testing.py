#!/usr/bin/env python
# coding: utf-8

# In[1]:


import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import h5py
from matplotlib import pyplot
from sklearn import preprocessing
from sklearn.model_selection import train_test_split

from tensorflow.keras.layers import Dense, Flatten, Dropout, Input, Conv1D, MaxPooling1D, BatchNormalization
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint, CSVLogger
from tensorflow.keras.initializers import GlorotUniform
from tensorflow.keras.models import Sequential
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.regularizers import l2
from tensorflow.keras import backend

import tensorflow as tf
import keras


# In[2]:


with tf.device("CPU"):
    labels = h5py.File('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/ML_training_sets/labels.h5', 'r')
    labels = labels['label'][:]
    labels.shape
#so this is whether or not a specific peak is present in the different cell_types

with tf.device("CPU"):
    data = h5py.File('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/ML_training_sets/sequences.h5', 'r')
    data = data['data'][:]
    data.shape
#so this is the actual hot encode format of the sequences (1000bp seqs) per peak


# In[3]:


with tf.device("CPU"):
    X_train, X_test, y_train, y_test = train_test_split(data,   labels, test_size=0.4, random_state=42, shuffle=True)
    X_val,   X_test, y_val,   y_test = train_test_split(X_test, y_test, test_size=0.5, random_state=42, shuffle=True)
    print("Train:     ", X_train.shape, y_train.shape)
    print("Validation:", X_val.shape,   y_val.shape)
    print("Test:      ", X_test.shape,  y_test.shape)


# In[4]:


with tf.device("CPU"):
    obj = sc.read_h5ad("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/cm.h5ad")
    df_train = pd.DataFrame([np.sum(y_train, axis=0)], index=['count'], columns=obj.var.index.tolist()).T
    df_train['set'] = 'train_label'
    df_val   = pd.DataFrame([np.sum(y_val, axis=0)],   index=['count'], columns=obj.var.index.tolist()).T
    df_val['set'] = 'val_label'
    df_test  = pd.DataFrame([np.sum(y_test, axis=0)],  index=['count'], columns=obj.var.index.tolist()).T
    df_test['set'] = 'test_label'
    df = pd.concat([df_train, df_val, df_test])
    df['cell_type'] = df.index.tolist()
    fig, ax = pyplot.subplots(figsize=(16,9))
    ax = sns.barplot(data=df, x="count", y="cell_type", hue="set", ax=ax)


# In[5]:


with tf.device("CPU"):
    batch_size = 100
    epochs     = 100

with tf.device("CPU"):
    input_shape  = (X_train.shape[-2], X_train.shape[-1])
    output_shape = labels.shape[-1]

with tf.device("CPU"):
    training_seqs   = tf.data.Dataset.from_tensor_slices((X_train, y_train)).batch(batch_size)
    validation_seqs = tf.data.Dataset.from_tensor_slices((X_val, y_val)).batch(batch_size)


# In[6]:


model = Sequential()
model.add(Input(shape=input_shape))
model.add(Conv1D(filters=4, kernel_size=(1), strides=(1), activation='relu', kernel_initializer=GlorotUniform(seed=42), padding="same"))
model.add(Dropout(0.2))
model.add(Conv1D(filters=300, kernel_size=(20), strides=(1), activation='relu', kernel_initializer=GlorotUniform(seed=42), padding="same"))
model.add(MaxPooling1D(pool_size=(3), strides=(3), padding="same"))
model.add(Dropout(0.2))
model.add(Conv1D(filters=600, kernel_size=(10), strides=1, activation='relu', kernel_initializer=GlorotUniform(seed=42), padding="same"))
model.add(MaxPooling1D(pool_size=(4), strides=(4), padding="same"))
model.add(Dropout(0.2))
model.add(Conv1D(filters=900, kernel_size=(8),  strides=1, activation='relu', kernel_initializer=GlorotUniform(seed=42), padding="same"))
model.add(MaxPooling1D(pool_size=(5), strides=(5), padding="same"))
model.add(Dropout(0.2))
model.add(Conv1D(filters=900, kernel_size=(4),  strides=1, activation='relu', kernel_initializer=GlorotUniform(seed=42), padding="same"))
model.add(MaxPooling1D(pool_size=(10), strides=(10), padding="same"))
model.add(Dropout(0.2))
model.add(Conv1D(filters=900, kernel_size=(8),  strides=1, activation='relu', kernel_initializer=GlorotUniform(seed=42), padding="same"))
model.add(MaxPooling1D(pool_size=(4), strides=(4), padding="same"))
model.add(Dropout(0.2))
model.add(Flatten())
model.add(Dense(output_shape, activation='softmax', kernel_regularizer=l2(0.001)))
model.summary()


# In[7]:


with tf.device("CPU"):    
    early_stop = EarlyStopping(monitor='val_loss', patience=10, verbose=1)
    wBestModel = '/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/weights_new_softmax/model.h5'
    model.load_weights(wBestModel)
    best_model = ModelCheckpoint(wBestModel, verbose=1, save_best_only=True)
    hBestModel     = '/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/histories_new_softmax/model.csv'
    history_logger = CSVLogger(hBestModel, separator=",", append=True)


# In[8]:


with tf.device("CPU"):
    model.compile(Adam(learning_rate=0.0001),
                      loss    = "binary_crossentropy", 
                      metrics = [tf.keras.metrics.AUC()])


# In[9]:


# with tf.device('CPU'):
#     model.evaluate(X_val,  y_val,  verbose=1)
#     model.evaluate(X_test, y_test, verbose=1)

#The model.evaluate function predicts the output for the given input and then computes the metrics function specified in the model.compile and based on y_true and y_pred and returns the computed metric value as the output.
#The keras.evaluate() function will give you the loss value for every batch. The keras.predict() function will give you the actual predictions for all samples in a batch, for all batches. So even if you use the same data, the differences will be there because the value of a loss function will be almost always different than the predicted values. These are two different things.


# In[10]:


with tf.device('CPU'):
    eval_val = model.evaluate(X_val,  y_val,  verbose=1)
    eval_test = model.evaluate(X_test, y_test, verbose=1)


# In[11]:


print(eval_val)
print(eval_test)

# >>> print(eval_val)
# [0.25856637954711914, 0.7514680027961731]
# >>> print(eval_test)
# [0.2572687268257141, 0.7519823908805847]

#auc: if the number of test samples is 1000 and model classifies 952 of those correctly, then the model's accuracy is 95.2%.
#If the model's prediction is perfect, the loss is zero; otherwise, the loss is greater.


# In[12]:


with tf.device('CPU'):
    pred_val  = model.predict(X_val, verbose=1)
    pred_test = model.predict(X_test, verbose=1)

#evaluate() will predict + compute metrics on your test set and trainer. predict() will only predict labels on your test set


# In[13]:


print(pred_val)
print(pred_test)
print(len(pred_val[0])) #so for each file, the prediction is given per however much input data was given
print(len(pred_val))
# >>> print(len(pred_val[0])) #so for each file, the prediction is given per however much input data was given
# 65
# >>> print(len(pred_val))
# 64256


# In[14]:


np.savetxt('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/prediction/pred_val.csv',  pred_val, delimiter='\t')
np.savetxt('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/prediction/true_val.csv',  y_val,    delimiter='\t')


# In[15]:


np.savetxt('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/prediction/pred_test.csv', pred_test, delimiter='\t')
np.savetxt('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/prediction/true_test.csv', y_test,    delimiter='\t')


# In[16]:


import random
colors =['black', 'k', 'dimgray', 'dimgrey', 'tab:gray', 'tab:grey', 'gray', 'grey', 'darkgray', 'darkgrey', 'silver', 'lightgray', 'lightgrey', 
         'gainsboro', 'whitesmoke', 'w', 'white', 'snow', 'rosybrown', 'lightcoral', 'indianred', 'brown', 'firebrick', 'maroon', 'darkred', 'r', 
         'red', 'mistyrose', 'salmon', 'tomato', 'tab:brown', 'darksalmon', 'coral', 'orangered', 'lightsalmon', 'sienna', 'seashell', 'chocolate', 
         'saddlebrown', 'sandybrown', 'tab:orange', 'peachpuff', 'peru', 'linen', 'bisque', 'darkorange', 'burlywood', 'antiquewhite', 'tan', 
         'navajowhite', 'blanchedalmond', 'papayawhip', 'moccasin', 'orange', 'wheat', 'oldlace', 'floralwhite', 'darkgoldenrod', 'goldenrod',
         'cornsilk', 'gold', 'lemonchiffon', 'khaki', 'palegoldenrod', 'darkkhaki', 'ivory', 'beige', 'lightyellow', 'lightgoldenrodyellow', 
         'olive', 'y', 'yellow', 'tab:olive', 'olivedrab', 'yellowgreen', 'darkolivegreen', 'greenyellow', 'chartreuse', 'lawngreen', 'honeydew', 
         'darkseagreen', 'palegreen', 'lightgreen', 'tab:green', 'forestgreen', 'limegreen', 'darkgreen', 'g', 'green', 'lime', 'seagreen', 
         'mediumseagreen', 'springgreen', 'mintcream', 'mediumspringgreen', 'mediumaquamarine', 'aquamarine', 'turquoise', 'lightseagreen', 
         'mediumturquoise', 'azure', 'lightcyan', 'paleturquoise', 'darkslategray', 'darkslategrey', 'teal', 'darkcyan', 'c', 'aqua', 'cyan', 
         'darkturquoise', 'cadetblue', 'tab:cyan', 'powderblue', 'lightblue', 'deepskyblue', 'skyblue', 'lightskyblue', 'tab:blue', 'steelblue', 
         'aliceblue', 'dodgerblue', 'lightslategray', 'lightslategrey', 'slategray', 'slategrey', 'lightsteelblue', 'cornflowerblue', 'royalblue', 
         'ghostwhite', 'lavender', 'midnightblue', 'navy', 'darkblue', 'mediumblue', 'b', 'blue', 'slateblue', 'darkslateblue', 'mediumslateblue', 
         'mediumpurple', 'rebeccapurple', 'blueviolet', 'tab:purple', 'indigo', 'darkorchid', 'darkviolet', 'mediumorchid', 'thistle', 'plum', 
         'violet', 'purple', 'darkmagenta', 'm', 'fuchsia', 'magenta', 'orchid', 'tab:pink', 'mediumvioletred', 'deeppink', 'hotpink', 'lavenderblush', 
         'palevioletred', 'crimson', 'pink', 'lightpink', 'tab:red']
print(len(colors))
colors.remove('deeppink')
colors.remove('navy')
for i in colors:
    if ("white" in i) | ("light" in i):
        colors.remove(i)
random.shuffle(colors)


# In[ ]:


pred_val = np.loadtxt("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/prediction/pred_val.csv", delimiter='\t')
y_val    = np.loadtxt("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/prediction/true_val.csv", delimiter='\t')


# In[ ]:


pred_test = np.loadtxt("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/prediction/pred_test.csv", delimiter='\t')
y_test    = np.loadtxt("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/prediction/true_test.csv", delimiter='\t')


# In[17]:


from sklearn.metrics import roc_auc_score
from sklearn.metrics import RocCurveDisplay
from itertools import cycle
import matplotlib.pyplot as plt
from sklearn.metrics import auc, roc_curve


# In[18]:


target_names = obj.var.index.tolist() #classifiers
n_classes = len(obj.var.index.tolist()) #number of classifiers ie files


# In[19]:


# store the fpr, tpr, and roc_auc for all averaging strategies
fpr, tpr, roc_auc = dict(), dict(), dict()
# Compute micro-average ROC curve and ROC area
fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), pred_test.ravel())
roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

print(f"Micro-averaged One-vs-Rest ROC AUC score:\n{roc_auc['micro']:.2f}")


# In[20]:


for i in range(n_classes):
    fpr[i], tpr[i], _ = roc_curve(y_test[:, i], pred_test[:, i])
    roc_auc[i] = auc(fpr[i], tpr[i])

fpr_grid = np.linspace(0.0, 1.0, 1000)

# Interpolate all ROC curves at these points
mean_tpr = np.zeros_like(fpr_grid)

for i in range(n_classes):
    mean_tpr += np.interp(fpr_grid, fpr[i], tpr[i])  # linear interpolation

# Average it and compute AUC
mean_tpr /= n_classes

fpr["macro"] = fpr_grid
tpr["macro"] = mean_tpr
roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

print(f"Macro-averaged One-vs-Rest ROC AUC score:\n{roc_auc['macro']:.2f}")


# In[21]:


print(fpr_grid)


# In[22]:


macro_roc_auc_ovr = roc_auc_score(
    y_test,
    pred_test,
    multi_class="ovr",
    average="macro"
)

print(f"Macro-averaged One-vs-Rest ROC AUC score:\n{macro_roc_auc_ovr:.2f}")


# In[23]:


fig, ax = plt.subplots(figsize=(10, 10))

# colors = cycle(["aqua", "darkorange", "cornflowerblue"])
for class_id, color in zip(range(n_classes), colors):
    RocCurveDisplay.from_predictions(
        y_test[:, class_id],
        pred_test[:, class_id],
        name=f"ROC curve for {target_names[class_id]}",
        color=color,
        ax=ax
    )

plt.plot(
    fpr["micro"],
    tpr["micro"],
    label=f"micro-average ROC curve (AUC = {roc_auc['micro']:.2f})",
    color="deeppink",
    linestyle=":",
    linewidth=4,
)

plt.plot(
    fpr["macro"],
    tpr["macro"],
    label=f"macro-average ROC curve (AUC = {roc_auc['macro']:.2f})",
    color="navy",
    linestyle=":",
    linewidth=4,
)

plt.axis("square")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Extension of Receiver Operating Characteristic\nto One-vs-Rest multiclass")
plt.legend(bbox_to_anchor =(1.8,0.0), loc='lower right')
#plt.axhline(y = 0.8, color = 'black', linestyle = '--') 
#plt.axvline(x = 0.2, color = 'black', linestyle = '--') 
plt.show()
fig.savefig('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/training_model_auc_roc_classifier_plots_testing_data_190624.png', bbox_inches='tight')

