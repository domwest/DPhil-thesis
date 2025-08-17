#IN TERMINAL!!!

#go to /well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/run_training_gpu
#run this via command:
#sbatch -p gpu_short --gres gpu:1 submit_04_training_310524.sh

import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns
import h5py
from matplotlib import pyplot
from matplotlib import colormaps
from sklearn import preprocessing
from sklearn.model_selection import train_test_split

from tensorflow.keras.layers import Dense, Flatten, Dropout, Input, Conv1D, MaxPooling1D, BatchNormalization
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint, CSVLogger
from tensorflow.keras.initializers import GlorotUniform
from tensorflow.keras.models import Sequential
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.regularizers import l2
from tensorflow.keras import backend

from keras.utils import plot_model

import tensorflow as tf
import keras

import pydot
import graphviz

gpu_device = 'GPU:0'

with tf.device(gpu_device):
    labels = h5py.File('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/ML_training_sets/labels.h5', 'r')
    labels = labels['label'][:]
    labels.shape
# (321280, 65)

with tf.device(gpu_device):
    data = h5py.File('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/ML_training_sets/sequences.h5', 'r')
    data = data['data'][:]
    data.shape
# (321280, 1000, 4)

data
data.shape
#so this is the actual hot encode format of the sequences (1000bp seqs) per peak

with tf.device(gpu_device):
    X_train, X_test, y_train, y_test = train_test_split(data,   labels, test_size=0.4, random_state=42, shuffle=True)
    X_val,   X_test, y_val,   y_test = train_test_split(X_test, y_test, test_size=0.5, random_state=42, shuffle=True)
    print("Train:     ", X_train.shape, y_train.shape)
    print("Validation:", X_val.shape,   y_val.shape)
    print("Test:      ", X_test.shape,  y_test.shape)
# Train:      (192768, 1000, 4) (192768, 65)
# Validation: (64256, 1000, 4) (64256, 65)
# Test:       (64256, 1000, 4) (64256, 65)

#Random state: Controls the shuffling applied to the data before applying the split.
#why 60% of data for training and 40% for validation/testing? Why not 80% for training and 20% for testing, for example?
#So in the first row, testing data becomes 40% * 38882 = 15553 & training data becomes 60% * 38882 = 23329
#And then in the second row, your X_val and X_test each become half of the testing data
#End up with 60% training data, 20% validation data and 20% testing data

#1). X_train - This includes your all independent variables,these will be used to train the model, also as we have specified the test_size = 0.4, this means 60% of observations from your complete data will be used to train/fit the model and rest 40% will be used to test the model.
#2). X_test - This is remaining 40% portion of the independent variables from the data which will not be used in the training phase and will be used to make predictions to test the accuracy of the model.
#3). y_train - This is your dependent variable which needs to be predicted by this model, this includes category labels against your independent variables, we need to specify our dependent variable while training/fitting the model.
#4). y_test - This data has category labels for your test data, these labels will be used to test the accuracy between actual and predicted categories.

y_test
y_train 
np.sum(y_train, axis=0)
#so this gives you a count of how many peaks are present in the cardiac_muscle_cell and then the bulk tissue heart LV

with tf.device(gpu_device):    
    obj = sc.read_h5ad("/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/cm.h5ad")
    df_train = pd.DataFrame([np.sum(y_train, axis=0)], index=['count'], columns=obj.var.index.tolist()).T
    df_train['set'] = 'train_label'
    df_val   = pd.DataFrame([np.sum(y_val, axis=0)],   index=['count'], columns=obj.var.index.tolist()).T
    df_val['set'] = 'val_label'
    df_test  = pd.DataFrame([np.sum(y_test, axis=0)],  index=['count'], columns=obj.var.index.tolist()).T
    df_test['set'] = 'test_label'
    df = pd.concat([df_train, df_val, df_test])
    df.sort_values(by=['count'], ascending=False, inplace=True)
    df['cell_type'] = df.index.tolist()
    df['logarithm_base10_count'] = np.log10(df['count'])
    fig, ax = pyplot.subplots(figsize=(34,18)) #was 16, 9
    #color = (0.341, 0.678, 0.58) #https://rgbcolorpicker.com/0-1
    ax = sns.barplot(data=df, x="logarithm_base10_count", y="cell_type", hue="set", ax=ax, palette="deep")#color=color)
    fig.savefig('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/training_figure_190624.png')
    #So if all of these counts of whether or not a specific peak is present in the cell_types has been divided randomly into the training, validation and testing sets then what would we be expecting to see here?
    #I suppose the training set has a greater chance of having more 1s present in either of the cell_type so that is bound to be bigger than validation and testing
    #You'd also expect the bulk tissue train, val, and test sets to have more 1s because most of the data comes from bulk tissue

df['cell_type'].value_counts()
#should be 3 each for training, validation, and testing
# Pericyte_General_2_ATAC                             3
# Endothelial_Myocardial_ATAC                         3
# Endothelial_ATAC                                    3
# cardiac_muscle_cell_CHIP_HISTONE_H3K27ac_single     3
# cardiac_muscle_cell_CHIP_HISTONE_H3K4me3_single     3
#                                                    ..
# atrial_cm_ATAC                                      3
# Cardiac_Fibroblast_ATAC                             3
# Mast.1_ATAC                                         3
# Ery_Don003_hg38_ATAC                                3
# bulk_heart_LV_tissue_CHIP_HISTONE_H3K4me3_single    3

with tf.device(gpu_device):
    batch_size = 100
    epochs     = 100

with tf.device(gpu_device):
    input_shape  = (X_train.shape[-2], X_train.shape[-1])
    output_shape = labels.shape[-1]

with tf.device(gpu_device):
    training_seqs   = tf.data.Dataset.from_tensor_slices((X_train, y_train)).batch(batch_size)
    validation_seqs = tf.data.Dataset.from_tensor_slices((X_val, y_val)).batch(batch_size)

X_train.shape[-2] #1000
X_train.shape[-1] #4
#so input shape is the 1000 bp hot encode sequences which are each 4 factor based 
labels.shape[-1] #65
#so output shape is whether there is a peak or not per celltype

#Context:

#Computational graphs are core components of neural networks 
#Computational graphs are made up of nodes and edges where the nodes are variables and the edges are data flow
#Neural networks are algorithms designed based on the human brain which contains many layers.Each layer contains a node called neuron which performs a maths operation
#Input layer takes in the data and preprocesses it
#Hidden layer performs nonlinear transformation of input
#Output layer takes the results from hidden layer, transforms them, and gives final output
#Sequential models are linear stacks of layers where one layer leads to the next. So each previous layer is the input to the next layer
#Can also say that a sequential model is appropriate for a plain stack of layers where each layer has exactly one input tensor and one output tensor
#So you basically create layers and then test layers on a test input
#Sequential models are used for simple classifier or declassifier models
#A classifier, or classification model, predicts categorical labels (classes). 

#Here we create a Sequential model incrementally via the add()
#When you instantiate a Sequential model without an input shape, it isn't "built": it has no weights (and calling model.weights results in an error stating just this). The weights are created when the model first sees some input data
#Model weights are all the parameters (including trainable and non-trainable) of the model which are in turn all the parameters used in the layers of the model
#When building the Sequential model incrementally above, we pass an Input object to the model so that it knows its input shape from the start and can generate model.summary() etc. from the get-go
#Models built with a predefined input shape like this always have weights (even before seeing any data) and always have a defined output shape.

model = Sequential()
model.add(Input(shape=input_shape))
model.add(Conv1D(filters=4, kernel_size=(1), strides=(1), activation='relu', kernel_initializer=GlorotUniform(seed=42), padding="same"))
#Convolutional neural network in 1D; creates a nr of small windows which float over the picture and each one is its own neural network which categorises based on this window... filter = 4 means how many of these little windows are there. kernel size = 1 so looking at one data point (?). activation='relu' is kin dof your default which gives you a yes or no but not a full yes or no... most common activation. Kernel_initializers play a crucial role in training Deep Learning models. They define the distribution of weights and biases, which can significantly impact a model's accuracy and speed of convergence.
#Regarding activation='relu':
#Activation functions give out the final value given out from a neuron, but what is activation function and why do we need it? 
#So, an activation function is basically just a simple function that transforms its inputs into outputs that have a certain range. There are various types of activation functions that perform this task in a different manner, For example, the sigmoid activation function takes input and maps the resulting values in between 0 to 1.
#One of the reasons that this function is added into an artificial neural network in order to help the network learn complex patterns in the data. These functions introduce nonlinear real-world properties to artificial neural networks. Basically, in a simple neural network, x is defined as inputs, w weights, and we pass f (x) that is the value passed to the output of the network. This will then be the final output or the input of another layer.
#If the activation function is not applied, the output signal becomes a simple linear function. A neural network without activation function will act as a linear regression with limited learning power. But we also want our neural network to learn non-linear states as we give it complex real-world information such as image, video, text, and sound.
#The function returns 0 if it receives any negative input, but for any positive value x, it returns that value back. Thus it gives an output that has a range from 0 to infinity.
#def ReLU(x):
#  if x>0:
#	  return x
#  else: 
#    return 0
model.add(Dropout(0.2)) 
#The Dropout layer randomly sets input units to 0 with a frequency of rate at each step during training time, which helps prevent overfitting.
#Overfitting refers to a model that models the training data too well. Overfitting happens when a model learns the detail and noise in the training data to the extent that it negatively impacts the performance of the model on new data.
model.add(Conv1D(filters=300, kernel_size=(20), strides=(1), activation='relu', kernel_initializer=GlorotUniform(seed=42), padding="same"))
model.add(MaxPooling1D(pool_size=(3), strides=(3), padding="same")) 
#After your Conv1D step you pool this into a neural network
model.add(Dropout(0.2))
model.add(Conv1D(filters=600, kernel_size=(10), strides=1, activation='relu', kernel_initializer=GlorotUniform(seed=42), padding="same"))
#A lot of repeating layers but with different options specified, so each step looks at the data differently, so you get added filtering here
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
#Flatten is taking this multidimensional data and putting it all into a 1D array
model.add(Dense(output_shape, activation='sigmoid', kernel_regularizer=l2(0.001)))
#In any neural network, a dense layer is a layer that is deeply connected with its preceding layer which means the neurons of the layer are connected to every neuron of its preceding layer.
#As known, the main difference between the Convolutional layer and the Dense layer is that Convolutional Layer uses fewer parameters by forcing input values to share the parameters. The Dense Layer uses a linear operation meaning every output is formed by the function based on every input.
#ReLu vs sigmoid activation:
#Context of what their graphs look like: sigmoid is an S shape & relu is flat and then linear, so looks like a hockey stick
#Relu has a faster training time than sigmoid and an array of other pros: https://wandb.ai/ayush-thakur/dl-question-bank/reports/ReLU-vs-Sigmoid-Function-in-Deep-Neural-Networks--VmlldzoyMDk0MzI
#Maybe sigmoid was used at the final step because: The main reason why we use sigmoid function is because it exists between (0 to 1). Therefore, it is especially used for models where we have to predict the probability as an output.Since probability of anything exists only between the range of 0 and 1, sigmoid is the right choice.
model.summary()

#plot_model(model, to_file='/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/model_plot_190624.png', show_shapes=True, show_layer_names=True)

with tf.device(gpu_device):    
    early_stop = EarlyStopping(monitor='val_loss', patience=10, verbose=1)
    #
    wBestModel = '/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/weights_new/model.h5'
    #model.load_weights(wBestModel)
    #
    best_model = ModelCheckpoint(wBestModel, verbose=1, save_best_only=True)
    #
    hBestModel     = '/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/histories_new/model.csv'
    history_logger = CSVLogger(hBestModel, separator=",", append=True)
    
#Assuming the goal of a training is to minimize the loss. With this, the metric to be monitored would be 'loss', and mode would be 'min'. A model.fit() training loop will check at end of every epoch whether the loss is no longer decreasing, considering the min_delta and patience if applicable. Once it's found no longer decreasing, model.stop_training is marked True and the training terminates.
#patience: Number of epochs with no improvement after which training will be stopped.
#verbose=1; Verbosity mode, 0 or 1. Mode 0 is silent, and mode 1 displays messages when the callback takes an action.

# In[63]:


with tf.device(gpu_device):
    model.compile(Adam(learning_rate=0.0001),
                      loss    = "binary_crossentropy", 
                      metrics = [tf.keras.metrics.AUC()])

#Model compilation is an activity performed after writing the statements in a model and before training starts. It checks for format errors, and defines the loss function, the optimizer or learning rate, and the metrics. A compiled model is needed for training but not necessary for predicting.
#Optimizer that implements the Adam algorithm.
#Adam optimization is a stochastic gradient descent method that is based on adaptive estimation of first-order and second-order moments.
#According to Kingma et al., 2014, the method is "computationally efficient, has little memory requirement, invariant to diagonal rescaling of gradients, and is well suited for problems that are large in terms of data/parameters".
#learning_rate: Defaults to 0.001. An optimal learning rate value (default value 0.001) means that the optimizer would update the parameters just right to reach the local minima. Varying learning rate between 0.0001 and 0.01 is considered optimal in most of the cases. You can start with a small value, such as 0.01, and increase or decrease it by a factor of 10 until you find a value that works well for your network.
#loss = binary_crossentropy: The model gives the output, how can we evaluate how good (or bad) is the prediction. This is the whole purpose of the evaluation metrics. 
#^: The loss function tells how good your model is in predictions. If the model predictions are closer to the actual values the Loss will be minimum and if the predictions are totally away from the original values the loss value will be the maximum.
#^:Loss= abs(Y_pred – Y_actual)
#^: binary_crossentropy: It is designed to measure the dissimilarity between the predicted probability distribution and the true binary labels of a dataset. Binary cross entropy compares each of the predicted probabilities to actual class output which can be either 0 or 1. It then calculates the score that penalizes the probabilities based on the distance from the expected value. That means how close or far from the actual value.
#metrics = tf.keras.metrics.AUC(): This metric creates four local variables, true_positives , true_negatives , false_positives and false_negatives that are used to compute the AUC

with tf.device('/GPU:0'):    
    history = model.fit(training_seqs,
                        batch_size=batch_size,
                        epochs=epochs,
                        verbose=1,
                        validation_data=(validation_seqs),
                        callbacks=[early_stop, best_model, history_logger])
    
#model.fit() trains the model for a fixed number of epochs (dataset iterations) ie in this case 100 iterations
#why is it doing 234 things within each epoch? I thought it would be 100 things since batch size is 100
#model.fit(
#train_data,
#steps_per_epoch = train_samples//batch_size #this is total nr training samples (X_train.shape gives us this) which is 23329 nr peaks divided by batch_size which is 100 = 234
#epochs = epochs, #100
#validation_data = validation_data,
#verbose = 1, #verbose: 'auto', 0, 1, or 2. Verbosity mode. 0 = silent, 1 = progress bar, 2 = one line per epoch.
#validation_steps = validation_samples//batch_size) #7776/100 = 78
#callbacks #A callback is an object that can perform actions at various stages of training (e.g. at the start or end of an epoch, before or after a single batch, etc).

#AUC ie area under the curve tells you how the model is performing
#ie from top left corner to bottom right corner, the model goes from better to worse
#True positive rate is on the y axis and false positive rate on the x axis

#loss curves help us evaluate the model's performance during training and how the model learns
#The curve has loss on the y axis and epochs on the x axis
#Learning curve of a good fit model has a moderately high training loss at the beginning which gradually decreases upon adding training examples and flattens gradually, indicating addition of more training examples doesn't improve the model performance on training data.

history = pd.read_csv('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/histories_new/model.csv', sep=',')
#history.rename({0:'epoch', 1:'auc', 2:'loss', 3:'val_auc', 4:'val_loss'}, axis=1, inplace=True)
#['auc', 'loss', 'val_auc', 'val_loss']
#history = history.tail(55)
history['epoch'] = history['epoch'] + 1
hist_df = history

# history
# print(history.history.keys())
# #history.history['auc']
# hist_df = pd.DataFrame(history.history) 
# hist_df
#print(epoch)
#hist_df['epoch'] = list(range(1, 56)) #because early stopping epoch ie 55 plus 1 
print(hist_df)

auc1 = hist_df[['epoch', 'auc']]
auc1['type'] = 'train'
auc1.columns = ['epoch', 'auc', 'type']
auc2 = hist_df[['epoch', 'val_auc']]
auc2['type'] = 'val'
auc2.columns = ['epoch', 'auc', 'type']

loss1 = hist_df[['epoch', 'loss']]
loss1['type'] = 'train'
loss1.columns = ['epoch', 'loss', 'type']
loss2 = hist_df[['epoch', 'val_loss']]
loss2['type'] = 'val'
loss2.columns = ['epoch', 'loss', 'type']

auc = pd.concat([auc1, auc2])
loss = pd.concat([loss1, loss2])

auc_max   = np.max(auc[auc['type']=='val']['auc'])
loss_min = np.min(loss[loss['type']=='val']['loss'])

tmp = loss[loss['type']=='val']['loss'].tolist()
opt_epoch = tmp.index(min(tmp))

import matplotlib.pyplot as plt

f, axs = plt.subplots(1,2,figsize=(16,5))
axs[0] = sns.lineplot(data=auc, x='epoch', y='auc', hue='type', ax=axs[0], palette="Accent")
# axs[0] = sns.scatterplot(x=auc_epoch, y=auc_max, ax=axs[0])
axs[0].axhline(y = auc_max, xmin = 0, xmax = 1, color = "black", linestyle = "dashed")
axs[0].axvline(x = opt_epoch, ymin = 0, ymax = 1, color = "black", linestyle = "dashed")
axs[1] = sns.lineplot(data=loss, x='epoch', y='loss', hue='type', ax=axs[1], palette="Accent")
axs[1].axhline(y = loss_min, xmin = 0, xmax = 1, color = "black", linestyle = "dashed")
axs[1].axvline(x = opt_epoch, ymin = 0, ymax = 1, color = "black", linestyle = "dashed")
axs[0].legend().set_title('')
axs[1].legend().set_title('')

plt.tight_layout()
plt.show(block=False)
f.savefig('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/training_model_auc_loss_plots_190624.png')
plt.close("all")

#Comparing performance on training data vs out-of-sample data can give an idea of whether performance could be improved via bias reduction or variance reduction. If both have poor performance, you would suspect that you haven't captured the trends in the data; you need more parameters or perhaps additional variables to explain the outcome. (This is high bias.) If you have strong in-sample performance but weak out-of-sample performance, you would suspect that you have overfit to the training data. (This is high variance.)

#model.save('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/model')
#del model
#model = load_model('/well/PROCARDIS/domwest/deephaem_prep/input_files_for_deephaem/model')

