# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 21:30:31 2018

@authors: Barbiero Pietro and Ciravegna Gabriele
"""


# sklearn library
from sklearn import datasets
from sklearn import decomposition
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score
from sklearn.neural_network import MLPClassifier
from sklearn.linear_model import RidgeClassifierCV

from keras.models import Sequential
from keras.layers import Dense, Activation
from keras.utils import to_categorical

from anytree import NodeMixin, RenderTree, LevelOrderIter, search
import networkx as nx

import sys
import os
import datetime
import numpy as np
import matplotlib.pyplot as plt
import copy
from scipy import stats
from scipy.spatial.distance import pdist, euclidean
import csv
import time

def retrieve_n_class_color_cubic(N):
    '''
    retrive color code for N given classes
    Input: class number
    Output: list of RGB color code
    '''

    # manualy encode the top 8 colors
    # the order is intuitive to be used
    color_list = [
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1),
        (1, 1, 0),
        (0, 1, 1),
        (1, 0, 1),
        (0, 0, 0),
        (1, 1, 1)
    ]

    # if N is larger than 8 iteratively generate more random colors
    np.random.seed(1)  # pre-define the seed for consistency

    interval = 0.5
    while len(color_list) < N:
        the_list = []
        iterator = np.arange(0, 1.0001, interval)
        for i in iterator:
            for j in iterator:
                for k in iterator:
                    if (i, j, k) not in color_list:
                        the_list.append((i, j, k))
        the_list = list(set(the_list))
        np.random.shuffle(the_list)
        color_list.extend(the_list)
        interval = interval / 2.0

    return color_list[:N]

options = {
    'node_color': 'black',
    'node_size': 100,
    'width': 3,
}

class GhcoreNode(NodeMixin):
	
	def __init__(self, node_id, child_id=None, parent=None, X=[], y=[], sample_indeces=[], child_indeces=[], w=None, T=None, etime=1, graph=nx.Graph()):
		super().__init__()
		self.node_id = node_id
		self.child_id = child_id
		self.parent = parent
		self.X = X
		self.y = y
		self.sample_indeces = sample_indeces
		self.child_indeces = child_indeces
		self.w = w
		self.T = T
		self.etime = etime
		self.graph = graph
		
	def display(self):
		for pre, _, node in RenderTree(self):
			treestr = u"%s%s" % (pre, node.node_id)
			purity = 0 
			if len(node.y) > 0: 
				purity = np.max( np.unique(node.y, return_counts=True) ) / len(node.y)
			print(treestr.ljust(8), 'n_samples: %d; purity = %.2f' %(len(node.y), purity), end='')
			print('')
	
	def set_w(self, w):
		self.w = w
		
	def set_T(self, T):
		self.T = T
		
	def set_up_child_indeces(self):
		self.child_indeces = [-1 for i in range(0, len(self.y))]
		
	def increment_elapsed_time(self):
		self.etime = self.etime + 1
	
	def update_child_id(self, child_node):
		child_node.child_id = len(self.children)-1
	
	def add_sample(self, sample, target, j, child):
		
		old_child_idx = self.child_indeces[j]
		
#		print(child.child_id)
#		print(child.sample_indeces)
#		print(self.child_indeces)
#		print(j)
		
		if old_child_idx != child.child_id:
			
			if old_child_idx != -1:
				
				old_child = self.children[old_child_idx]
				sample_j = old_child.sample_indeces.index(j)
				
#				print(old_child.X.shape)
				old_child.X = np.delete(old_child.X, sample_j, axis=0)
				old_child.y = np.delete(old_child.y, sample_j, axis=0)
				old_child.sample_indeces.remove(j)
#				print(old_child.X.shape)
			
			if len(child.X) == 0:
				child.X, child.y = sample, target
			else:
				child.X = np.concatenate((child.X, sample))
				child.y = np.concatenate((child.y, target))
			
			child.sample_indeces.append(j)
			self.child_indeces[j] = child.child_id
			
#			print(child.sample_indeces)
#			print(self.child_indeces)
#			print()
		
			
	def soft_competitive_learning(self, epsilon, sample):
		Delta = epsilon * (sample - self.w) / self.etime
		self.w = self.w + Delta
		
#	def soft_competitive_learning_2(self, epsilon_w, epsilon_n, sample, winner_1):
#		Delta = epsilon_w * (sample - self.w) / self.etime
#		self.w = self.w + Delta
		
	def update_threshold(self, child_node):
		neighbours_id = self.get_graph_neighbours(child_node)
		neighbours = search.findall(self, filter_=lambda node: node.node_id in neighbours_id, maxlevel=2)
		neighbours_W = np.array([node.w.squeeze() for node in neighbours])
		
		distances = np.sum( (child_node.w - neighbours_W)**2 , 1)
		
		if len(distances) > 1: average_distance = np.mean(distances)
		else: average_distance = distances[0]
		
		max_distance = np.max((self.T, average_distance))
		
		if max_distance == None: 
			print(1)
			
		child_node.set_T(max_distance)
		
		return neighbours
	
	
		
		
	def add_graph_node(self, node_id):
		self.graph.add_node(node_id)
		
	def add_graph_edge(self, winner_1, winner_2):
		self.graph.add_edge(winner_1.node_id, winner_2.node_id)
	
	def get_graph_neighbours(self, node):
		return list(self.graph.adj[node.node_id])
	
	def draw_graph(self):
		plt.figure()
		nx.draw(self.graph, with_labels=True, font_weight='bold')
		plt.show()
		
		
		
	def plot_local_quantization(self, accuracy, n_leaves):
		
		nclass = len(np.unique(self.root.y))
		colors = np.array(retrieve_n_class_color_cubic(N=nclass))
		cy = np.array([colors[i].squeeze() for i in self.root.y-1])
		
		W = np.array([child.w.squeeze() for child in self.children])
		
		plt.figure()
		plt.scatter(self.root.X[:, 0], self.root.X[:, 1], c=cy, marker='.', alpha=0.3, label='voronoi set')
		plt.scatter(W[:, 0], W[:, 1], c='k', marker='o', label='gexin')
		plt.title('Ghcore - h=%d - #L=%d - acc.=%.2f' %(self.root.height, n_leaves, accuracy))
		plt.legend()
		plt.show()
		
	def plot_quantization(self, accuracy, leaves):
		
		nclass = len(np.unique(self.root.y))
		colors = np.array(retrieve_n_class_color_cubic(N=nclass))
		cy = np.array([colors[i].squeeze() for i in self.root.y-1])
		
		X_arch_core = np.array([leaf.w.squeeze() for leaf in leaves])
		y_arch_core = []
		for leaf in leaves: 
			if len(leaf.y) > 0:
				unique_y, count_y = np.unique(leaf.y, return_counts=True)
				amax = np.argmax(count_y)
				y_arch_core.append(unique_y[amax])
		y_arch_core = np.array(y_arch_core)
				
		ccore = np.array([colors[i].squeeze() for i in y_arch_core-1])
		
		plt.figure()
		plt.scatter(self.root.X[:, 0], self.root.X[:, 1], c=cy, marker='.', alpha=0.2, label='voronoi set')
		plt.scatter(X_arch_core[:, 0], X_arch_core[:, 1], c=ccore, marker='o', label='gexin')
		plt.title('Ghcore - h=%d - #L=%d - acc.=%.2f' %(self.root.height, len(leaves), accuracy))
		plt.legend()
		plt.show()
		
def predict_by_core(root, X_test, y_test):
	
	leaves = [node for node in LevelOrderIter(root) if node.is_leaf and len(node.y) > 0]
	
	X_arch_core = np.array([leaf.w.squeeze() for leaf in leaves if len(leaf.y) > 0])
	y_arch_core = []
	for leaf in leaves: 
		if len(leaf.y) > 0:
			unique_y, count_y = np.unique(leaf.y, return_counts=True)
			amax = np.argmax(count_y)
			y_arch_core.append(unique_y[amax])
	y_arch_core = np.array(y_arch_core)
		
	model = RidgeClassifierCV()
	model.fit(X_arch_core, y_arch_core)
	accuracy = model.score(X_test, y_test)
	
#	print ('# samples: %d' %(len(y_arch_core)))
#	
#	model = Sequential()
#	model.add(Dense(3, input_dim=X_arch_core.shape[1]))
#	model.add(Activation('relu'))
#	model.add(Dense(3))
#	model.add(Activation('relu'))
#	
#	# For a multi-class classification problem
#	model.compile(optimizer='rmsprop',
#	              loss='categorical_crossentropy',
#	              metrics=['accuracy'])
#	
#	# Convert labels to categorical one-hot encoding
#	one_hot_y = to_categorical(y_arch_core, num_classes=len(np.unique(y_arch_core)))
#	# Train the model, iterating on the data in batches of 32 samples
#	model.fit(X_arch_core, one_hot_y, epochs=10, batch_size=32)
#	
#	y_pred = model.predict(X_test, batch_size=100)
#	y_pred1D = y_pred.argmax(1)
#	
#	print ('Accuracy on validation data: %.2f' %(accuracy))
	
	return accuracy, leaves

def Ghcore(X, y, X_test, y_test, max_height, min_epochs, min_purity, epsilon_w, epsilon_n, min_size, min_accuracy):
	
	y = np.reshape(y, (len(y), 1))
	centroid_X = np.mean(X, axis=0)
	centroid_X = np.reshape(centroid_X, (1, len(centroid_X)))
	
	n_nodes = 0
	root = GhcoreNode('Node_' + str(n_nodes), parent=None, X=X, y=y, sample_indeces=[], w=centroid_X, T=np.Inf)
	root.set_up_child_indeces()
	n_nodes = n_nodes + 1
	
	k = 1
	parent = root
	accuracy = 0
	
	while k < max_height and accuracy < min_accuracy:
#		print("Vertical growth - height = %d" %(k))
		
		leaves = [node for node in LevelOrderIter(root) if node.is_leaf and len(node.y) > min_size]
		
		n_leaves = len(leaves)
		if n_leaves == 0:
			break;
		
		for i in range(0, n_leaves):
			
			parent = leaves[i]
			
			counter = 0
			epoch = 0
			purity = 0
			
			noise = np.random.uniform(0, 0.0001, parent.w.shape)
			n = GhcoreNode('Node_' + str(n_nodes), parent=parent, X=[], y=[], sample_indeces=[], w=parent.w+noise)
			parent.update_child_id(n)
			parent.add_graph_node('Node_' + str(n_nodes))
			n_nodes = n_nodes + 1
			n = GhcoreNode('Node_' + str(n_nodes), parent=parent, X=[], y=[], sample_indeces=[], w=parent.w-noise)
			parent.update_child_id(n)
			parent.add_graph_node('Node_' + str(n_nodes))
			n_nodes = n_nodes + 1
			
			while epoch < min_epochs or purity < min_purity:
				
				first_time = True
				
				# learning process
				for j in range(0, len(parent.y)):
						
#					if k > 2 and epoch > 0 and j > 3:
#						print(epoch)
#						print(j)
					
					sample = parent.X[j, :]
					sample = np.reshape(sample, (1, len(sample)))
					target = parent.y[j]
					target = np.reshape(target, (1, len(target)))
					
					W = np.array([leaf.w.squeeze() for leaf in parent.children])
					distances = np.sum( (sample - W)**2 , 1)
					winner_1_idx = np.argmin(distances)
					distances[winner_1_idx] = np.inf
					winner_2_idx = np.argmin(distances)
					
					winner_1 = parent.children[winner_1_idx]
					winner_2 = parent.children[winner_2_idx]
					
					if first_time:
						first_time = False
						avgT = np.mean( pdist(parent.X) )
						
						if epoch == 0:
							winner_1.set_T(avgT)
							winner_2.set_T(avgT)
							parent.set_T(avgT)
							
						parent.add_sample(sample, target, j, winner_1)
						winner_1.increment_elapsed_time()
						winner_1.soft_competitive_learning(epsilon_w, sample)
						
						parent.add_graph_edge(winner_1, winner_2)
						
#						parent.draw_graph()
						
					else:
						
						if False: #parent.get_graph_neighbours(winner_1) >= parent.X.shape[1]:
							
							# use convex hull
							1
							
						else:
							if winner_1.T == None: 
								print(1)
							explainable = euclidean( winner_1.w, sample ) < winner_1.T
							
						if explainable:
							
							parent.add_sample(sample, target, j, winner_1)
							winner_1.increment_elapsed_time()
							winner_1.soft_competitive_learning(epsilon_w, sample)
							
							parent.add_graph_edge(winner_1, winner_2)
							
#							parent.draw_graph()
							
							neighbours = parent.update_threshold(winner_1)
							
							for neighbour in neighbours:
								neighbour.soft_competitive_learning(epsilon_n, sample)
								parent.update_threshold(neighbour)
								
						else:
							
							new_node = GhcoreNode('Node_' + str(n_nodes), parent=parent, X=[], y=[], sample_indeces=[], w=sample)
							parent.update_child_id(new_node)
							parent.add_graph_node('Node_' + str(n_nodes))
							n_nodes = n_nodes + 1
							
							parent.add_sample(sample, target, j, new_node)
							
							new_node.set_T(parent.T)
							counter = 0
							
							if new_node.T == None: 
								print(1)
				
				purities = [np.max( np.unique(node.y, return_counts=True) ) / len(node.y) for node in parent.children if len(node.y) > 0]
				purity = np.mean(purities)
				
				epoch = epoch + 1
				counter = counter + 1
			
			for child in parent.children:
				child.set_up_child_indeces()
#			parent.draw_graph()
		accuracy, leaves = predict_by_core(root, X_test, y_test)
		parent.plot_quantization(accuracy, leaves)
			
#		accuracy = predict_by_core(root, X_test, y_test)
				
		k = k + 1

	root.display()
	leaves = [node for node in LevelOrderIter(root) if node.is_leaf and len(node.y) > 0]
	
	return root, leaves



print("Loading datasets...")
X, y = datasets.load_iris(return_X_y=True)
X = X[:, 2:4]

skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)
list_of_splits = [split for split in skf.split(X, y)]
train_index, test_index = list_of_splits[0]
X_train, y_train = X[train_index], y[train_index]
X_test, y_test = X[test_index], y[test_index]

nn = Ghcore(X_train, y_train, X_test, y_test, max_height=6, min_epochs=3, min_purity=0.6, epsilon_w=0.2, epsilon_n=0.01, min_size=5, min_accuracy=0.8)

