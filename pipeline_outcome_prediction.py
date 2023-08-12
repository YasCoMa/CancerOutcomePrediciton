import os
import networkx as nx

from sklearn.base import clone
from sklearn.metrics import accuracy_score
from sklearn.model_selection import KFold
from sklearn.model_selection import ParameterGrid
from sklearn.model_selection import StratifiedKFold
from sklearn.svm import SVC

import graphkernels.kernels as gk # use this version> https://github.com/eghisu/GraphKernels
import igraph as ig
import numpy as np

class KernelGridSearchCV:
    """
    A simple class for performing a grid search for kernel matrices with
    a cross-validation strategy. At present, the class interface follows
    the default interface of `scikit-learn`. However, the class is *not*
    yet inheriting from any base class.

    The class only supports scoring based on accuracy so far.
    """

    def __init__(self, clf, param_grid, cv = None, random_state = None, refit = True):
        self.clf_             = clf
        self.grid_            = param_grid
        self.cv_              = cv
        self.random_state_    = random_state
        self.refit_           = refit
        self.best_estimator_  = None
        self.best_score_      = None

    def fit(self, X, y):

        # Use stratified k-folds with a user-specified number
        if self.cv_ is None:
            cv = KFold(
                    n_splits = 3,
                    shuffle = True,
                    random_state = self.random_state_
            )
        elif isinstance(self.cv_, int):
            cv = StratifiedKFold(
                    n_splits = self.cv_,
                    shuffle = True,
                    random_state = self.random_state_
            )
        else:
            cv = self.cv_

        grid = ParameterGrid(self.grid_)
        for parameters in grid:
            clf = self.clf_
            clf.set_params(**parameters)

            scores = []
            for train, test in cv.split(np.zeros(len(y)), y):
                X_train = X[train][:, train]
                #print(y)
                y_train = y[train]
                X_test  = X[test][:, train]
                y_test  = y[test]

                clf.fit(X_train, y_train)

                # The class only supports the accuracy score for now.
                ac = accuracy_score(y_test, clf.predict(X_test))
                scores.append(ac)

            score = np.mean(scores)
            if self.best_score_ is None or score > self.best_score_:
                self.best_estimator_ = clone(clf)
                self.best_score_     = score
                self.best_params_    = parameters

# Step 1
class Generate_personalised_networks:
    def load_interactome(self, folder, interactome_file):
        interactome=[]

        f=open(folder+interactome_file,"r")
        for line in f:
            l=line.replace("\n","").split("\t")
            interactome.append([l[0], l[1]])
        f.close()

        return interactome

    def load_expression_table(self, folder, expr_table_file):
        expr_data={}
        c=0
        f=open(folder+expr_table_file,"r")
        for line in f:
            l=line.replace("\n","").split("\t")
            if(c==0):
                for sample in l[1:]:
                    expr_data[sample]={}
            else:
                j=1
                for sample in expr_data.keys():
                    value=0.0
                    if(j<len(l)):
                        value=float(l[j])
                    expr_data[sample][l[0]]=value
                    j+=1
            c+=1
        f.close()

        return expr_data

    def make_patients_graph(self, folder, interactome_file, expr_table_file):
        print("Loading interactome...")
        positive_pairs=self.load_interactome(folder, interactome_file)
        
        print("Loading expression table...")
        expr_data=self.load_expression_table(folder, expr_table_file)

        print("Generating the patients graphs and exporting to graphml...")
        if(not os.path.isdir(folder+"patient_graphs")):
            os.system("mkdir "+folder+"patient_graphs")
        else:
            os.system("rm "+folder+"patient_graphs/*")

        for k in expr_data.keys():
            if(k!=""):
                personal_net=[]
                for p in positive_pairs:
                    p1=p[0]
                    p2=p[1]
                    if(p1 in expr_data[k].keys() and p2 in expr_data[k].keys()):
                        if( ( expr_data[k][p1]>=1 and expr_data[k][p2]>=1) or ( expr_data[k][p1]<=-1 and expr_data[k][p2]<=-1) ):
                            personal_net.append([p1, p2])
                
                if(len(personal_net)>0):
                    gr=nx.Graph()
                    for p in personal_net:
                        gr.add_edge(p1,p2)
                    nx.write_graphml(gr, folder+"patient_graphs/graph_"+k+".graphml")

        
# Step 2
class Import_patient_graphs:
    def load_dataset_networks(self, folder, labels_file):
        data_directory = folder+'patient_graphs'
        y_real={}
        with open(folder+labels_file) as f:
            for line in f:
                l=line.replace("\n","").split("\t")
                y_real[l[0]]=float(l[1])
        
        X = []
        y = []
        for filename in os.listdir(data_directory):
            full_path=os.path.join(data_directory, filename)
            
            sample=filename.replace("graph_","").replace(".graphml","")
            y.append(y_real[sample])
            
            if(os.path.isfile(full_path)):
                gr=ig.Graph.Read_GraphML(full_path)
                el=[]
                for j in range(gr.ecount()):
                    el.append(1)
                vl=[]
                for j in range(gr.vcount()):
                    vl.append(1)
                gr.vs["v_label"]=vl
                gr.es["e_label"]=el
                X.append(gr)
        
        y=np.array(y)
        
        return [X,y]

class Graph_kernel_classification:

    def execute_classification(self, folder, labels_file):
        print("Loading Patients graphs...")
        ipg = Import_patient_graphs()
        data = ipg.load_dataset_networks(folder, labels_file)

        X=data[0]
        y=data[1]

        kernels = {
            #'vertex_hist_gauss': gk.CalculateVertexHistGaussKernel,
            #'edge_hist_gauss': gk.CalculateEdgeHistGaussKernel,
            #'vertex_edge_hist_gauss': gk.CalculateVertexEdgeHistGaussKernel
        
            'vertex_histogram': gk.CalculateVertexHistKernel,
            'edge_histogram': gk.CalculateEdgeHistKernel,
            'vertex_edge_histogram': gk.CalculateVertexEdgeHistKernel,
            'vertex_vertex_edge_histogram': gk.CalculateVertexVertexEdgeHistKernel
            
            #'shortest_path': gk.CalculateShortestPathKernel, 
            #'graphlet': gk.CalculateGraphletKernel,
            #'connected_graphlet': gk.CalculateConnectedGraphletKernel,
            #'weisfeiler_lehman': gk.CalculateWLKernel,
            #'geometric_random_walk': gk.CalculateGeometricRandomWalkKernel, 
            #'exponential_random_walk': gk.CalculateExponentialRandomWalkKernel, 
            #'k_stop_random_walk': gk.CalculateKStepRandomWalkKernel
        }

        print("Calculating kernel matrix...")
        kernel_matrices = dict()
        for kernel_name, kernel_function in sorted(kernels.items()):
            print('\tKernel matrix for', kernel_name)
            kernel_matrices[kernel_name] = kernel_function(X)

        print("Getting accuracy results...")
        results={}
        for kernel_name, kernel_matrix in sorted(kernel_matrices.items()):
            print('\tExecuting Kernel ', kernel_name)
            grid = {
                    'C': 10. ** np.arange(-2,3)
            }

            clf = SVC(kernel='precomputed')

            grid_search = KernelGridSearchCV(
                clf,
                param_grid = grid,
                cv = 10, 
                random_state = 42,
            )

            grid_search.fit(kernel_matrix,y)

            clf = grid_search.best_estimator_
            results[kernel_name] = grid_search.best_score_ * 100.0

        return results

# Step 3
class Export_results:
    def make_report(self, folder, results):
        print("Exporting report...")
        f=open(folder+"outcome_prediction_report.tsv","w")
        for k in results.keys():
            f.write( "%s\t%.4f\n" %(k, results[k]) )
        f.close()

class Running:
    def run(self, folder, interactome_file, expr_table_file, labels_file):
        step1=Generate_personalised_networks()
        step1.make_patients_graph(folder, interactome_file, expr_table_file)

        step2=Graph_kernel_classification()
        results=step2.execute_classification(folder, labels_file)

        step3=Export_results()
        step3.make_report(folder, results)

class Running_config:
    def run(self, args):
        if(args.folder!="" and args.interactome_file!="" and args.expression_table_file!="" and args.label_file!=""):
            if(os.path.isdir(args.folder) and os.path.isfile(args.folder+args.interactome_file) and os.path.isfile(args.folder+args.expression_table_file) and os.path.isfile(args.folder+args.label_file) ):
                r=Running()
                r.run(args.folder, args.interactome_file, args.expression_table_file, args.label_file)
            else:
                print("Error: There are invalid folder or files, some of them were not found")
        else:
            print("Error: There are missing parameters")

if __name__ == '__main__':
    import argparse
    from argparse import RawTextHelpFormatter
    parser = argparse.ArgumentParser(description='Pipeline for outcome disease prediction', formatter_class=RawTextHelpFormatter)
    parser.add_argument("-fo", "--folder", action="store", help="Required - Folder to store the files (use the folder where the other files required can be found, ex.: /home/user/experiment/ )")
    parser.add_argument("-if", "--interactome_file", action="store", help="Required - File with the pairs (two columns with uniprot identifiers in tsv format)")
    parser.add_argument("-etf", "--expression_table_file", action="store", help="Required - File with the expression values for the genes by sample/patient in tsv format")
    parser.add_argument("-lf", "--label_file", action="store", help="Required - File with the prognosis label for each sample")
    args = parser.parse_args()
    r=Running_config()
    r.run(args)
