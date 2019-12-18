import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import collections
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

def get_strains_clusters(linear_params, pca_data):
    samples_count = len(pca_data)
    strains_id = []
    for i in range(samples_count):
        distance_from_real_y = []
        for j in range(3):
            distance_from_real_y.append(
                abs(pca_data[i][0] * linear_params[j][0] + linear_params[j][1] - pca_data[i][1]))
        strains_id.append(float(distance_from_real_y.index(min(distance_from_real_y))))
    return strains_id

def save_fig_strains_clusters(pca_data, strains_id):
    colors = ['b','r','g']
    for i in range(3):
        v = pca_data[np.where(np.array(strains_id) == i)]
        plt.scatter(v[:, 0], v[:, 1], c=colors[i], s=10, alpha=0.5, label='strain {}'.format(i))
    plt.legend(loc=2)
    plt.title('2 component PCA of strains')
    plt.xlabel('pc1', fontsize=12)
    plt.ylabel('pc2', fontsize=12)
    plt.savefig('strains_clusters.png')
    plt.cla()

def save_fig_genes_clusters(pca_data, genes_id):
    colors = ['c','m','k','y']
    for i in range(4):
        v = pca_data[np.where(np.array(genes_id) == i)]
        plt.scatter(v[:, 0], v[:, 1], c=colors[i], s=10, alpha=0.5, label='group {}'.format(chr(i + 97)))
    plt.legend(loc=2)
    plt.title('2 component PCA of genes')
    plt.xlabel('pc1', fontsize=12)
    plt.ylabel('pc2', fontsize=12)
    plt.savefig('genes_clusters.png')
    plt.cla()

def save_fig_heatmap(data, order_rows, order_cols):
    sorted_data = data[order_rows].transpose()[order_cols].transpose()
    plt.subplot(2, 1, 1)
    plt.title('raw data')
    plt.imshow(data, cmap='hot')
    plt.subplot(2, 1, 2)
    plt.title('clustered data')
    plt.imshow(sorted_data, cmap='hot')
    plt.savefig('heatmap.png')
    plt.cla()


genes_count = 500
samples_count = 100
gene_groups_count = 4
csvfile = #file name here#
gene_mapping = pd.read_csv(csvfile, header=0).values[:,range(1,genes_count+1)]
pca = PCA(n_components=2)
pc_strains = pca.fit_transform(gene_mapping)
pc_genes = pca.fit_transform(gene_mapping.transpose())

# this is a very rought estimate, but its enough for our purpose
linear_params = [[0.25, 50], [0,0], [-0.25, -50]]
strains_id = get_strains_clusters(linear_params, pc_strains)
save_fig_strains_clusters(pc_strains, strains_id)

kmeans = KMeans(n_clusters = gene_groups_count).fit(pc_genes)
genes_id = kmeans.labels_
save_fig_genes_clusters(pc_genes, genes_id)

strains_order = np.array(strains_id).argsort()
genes_order = np.array(genes_id).argsort()
save_fig_heatmap(gene_mapping, strains_order, genes_order)

group_label = ['a', 'b', 'c', 'd']
f= open("genes_groups.tsv","w+")
f.write('gene\tgroup\n')
for i in range(genes_count):
    f.write('{}\t{}\n'.format(i,group_label[genes_id[i]]))
f.close()

















