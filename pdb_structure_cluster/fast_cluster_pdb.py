from pathlib import Path
import numpy as np
import mdtraj as md
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, PDBIO
from multiprocessing import Pool
import pickle
import itertools
import os
import shutil
# Step 1: Load PDB files using pathlib
def load_pdb_files(pdb_directory, N=None, ca_only=True):
    """Load PDB files using pathlib and select Cα atoms for efficiency."""
    pdb_directory = Path(pdb_directory)
    pdb_files = list(pdb_directory.glob('*.cif'))
    if N is not None:
        pdb_files = pdb_files[:N]
    trajectories = []
    for pdb_file in pdb_files:
        traj = md.load(str(pdb_file))  # Convert Path to string for mdtraj
        if ca_only:
            traj = traj.atom_slice(traj.topology.select('name CA'))  # Select Cα atoms only
        trajectories.append(traj)
    return trajectories, [str(p) for p in pdb_files]  # Convert Path objects to strings for compatibility

# Step 2: Parallel RMSD calculation with superimposition
def compute_pairwise_rmsd(args):
    """Compute RMSD between two structures after superimposition."""
    i, j, trajectories = args
    # Superimpose the second structure onto the first based on Cα atoms
    traj_i = trajectories[i]
    traj_j = trajectories[j]
    traj_j.superpose(traj_i, atom_indices=traj_i.topology.select('name CA'))
    # Compute RMSD
    rmsd = md.rmsd(traj_i, traj_j, precentered=True)[0]
    return i, j, rmsd

def compute_rmsd_matrix(trajectories, n_processes=-1, cache_file='rmsd_matrix.pkl'):
    """Compute RMSD distance matrix using multiprocessing."""
    cache_file = Path(cache_file)
    if cache_file.exists():
        print("Loading cached RMSD matrix...")
        with open(cache_file, 'rb') as f:
            return pickle.load(f)
    
    n_structures = len(trajectories)
    rmsd_matrix = np.zeros((n_structures, n_structures))
    
    # Generate pairs for parallel computation
    pairs = list(itertools.combinations(range(n_structures), 2))
    
    with Pool(processes=n_processes) as pool:
        results = pool.map(compute_pairwise_rmsd, [(i, j, trajectories) for i, j in pairs])
    
    # Fill the symmetric matrix
    for i, j, rmsd in results:
        rmsd_matrix[i, j] = rmsd
        rmsd_matrix[j, i] = rmsd
    
    # Cache the RMSD matrix
    with open(cache_file, 'wb') as f:
        pickle.dump(rmsd_matrix, f)
    
    return rmsd_matrix

# Step 3: Perform K-Means clustering with specified number of clusters
def perform_clustering(rmsd_matrix, n_clusters=5):
    """Perform K-Means clustering to generate a specified number of clusters."""
    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    clusters = kmeans.fit_predict(rmsd_matrix) + 1  # Add 1 to start cluster IDs from 1
    # print(clusters)
    return clusters

# Step 4: Identify representative structures
def get_representative_structures(clusters, pdb_files, output_dir='Cluster_PDB'):
    """Identify and save representative structure for each cluster, return sorted representatives."""
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    unique_clusters = np.unique(clusters)
    pdb_file_name = []
    
    for cluster_id in unique_clusters:
        output_dir_cluster = output_dir / f"cluster_{cluster_id}"
        output_dir_cluster.mkdir(exist_ok=True)
        cluster_indices = np.where(clusters == cluster_id)[0]
        for cluster_indice in cluster_indices:
            cif_file = Path(pdb_files[cluster_indice])
            shutil.copy(cif_file,output_dir_cluster)
            pdb_file_name.append(cif_file.name)
        
        print(f"Cluster {cluster_id} (size {len(cluster_indices)}): {','.join(pdb_file_name)}")
    
    return

# Main function to run the clustering pipeline
def main(pdb_directory, N=None, n_clusters=5, n_processes=None, ca_only=True):
    # Load PDB files
    print("Loading PDB files...")
    trajectories, pdb_files = load_pdb_files(pdb_directory, N, ca_only=ca_only)
    
    # Compute RMSD matrix
    print("Computing RMSD distance matrix...")
    global rmsd_matrix  # Make rmsd_matrix global for use in get_representative_structures
    rmsd_matrix = compute_rmsd_matrix(trajectories, n_processes=n_processes)
    
    # Perform clustering
    print(f"Performing K-Means clustering to generate {n_clusters} clusters...")
    clusters = perform_clustering(rmsd_matrix, n_clusters=n_clusters)
    
    # Get representative structures
    print("Extracting representative structures...")
    get_representative_structures(clusters, pdb_files)
    
    # # Output sorted representative structures
    # print(f"\nFound {len(np.unique(clusters))} clusters.")
    # print("Representative structures (sorted by cluster ID):")
    # for cluster_id, rep_file in representatives:
    #     print(f"Cluster {cluster_id}: {rep_file}")

if __name__ == "__main__":
    # 两个参数，1. 选择前多少个结构聚类，2. 聚成多少类
    pdb_directory = "pdbs"
    N = 10 #选择前10个结构聚类
    N = None #选择所有的结构聚类
    n_clusters = 5  # Number of desired clusters
    n_processes = 4 # 使用多少个节点进行计算
    ca_only = True  # Use only Cα atoms for RMSD calculation
    main(pdb_directory, N, n_clusters, n_processes, ca_only)


