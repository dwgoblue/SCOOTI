#%load_ext autoreload
#%autoreload 2
# import essential packages
import numpy as np
import pandas as pd
import anndata as ad
import os
import scanpy as sc
from scipy.sparse import coo_matrix
from scipy.io import mmwrite, mmread
import gzip
from tqdm import tqdm
from SCOOTI.regressionAnalyzer import *


# for MOOMIN
def posterior_probability(p_values, prior=0.05, alpha_H1=1, beta_H1=19, alpha_H0=1, beta_H0=1):
    """
    Calculate the posterior probability that a gene is differentially expressed.
    
    Parameters:
    - p_value (float): The p-value from the differential expression analysis.
    - prior (float): The prior probability of a gene being differentially expressed.
    
    Returns:
    - posterior_prob (float): The posterior probability of differential expression.
    """
    #numerator = p_value * prior
    #denominator = (p_value * prior) + ((1 - p_value) * (1 - prior))
    #posterior_prob = numerator / denominator
    import numpy as np
    from scipy.stats import beta

    # Define prior probabilities
    P_H1 = prior # alternative
    P_H0 = 1 - P_H1 # null
    
    # Calculate likelihoods
    likelihood_H1 = beta.pdf(p_values, alpha_H1, beta_H1)
    likelihood_H0 = beta.pdf(p_values, alpha_H0, beta_H0)
    
    # Calculate marginal likelihood P(D)
    marginal_likelihood = likelihood_H1 * P_H1 + likelihood_H0 * P_H0
    
    # Calculate posterior probability P(H1 | D)
    posterior_H1 = (likelihood_H1 * P_H1) / marginal_likelihood
    
    return posterior_H1

# class for single-cell datasets
class findSigGenes:
    """Acquiring up- and down-regulated genes from single-cell or bulk omics data

    The object consists the following methods.
    1) Preprocessing single-cell data: 
        read single-cell data with 10xformat (a folder with matrix.mtx, genes.tsv, and barcoodes.tsv)
    2) Convert tables of single-cell data into 10x format:
        read tables saved with .csv/.xlsx files (with chunks) and convert them into 10x format 
    3) Get significant genes via the comparisons among unique labels:
        read tables and statistically compare the differences of the gene/protein expression
    4) Get significant genes in a sample:
        instead of comparing the data across columns, the method focuses on comparing expression levels
        with z-scores or percentiles.
    
    
    Parameters
    ----------
    folder_path : {string},
        path to access the directory where the single cell data is saved.
    
    Attributes
    ----------
    folder_path : {string},
        the path to access the single-cell datasets
    adata : {anndata format},
        the data table with the format of anndata
    genedf : {pandas.DataFrame},
        the data table of genes by cells with expression levels as values

    
    Examples
    --------
    >>> from SCOOTI.GeneralMethods.findSigGenes import findSigGenes
    >>> # path to access single-cell data
    >>> path = './data/scEmbryo/GSE136714/single_cell/'
    >>> # initiate the class objective
    >>> fr = findSigGenes(path)
    >>> # read data
    >>> fr.read_scRNAseq()
    >>> # get gene expression table 
    >>> # transpose: False, index: cells, columns: genes
    >>> genedf = fr.get_genedf(transpose=False)
    >>> # Zygote
    >>> ref_cells = pd.Series(
    >>>                 genedf.index
    >>>             ).apply(lambda x: x.split('_')[0]).str.contains('Zygote').to_numpy()
    >>> # 2C
    >>> exp_cells = pd.Series(
    >>>                 genedf.index
    >>>             ).apply(lambda x: x.split('_')[0]).to_numpy()=='2cell'
    >>> # Get significant genes by the comparisons between zygote and 2cell
    >>> # method "AVGSTD" relies on 2 STD as thresholds to filter significant genes
    >>> upgene, dwgene = fr.get_transition_genes(
    >>>                 ref_cells, exp_cells, split_str='_', method='AVGSTD', save_files=True
    >>>             )
    """
    def __init__(self, folder_path):

        self.folder_path = folder_path
        self.adata = []
        self.genedf = []
    
    def rename_scFiles(self, gene_col_ind=[], barcode_col_ind=[], sampling_ind=[]):
        # get mother folder path
        fpath = self.folder_path
        for folder in os.listdir(fpath):
            if os.path.isdir(fpath+folder):
                files = os.listdir(fpath+folder)
                for file in files:
                    #print('Processing...')
                    #print(folder)
                    #print(file)
                    # process the expression table
                    if '.mtx.gz' in file:
                        try:
                            os.rename(fpath+folder+'/'+file, fpath+folder+'/matrix.mtx.gz')
                            with gzip.open(fpath+folder+'/matrix.mtx.gz', 'rb') as fmtx:
                                mtx = mmread(fmtx)
                                #print(mtx)
                                if len(sampling_ind)>0:
                                    # Convert to CSC format for efficient column slicing
                                    csc = mtx.tocsc()
                                    # select sampled columns
                                    sliced_matrix = csc[:, sampling_ind]
                                    # Convert sliced matrix back to COO or other format if needed
                                    mtx = sliced_matrix.tocoo()
                                # save the sparse matrix with .mtx format
                                mmwrite(fpath+folder+'/matrix.mtx', mtx)
                        except:
                            print('No change for mtx.gz file')
                    # process the list of genes
                    if 'genes' in file and '.tsv.gz' in file:
                        #print('process for genes', file)
                        try:
                            os.rename(fpath+folder+'/'+file, fpath+folder+'/genes.tsv.gz')
                            gdf = pd.read_csv(
                                    fpath+folder+'/genes.tsv.gz',
                                    sep='\t'
                                    )
                            # only keep a certain range of columns
                            if len(gene_col_ind)>0:
                                gdf = gdf[gdf.columns[gene_col_ind]]
                            # save genes into a tsv file
                            #print(gdf)
                            gdf.to_csv(fpath+folder+'/genes.tsv', index=False, header=True, sep='\t')
                        except:
                            print('No change for gene file...')
                    # process the list of barcodes
                    if 'barcodes' in file and '.tsv.gz' in file:
                        #print('process for barcodes', file)
                        try:
                            os.rename(fpath+folder+'/'+file, fpath+folder+'/barcodes.tsv.gz')
                            bdf = pd.read_csv(
                                    fpath+folder+'/barcodes.tsv.gz',
                                    sep='\t'
                                    )
                            # only keep a certain range of columns
                            if len(barcode_col_ind)>0:
                                bdf = bdf[bdf.columns[barcode_col_ind]]
                            # if sampling is needed, the code will only keep selected rows
                            if len(sampling_ind)>0:
                                bdf = bdf[bdf.index.isin(bdf.index[sampling_ind])]
                            # save cell barcodes into a tsv file
                            #print(bdf)
                            bdf.to_csv(fpath+folder+'/barcodes.tsv', index=False, header=True, sep='\t')
                    # not processing
                        except:
                            print(bdf[bdf.columns[barcode_col_ind]])
                            print('No change for barcode files...')
    
    # convert a table to a single-cell folder
    def table_to_10xformat(
            self,
            gene_cols=[],
            barcode_cols=[],
            suffix='',
            sep='\t',
            transpose=False,
            chunksize=1000,
            column_slice=0,
            column_slice_func=None
            ):
        
        fpath = self.folder_path
        gene_cols = np.arange(gene_cols[0], gene_cols[1])
        # iterate thru all the files under the path
        for file in os.listdir(fpath):
            if os.path.isdir(fpath+file)==0 and '.tar' not in file:
                print('Processing...')
                print(file)
                subfolder = file.split('.')[0]
                print(subfolder)
                if not os.path.exists(fpath+subfolder):
                   os.makedirs(fpath+subfolder)

                if chunksize>0:
                    from scipy.sparse import vstack
                    df_gen = pd.read_csv(fpath+file, sep=sep, chunksize=chunksize)
                    mtx_collect = []
                    gene_collect = []
                    print('Iterate thru the object')
                    print('column_slice', column_slice)
                    for df in tqdm(df_gen):
                        print(df)
                        if column_slice:
                            df = column_slice_func(df)
                            print(df.shape)
                        if transpose:
                            df = df.T
                        genes = df.iloc[:, 0]
                        barcodes = pd.DataFrame(df.iloc[:, 1:].columns)
                        exp = df.iloc[:, 1:].values
                        print(genes.shape, barcodes.shape, exp.shape)
                        # process the expression table
                        mtx = coo_matrix(exp)
                        # collect genes and mtx
                        gene_collect.append(pd.Series(genes))
                        mtx_collect.append(mtx)
                    stack_mtx = vstack(mtx_collect)
                    mmwrite(fpath+subfolder+'/matrix.mtx', stack_mtx)
                    # save genes into a dataframe and a .tsv file
                    genes = pd.concat(gene_collect)
                    print(genes, stack_mtx)
                    genes.to_csv(fpath+subfolder+'/genes.tsv', index=0, header=False, sep='\t')
                    # save barcodes into a dataframe and a .tsv file
                    barcodes.to_csv(fpath+subfolder+'/barcodes.tsv', index=0, header=False, sep='\t')

                
                else:
                    df = pd.read_csv(fpath+file, sep=sep)

                    if transpose:
                        df = df.T
                    print(df.head(5))
                    genes = df.iloc[:, gene_cols]
                    #if len(genes_cols)==1:
                    #    # convert ensemble ids to symbols
                    #    import mygene
                    #    mg = mygene.MyGeneInfo()
                    #    mg.getgenes(['ENSG00000260972.1'], scopes='entrezgene', fields='symbol')

                    #    mg.query('ENSG00000260972')#, scopes='entrezgene', fields='symbol')
                    #    mg.querymany(['ENSG00000185070.10'], scopes='entrezgene', fields='symbol')


                    #    from pyensembl import EnsemblRelease
                    #    data = EnsemblRelease(75, species='human')
                    #    geneids = genes.iloc[:,0].to_numpy()
                    #    symbols = []
                    #    for gene_id in tqdm(geneids):
                    #        try:
                    #            symbols.append(data.gene_by_id(gene_id).gene_name)
                    #        except:
                    #            symbols.append(gene_id)
                    #    # replace ids with symbols
                    #    genes['gene_symbols'] = symbols

                    barcode_sel = np.arange(barcode_cols[0], len(df.columns)+barcode_cols[1])
                    barcodes = pd.DataFrame(df.iloc[:, barcode_sel].columns)
                    exp = df.iloc[:, barcode_sel].values
                    print(genes.shape, barcodes.shape, exp.shape)
                    # process the expression table
                    mtx = coo_matrix(exp)
                    mmwrite(fpath+subfolder+'/matrix.mtx', mtx)
                    # save genes into a dataframe and a .tsv file
                    genes.to_csv(fpath+subfolder+'/genes.tsv', index=0, header=False, sep='\t')
                    # save barcodes into a dataframe and a .tsv file
                    barcodes.to_csv(fpath+subfolder+'/barcodes.tsv', index=0, header=False, sep='\t')

    # convert a processed dataframe to a single-cell folder
    def df_to_10xformat(
            self,
            df,
            prefix='',
            ):
        # get essential values
        genes = df.index
        barcodes = df.columns
        expression = df.values
        # get mother folder path
        fpath = self.folder_path
        # subfolder name
        dir10x = self.folder_path+'/'+prefix+'_10x/'
        # make sure the subfolder is existing ornot
        isExist = os.path.exists(dir10x)
        if not isExist:
            # create a new dir if not existing
            os.makedirs(dir10x)
            print(f'Create a folder {dir10x}')
        # process the expression table
        try:
            # Convert sliced matrix back to COO or other format if needed
            mtx = coo_matrix(expression)
            # save the sparse matrix with .mtx format
            mmwrite(dir10x+'/matrix.mtx', mtx)
        except:
            print('No change for mtx.gz file')
        # process the list of genes
        try:
            # save genelist
            genes = pd.DataFrame(genes)
            genes.to_csv(dir10x+'/genes.tsv', index=False, header=False, sep='\t')
        except:
            print('No change for gene file...')
        # process the list of barcodes
        try:
            # save cell barcodes into a tsv file
            barcodes = pd.DataFrame(barcodes)
            barcodes.to_csv(dir10x+'/barcodes.tsv', index=False, header=False, sep='\t')
        except:
            print('No change for barcode files...')
    


    # fuion to read single scRNAseq data
    def read_scRNAseq(self, folder='', batch_name='', rename_cellNames=False):
    
        import scanpy as sc

        if folder=='':
            access_path = self.folder_path
        else:
            access_path = folder

        print('initiation', access_path)
        # loading
        adata = sc.read_10x_mtx(
            access_path,  # the directory with the `.mtx` file
            var_names='gene_symbols', # use gene symbols for the variable names (variables-axis index)
            cache=False,
            ) # write a cache file for faster subsequent reading
        print(adata.shape)
        adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
        
        # 1st filtering
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        
        adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(
                adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True
                )
        
        # 2nd filtering
        #adata = adata[adata.obs.n_genes_by_counts < 5000, :]
        adata = adata[adata.obs.pct_counts_mt < 5, :]
        
        # normalize data
        sc.pp.normalize_total(adata, target_sum=1e4)
    
        # Logarithmize the data
        sc.pp.log1p(adata)
            
        # update
        adata.obs['batch'] = [batch_name]*adata.shape[0]
        if rename_cellNames == True:
            adata.obs.index = np.arange(len(adata.obs.index))
            adata.obs.index = pd.Series(adata.obs.index).apply(lambda x: f'{batch_name}_{x}')
        self.adata = adata
        print('Batch name assigned')

        return adata
    

    # read multiple rna seq datasets
    def batch_read_scRNAseq(self, folders, batch_names, rename_cellNames):
        # get the path from class attribute
        #folders = os.listdir(self.folder_path)
        # read files
        for folder, batch_name in zip(folders, batch_names):
            yield (
                    self.read_scRNAseq(
                        folder+'/',
                        batch_name=batch_name,
                        rename_cellNames=rename_cellNames,
                        )
                    )
        
    # merge adata
    def merge_adata(self, rename_cellNames=False):
        # get files from the path
        folders = [
                self.folder_path+f for f in os.listdir(self.folder_path) if os.path.isdir(self.folder_path+f) and f!='sigGenes'
                ]
        print(folders)
        batches = [folder.split('/')[-1] for folder in folders]
        # concatenate anndata
        scRes = ad.concat(
                self.batch_read_scRNAseq(
                    folders,
                    batch_names=batches,
                    rename_cellNames=rename_cellNames
                    ),
                axis=0,
                join='inner'
                )
        print(batches)
        print(scRes.to_df())
        print('Successfully merge the scRNAseq datasets!')
        sc.pp.combat(scRes)
        self.genedf = scRes.to_df()

    # convert adata to table (gene by cells)
    def get_genedf(self, transpose=False):

        # load data 
        res = self.adata

        marker_genes = res.var.index.to_numpy()
        genedf = sc.get.obs_df(
                res,
                keys=[*marker_genes]
            )
        if transpose==True:
            genedf = genedf.T
            self.genedf = genedf
        else:
            self.genedf = genedf


        return self.genedf
    
    # apply for tables come with metadata
    def get_genes_from_processed_data(self, csv_suffix, metadata_suffix):

        """
        Parameters
        ----------
        csv_suffix (str): like '_magic_imputed.csv'
        metadata_suffix (str): like '_prep_metadata.csv' paired with csv_suffix
        """

        # get path
        path = self.folder_path
        # get files under the path
        files = os.listdir(path)
        df_collect = [] # empty list for collecting anndata
        batches = [] # empty list for collecting batch names
        # iterate thru the files
        for file in files:
            if csv_suffix in file:
                fname = file.split(csv_suffix)[0]
                # read metadata 
                metadata = pd.read_csv(path+fname+metadata_suffix, index_col=0)
                id_map = dict(zip(metadata.index, np.arange(len(metadata.index))))
                induction_map = dict(zip(metadata.index, metadata.Time))
                # read csv files
                exp = pd.read_csv(path+file, index_col=0)
                exp.columns = pd.Series(exp.columns).apply(
                        lambda x: '{0}_{1}_{2}'.format(fname, id_map[x], induction_map[x])
                        )
                # convert dataframe into anndata
                adata = ad.AnnData(exp.T) 
                batches = batches+[fname]*len(exp.columns)
                #adata = sc.AnnData(exp, exp.index.to_frame(), exp.columns.to_frame())
                #display(adata.obs)
                df_collect.append(adata)


        # ![issue] inner-join drops genes that are not cross-listed in different cells
        # Think of imputation in the future
        scRes = ad.concat(df_collect, axis=0, join='inner')
        scRes.obs['batch'] = batches
        sc.pp.combat(scRes)
        self.genedf = scRes.to_df()

        return scRes

    # apply to tables without metadata and remove batch effects with pyCombat
    def get_genes_from_tables_pyCombat(self, dfs, transpose=True, plot=False):

        """
        Parameters
        ----------
        dfs (dictionary): e.g. {'tissue1':df1, 'tissue2':df2}
        """

        # import libraries
        from combat.pycombat import pycombat
       
        # get path
        df_collect = [] # empty list for collecting anndata
        batches = []
        # iterate thru the files
        for k, df in dfs.items():
            # convert dataframe into anndata
            if transpose==True:
                data = df
                batches = batches+[k]*len(df.columns)
            else:
                data = df.T
                batches = batches+[k]*len(df.index)
            print(data)
            df_collect.append(data)

        # ![issue] inner-join drops genes that are not cross-listed in different cells
        # Think of imputation in the future
        res = pd.concat(df_collect, axis=1, join='inner')
        print(res)
        # run pyComBat
        df_corrected = pycombat(res, batches)
        self.genedf = df_corrected.T

        if plot:
            import seaborn as sns
            import matplotlib.pyplot as plt
            # plot raw data
            fig, ax = plt.subplots(1,1,figsize=(50, 8))
            #plot_df = res.copy().melt()
    
            sns.boxplot(data=res)#, x=batches, hue=batches)
            plt.xticks(rotation=90)
            plt.savefig('/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/combined_expression.png')

            # visualise results
            fig, ax = plt.subplots(1,1,figsize=(50, 8))
            sns.boxplot(data=df_corrected)#, hue=batches)
            plt.xticks(rotation=90)
            plt.savefig('/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/combined_expression_combat.png')

        return res

    # apply to tables without metadata
    def get_genes_from_tables(self, dfs, transpose=True):

        """
        Parameters
        ----------
        dfs (dictionary): e.g. {'tissue1':df1, 'tissue2':df2}
        """

        # get path
        df_collect = [] # empty list for collecting anndata
        batches = []
        # iterate thru the files
        for k, df in dfs.items():
            # convert dataframe into anndata
            if transpose==True:
                adata = ad.AnnData(df.T)
                adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
                batches = batches+[k]*len(df.columns)
            else:
                adata = ad.AnnData(df)
                adata.var_names_make_unique()  # this is unnecessary if using `var_names='gene_ids'` in `sc.read_10x_mtx`
                batches = batches+[k]*len(df.index)
            print(adata)
            df_collect.append(adata)


        # ![issue] inner-join drops genes that are not cross-listed in different cells
        # Think of imputation in the future
        scRes = ad.concat(df_collect, axis=0, join='inner')
        print(scRes)
        scRes.obs['batch'] = batches
        sc.pp.combat(scRes)
        self.genedf = scRes.to_df()

        return scRes
        #return df_collect

    def regulator_clustering(self, upgenedf, dwgenedf, method, correction, alpha):

        import seaborn as sns
        import matplotlib.pyplot as plt
        # check the similarity of the genes
        glist = pd.read_csv('/home/daweilin/StemCell/glist.csv')
        cols = dwgenedf.columns[pd.Series(dwgenedf.columns).isin(glist['Var1'])]
        
        plot_dw = dwgenedf[cols].copy()
        plot_up = upgenedf[cols].copy()

        labels = pd.Series(dwgenedf.index).apply(lambda x: x.split('_')[0])

        cf = clustering_func(
                plot_dw.T,
                '/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/',
                f'cellCycle_dwGene_{method}_{alpha}_{correction}',
                mets_category(),
            )
        # show correlations
        cf.corr_clustermap(labels)

        cf = clustering_func(
                plot_up.T,
                '/home/daweilin/StemCell/Project_mESC_JinZhang/regressor_results_new/',
                f'cellCycle_upGene_{method}_{alpha}_{correction}',
                mets_category(),
            )
        # show correlations
        cf.corr_clustermap(labels)

    
    # get up- and down-regulated genes with Z-scores
    def get_regulators_by_zscores(self, th=2, split_str='_'):

        import scipy.stats as st
        # get expression table
        genedf = self.genedf
        print('shape of gene expression data:', genedf.shape)
        # create 95% confidence interval for population mean weight
        genedf = (genedf.sub(genedf.mean(axis=1), axis=0)).div(genedf.std(axis=1), axis=0)
        #print(genedf)
        #print('shape of limits:', ub.shape)
        # get up- or down-regulated genes for each cells
        dwgenedf = genedf.apply(lambda x: x<-th, axis=1)
        upgenedf = genedf.apply(lambda x: x>th, axis=1)
        
        # save files
        for i in range(len(upgenedf)):

            # get cell names
            prefix = upgenedf.index[i]
            cell_dir = self.folder_path+'sigGenes/'+prefix.split(split_str)[0]+'/'
            isExist = os.path.exists(cell_dir)
            if not isExist:
                # create a new dir if not existing
                os.makedirs(cell_dir)
                print(f'Create a folder for {prefix.split("_")[0]}')
            # save up-/down-regulated genes for each cells
            pd.DataFrame({
                'upgenes':upgenedf.iloc[i,:][upgenedf.iloc[i, :]==True].index.to_numpy()
                }).to_csv(cell_dir+prefix+'_upgenes.csv')

            pd.DataFrame({
                'dwgenes':dwgenedf.iloc[i,:][dwgenedf.iloc[i, :]==True].index.to_numpy()
                }).to_csv(cell_dir+prefix+'_dwgenes.csv')


        return genedf

    # methods to get sig genes during transitions
    def get_compared_genes(self, alpha=0.95, split_str='_', prefix_define='', method='CI', std_num=2, correction=False, save_files=False):
        
        import scipy.stats as st
        import statsmodels.stats.multitest as mt
        # get expression table
        genedf = self.genedf
        print(genedf)
        print('shape of gene expression data:', genedf.shape)
        if method=='CI':
            # create 95% confidence interval for population mean weight
            ub, lb = st.t.interval(
                        alpha=alpha,
                        df=len(genedf)-1,
                        loc=np.mean(genedf),
                        scale=st.sem(genedf)
                    )
            print('shape of limits:', ub.shape)

        else: # AVGSTD
            # create mean+2std for population mean weight
            ub = genedf.mean(axis=0)+std_num*genedf.std(axis=0)
            lb = genedf.mean(axis=0)-std_num*genedf.std(axis=0)
            print('shape of limits:', ub.shape)
        # get up- or down-regulated genes for each cells
        dwgenedf = genedf.apply(lambda x: x<lb, axis=1)
        upgenedf = genedf.apply(lambda x: x>ub, axis=1)

        if correction=='True':
            zscores = genedf.sub(genedf.mean(axis=0), axis=1).div(genedf.std(axis=0), axis=1)
            pvalues = zscores.apply(lambda x: st.norm.sf(abs(x))*2)
            corrected_p = {}
            for i in range(len(pvalues)):
                corrected_p[pvalues.index[i]] = mt.fdrcorrection(pvalues.iloc[i,:])[1]
            corrected_p = pd.DataFrame(corrected_p).T
            print('FDR')
            print(corrected_p)
            dwgenedf[corrected_p>=0.05] = 0
            upgenedf[corrected_p>=0.05] = 0

        print(upgenedf)
        
        self.regulator_clustering(upgenedf, dwgenedf, method, correction, alpha)
        
        if save_files==True:
            # save files
            for i in tqdm(range(len(upgenedf))):
                # get cell names
                prefix = upgenedf.index[i] if prefix_define=='' else prefix_define
                cell_dir = self.folder_path+'sigGenes/'+prefix.split(split_str)[0]+'/'
                isExist = os.path.exists(cell_dir)
                if not isExist:
                    # create a new dir if not existing
                    os.makedirs(cell_dir)
                    print(f'Create a folder for {prefix.split("_")[0]}')
                # save up-/down-regulated genes for each cells
                pd.DataFrame({
                    'upgenes':upgenedf.iloc[i,:][upgenedf.iloc[i, :]==True].index.to_numpy()
                    }).to_csv(cell_dir+prefix+'_upgenes.csv')

                pd.DataFrame({
                    'dwgenes':dwgenedf.iloc[i,:][dwgenedf.iloc[i, :]==True].index.to_numpy()
                    }).to_csv(cell_dir+prefix+'_dwgenes.csv')


        return upgenedf, dwgenedf


    # methods to get sig genes during transitions
    def get_transition_genes(self, ref_cells, exp_cells, alpha=0.95, split_str='_', prefix_define='', method='CI', std_num=2, correction=False, save_files=False):
        
        """
        Parameter
        ---------
        ref_cells (numpy.array): a boolean array that indicates where is the reference rows
        exp_cells (numpy.array): a boolean array that indicates where is the compared rows
        alpha (float): threshold for statistical tests
        split_str (str): split the name of cells
        
        Return
        ------
        genedf (pandas.DataFrame): gene expression table

        """

        import scipy.stats as st
        # get expression table
        genedf = self.genedf
        print('shape of gene expression data:', genedf.shape)
        if method=='CI':
            # create 95% confidence interval for population mean weight
            lb, ub = st.t.interval(
                        alpha=alpha,
                        df=len(genedf)-1,
                        loc=np.mean(genedf),
                        scale=st.sem(genedf)
                    )
            print('shape of limits:', ub.shape)
        else: # AVGSTD
            # create mean+2std for population mean weight
            ub = genedf.mean(axis=0)+std_num*genedf.std(axis=0)
            lb = genedf.mean(axis=0)-std_num*genedf.std(axis=0)
            print('shape of limits:', ub.shape)
        # get up- or down-regulated genes for each cells
        dwgenedf = genedf[exp_cells].apply(lambda x: x<lb, axis=1)
        upgenedf = genedf[exp_cells].apply(lambda x: x>ub, axis=1)

        if correction=='True':
            zscores = genedf.sub(genedf.mean(axis=0), axis=1).div(genedf.std(axis=0), axis=1)
            pvalues = zscores.apply(lambda x: st.norm.sf(abs(x))*2)
            corrected_p = {}
            for i in range(len(pvalues)):
                corrected_p[pvalues.index[i]] = mt.fdrcorrection(pvalues.iloc[i,:])[1]
            corrected_p = pd.DataFrame(corrected_p).T
            print('FDR')
            print(corrected_p)
            dwgenedf[corrected_p>=0.05] = 0
            upgenedf[corrected_p>=0.05] = 0

        print(upgenedf)


        self.regulator_clustering(upgenedf, dwgenedf, method, correction, alpha)
        
        if save_files==True:
            # save files
            for i in range(len(upgenedf)):
                # get cell names
                prefix = upgenedf.index[i] if prefix_define=='' else prefix_define
                cell_dir = self.folder_path+'sigGenes/'+prefix.split(split_str)[0]+'/'
                isExist = os.path.exists(cell_dir)
                if not isExist:
                    # create a new dir if not existing
                    os.makedirs(cell_dir)
                    print(f'Create a folder for {prefix.split("_")[0]}')
                # save up-/down-regulated genes for each cells
                pd.DataFrame({
                    'upgenes':upgenedf.iloc[i,:][upgenedf.iloc[i, :]==True].index.to_numpy()
                    }).to_csv(cell_dir+upgenedf.index[i]+'_upgenes.csv')

                pd.DataFrame({
                    'dwgenes':dwgenedf.iloc[i,:][dwgenedf.iloc[i, :]==True].index.to_numpy()
                    }).to_csv(cell_dir+upgenedf.index[i]+'_dwgenes.csv')


        return upgenedf, dwgenedf


    # methods to get sig genes during transitions
    def get_top_last_genes(self, ratio=0.1, split_str='_', prefix_define='', flip=False, save_files=True, zscores=True, th=2, chunk=0):
        
        import scipy.stats as st
        if chunk==1:
            dwgenedf_collect = []
            upgenedf_collect = []
            for i in range(len(self.adata.obs.index)):
                col = self.adata.obs.index[i]
                print(col)
                genedf = pd.DataFrame(sc.get.var_df(self.adata, keys=col))
                #print(genedf)
                #print('zscores', zscores)

                if zscores==True:
                    # create 95% confidence interval for population mean weight
                    genedf = genedf.sub(genedf.mean()).div(genedf.std())
                    # get top/last n values for ub and lb
                    lb, ub = -th, th
                else:
                    # get top/last n values for ub and lb
                    lb = np.nanpercentile(genedf[col], ratio*100)
                    ub = np.nanpercentile(genedf[col], (1-ratio)*100)
                    print('Bounds:', lb, ub, min(genedf.min()), max(genedf.max()))
                    #print(lb.shape)

                if min(genedf.min())==lb:
                    print('Special case')
                    # get up- or down-regulated genes for each cells
                    dwgenedf = pd.DataFrame(genedf.apply(lambda x: x<=lb))
                    upgenedf = pd.DataFrame(genedf.apply(lambda x: x>ub))
                else:
                    # get up- or down-regulated genes for each cells
                    dwgenedf = pd.DataFrame(genedf.apply(lambda x: x<lb))
                    upgenedf = pd.DataFrame(genedf.apply(lambda x: x>ub))
                # save binary results of gene types
                dwgenedf_collect.append(dwgenedf)
                upgenedf_collect.append(upgenedf)
                #print(dwgenedf.sum(), lb, ub)

                # save files
                if save_files==True:
                    # get cell names
                    prefix = col if prefix_define=='' else prefix_define
                    cell_dir = self.folder_path+'sigGenes/'+prefix.split(split_str)[0]+'/'
                    isExist = os.path.exists(cell_dir)
                    if not isExist:
                        # create a new dir if not existing
                        os.makedirs(cell_dir)
                        print(f'Create a folder for {prefix.split("_")[0]}')

                    if flip==True:
                        # save up-/down-regulated genes for each cells
                        pd.DataFrame({
                            'upgenes':dwgenedf[dwgenedf.iloc[:,0]==True].index.to_numpy()
                            }).to_csv(cell_dir+prefix+f'_{col}_upgenes.csv')

                        pd.DataFrame({
                            'dwgenes':upgenedf[upgenedf.iloc[:,0]==True].index.to_numpy()
                            }).to_csv(cell_dir+prefix+f'_{col}_dwgenes.csv')
                    else:
                        # save up-/down-regulated genes for each cells
                        pd.DataFrame({
                            'upgenes':upgenedf[upgenedf.iloc[:,0]==True].index.to_numpy()
                            }).to_csv(cell_dir+prefix+f'_{col}_upgenes.csv')

                        pd.DataFrame({
                            'dwgenes':dwgenedf[dwgenedf.iloc[:,0]==True].index.to_numpy()
                            }).to_csv(cell_dir+prefix+f'_{col}_dwgenes.csv')
                        #print(dwgenedf[col])
                        #print(dwgenedf[dwgenedf[col]==True].index.to_numpy())
                        #print(upgenedf[upgenedf[col]==True].index.to_numpy())

            return pd.concat(upgenedf_collect), pd.concat(dwgenedf_collect)

        else:

            # get expression table
            genedf = self.genedf # row: cells, columns: genes
            print('shape of gene expression data:', genedf.shape)
            print(genedf)
            
            if zscores==True:
                # create 95% confidence interval for population mean weight
                genedf = (genedf.sub(genedf.mean(axis=1), axis=0)).div(genedf.std(axis=1), axis=0)
                # get top/last n values for ub and lb
                lb, ub = -th, th
            else:
                # get top/last n values for ub and lb
                lb = genedf.apply(lambda x: np.nanpercentile(
                    x,
                    ratio*100
                    ), axis=1).to_numpy()
                ub = genedf.apply(lambda x: np.nanpercentile(
                    x,
                    (1-ratio)*100
                    ), axis=1).to_numpy()
                print(lb.shape)

            # get up- or down-regulated genes for each cells
            dwgenedf = genedf.apply(lambda x: x<lb, axis=0)
            upgenedf = genedf.apply(lambda x: x>ub, axis=0)

            # save files
            if save_files==True:
                # save files
                for i in range(len(upgenedf)):
                    # get cell names
                    col = upgenedf.index[i]
                    # get cell names
                    prefix = upgenedf.index[i] if prefix_define=='' else prefix_define
                    cell_dir = self.folder_path+'sigGenes/'+prefix.split(split_str)[0]+'/'
                    isExist = os.path.exists(cell_dir)
                    if not isExist:
                        # create a new dir if not existing
                        os.makedirs(cell_dir)
                        print(f'Create a folder for {prefix.split("_")[0]}')

                    if flip==True:
                        # save up-/down-regulated genes for each cells
                        pd.DataFrame({
                            'upgenes':dwgenedf.iloc[i,:][dwgenedf.iloc[i, :]==True].index.to_numpy()
                            }).to_csv(cell_dir+prefix+f'_{col}_upgenes.csv')

                        pd.DataFrame({
                            'dwgenes':upgenedf.iloc[i,:][upgenedf.iloc[i, :]==True].index.to_numpy()
                            }).to_csv(cell_dir+prefix+f'_{col}_dwgenes.csv')
                    else:
                        # save up-/down-regulated genes for each cells
                        pd.DataFrame({
                            'upgenes':upgenedf.iloc[i,:][upgenedf.iloc[i, :]==True].index.to_numpy()
                            }).to_csv(cell_dir+prefix+f'_{col}_upgenes.csv')

                        pd.DataFrame({
                            'dwgenes':dwgenedf.iloc[i,:][dwgenedf.iloc[i, :]==True].index.to_numpy()
                            }).to_csv(cell_dir+prefix+f'_{col}_dwgenes.csv')



            return upgenedf, dwgenedf

        ## save files
        #upgenes_dict = []
        #dwgenes_dict = []
        #for col in genedf.columns:
        #    # get cell names
        #    prefix = col if prefix=='' else prefix
        #    cell_dir = self.folder_path+'sigGenes/'+prefix.split(split_str)[0]+'/'
        #    isExist = os.path.exists(cell_dir)
        #    if not isExist:
        #        # create a new dir if not existing
        #        os.makedirs(cell_dir)
        #        print(f'Create a folder for {prefix.split(split_str)[0]}')

        #    if flip==True:
        #        # save up-/down-regulated genes for each cells
        #        dwgenes = pd.DataFrame({
        #            'dwgenes':genedf.nlargest(n, col, keep='all').index.to_numpy()
        #            })
        #        upgenes = pd.DataFrame({
        #            'upgenes':genedf.nsmallest(n, col, keep='all').index.to_numpy()
        #            })

        #    else:
        #        # save up-/down-regulated genes for each cells
        #        upgenes = pd.DataFrame({
        #            'upgenes':genedf.nlargest(n, col, keep='all').index.to_numpy()
        #            })
        #        dwgenes = pd.DataFrame({
        #            'dwgenes':genedf.nsmallest(n, col, keep='all').index.to_numpy()
        #            })

        #    if save_files==True:
        #        dwgenes.to_csv(cell_dir+prefix+f'_{col}_dwgenes.csv')
        #        upgenes.to_csv(cell_dir+prefix+f'_{col}_upgenes.csv')
        #    
        #    # output
        #    upgenes.columns = [col]
        #    dwgenes.columns = [col]
        #    upgenes_dict.append(upgenes)
        #    dwgenes_dict.append(dwgenes)

        #return upgenes_dict, dwgenes_dict

    # get random sampling of up- and down-gene
    def random_sample_regulators(self, rep=1000):
        # shuffle the index
        geneList = self.genedf.columns
        genedf = self.genedf.T
        # shuffle the index
        def df_shuffle(df, i):
            df = df.sample(frac=1)
            df.columns = pd.Series(df.columns).apply(lambda x: f'{x}-{i}')
            return df
        genedf = pd.concat([df_shuffle(genedf, i) for i in range(rep)], axis=1)
        genedf.index = geneList
        genedf = genedf.T
        self.genedf = genedf

