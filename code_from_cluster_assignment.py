import pandas as pd

def load_gene_expression_data(filepath):
    """
    This function reads a gene expression data file in TSV format and preprocesses it by dropping 
    genes with missing values , transposing the DataFrame to have genes in columns and samples in rows, and 
    setting  Sample ID on the first column with name "NAME".

    Parameters:
    -----------
    filepath : str
        The path to the txt file containing the gene expression data.

    Returns:
    --------
    pandas DataFrame
        A DataFrame containing the preprocessed gene expression data.
        The DataFrame has samples in rows  and genes in columns.
    """
    genes_expression = pd.read_csv(filepath, sep="\t", index_col=0)
    genes_expression = genes_expression.dropna().T
    genes_expression.reset_index(inplace=True)
    genes_expression = genes_expression.rename_axis(None, axis=1)
    genes_expression = genes_expression.rename(columns={"index": "NAME"})
    return genes_expression
def load_clinical_data(filepath):
    """
    This function reads an Excel file and extract  estrogen receptor status data,AJCC Stage and Sample ID.
    It  only takes samples with "Negative" or "Positive" ER status and dropping samples with missing values, or other value

    Parameters:
    -----------
    filepath : str
        The path to the Excel file

    Returns:
    --------
    pandas DataFrame
        A DataFrame containing the preprocessed estrogen receptor status data.
        The DataFrame has two columns sample IDs and ER status "Negative" or "Positive"  .
    """
    estrogen_receptor = pd.read_excel(filepath, header=1)[["Complete TCGA ID","ER Status","AJCC Stage"]].dropna()
    estrogen_receptor = estrogen_receptor[estrogen_receptor["ER Status"].isin(["Negative","Positive"])]
    return estrogen_receptor
def encode_estrogen_receptor(clinical_data):
    """
    Encodes the estrogen receptor status in a given DataFrame using one-hot encoding.

    Parameters:
    -----------
    estrogen_receptor : pandas DataFrame
        A DataFrame containing the estrogen receptor status of each sample.
        The DataFrame must have a column named "ER Status" with values "Positive" or "Negative".

    Returns:
    --------
    pandas DataFrame
        A new DataFrame with o from sklearn.preprocessing import StandardScalerne-hot encoded columns for the estrogen receptor status.
        The  "Negative" values are encoded as 1
        The "Positive" values are encoded as 0
    """
    dummy_data= pd.get_dummies(clinical_data["ER Status"])
    estrogen_receptor= pd.concat((dummy_data,clinical_data),axis=1)
    estrogen_receptor =estrogen_receptor.drop(["ER Status"],axis=1)
    estrogen_receptor = estrogen_receptor.drop(["Positive"],axis=1)
    estrogen_receptor = estrogen_receptor.rename(columns={"Negative":"ER Status"})

    return estrogen_receptor[["Complete TCGA ID","ER Status","AJCC Stage"]]
def create_data_structures_with_matching_sample(genes_expression, estrogen_receptor):
    """
    This function first matches the samples in the estrogen receptor DataFrame with the samples in the
    gene expression DataFrame by finding the common sample IDs. The matching is done by modifying the
    "NAME" column in the gene expression DataFrame to match the sample IDs in the estrogen receptor
    DataFrame.

    Once the samples are matched, this function sets the index of the gene expression DataFrame to the
    matched sample IDs and reorders the rows according to the order of samples in the estrogen receptor
    DataFrame.

    Parameters:
    -----------
    genes_expression : pandas DataFrame
        A DataFrame containing the gene expression values for each sample.
        The DataFrame must have a column named "NAME" with sample IDs.

    estrogen_receptor : pandas DataFrame
        A DataFrame containing the estrogen receptor status of each sample.
        The DataFrame must have a column named "Complete TCGA ID" with sample IDs.

    Returns:
    --------
    tuple of two pandas DataFrames
        A tuple containing two DataFrames with matching samples.
        The first DataFrame contains the gene expression values, and its index is the matched sample IDs.
        The second DataFrame contains the estrogen receptor status, and its index is also the matched sample IDs.
    """
    sample_ID = estrogen_receptor["Complete TCGA ID"].to_list()
    for id in sample_ID:
        genes_expression.loc[genes_expression["NAME"].str.contains(id),"NAME"]=id
    TCGA_ID = genes_expression["NAME"].to_list()
    estrogen_receptor=estrogen_receptor.loc[estrogen_receptor["Complete TCGA ID"].isin(TCGA_ID)]
    genes_expression = genes_expression.set_index("NAME")
    genes_expression = genes_expression.reindex(index=estrogen_receptor["Complete TCGA ID"])
    return genes_expression, estrogen_receptor

def get_std(genes_expression):
    """
    Calculate the standard deviation of each gene in the expression dataset.

    Parameters:
    -----------
    genes_expression : pandas DataFrame
        A DataFrame containing gene expression data.
        Each row represents a sample, and each column represents a gene.

    Returns:
    --------
    pandas Series
        A Series containing the standard deviation of each gene in the input DataFrame.
    """
    std = genes_expression.std()
    return std

def plot_std_histogram(std, bins=100, alpha=0.5):
    """
    Plot a histogram of the pandas series.

    This function plots a histogram of a give pandas series . The histogram is saved to a file.
    
     To display it in the current figure window, you need to uncomment the last line of this function.

    Parameters:
    -----------
    std: pandas Series

    bins : int, optional
        The number of bins to use in the histogram. Default is 100.
    alpha : float, optional
        The alpha value for the histogram bars. Default is 0.5.

    Returns:
    --------
    None
    """
    ax = std.plot.hist(bins=bins, alpha=alpha)
    ax.set_title("Histogram of standard deviations of gene expression data")
    ax.set_xlabel("Standard deviations")
    ax.set_ylabel("Gene frequency")