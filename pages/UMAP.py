import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import os, glob, atexit

# Set the page config to wide mode
# st.set_page_config(layout="wide")
    
## Function to clean up temporary files
def cleanup_temp_folder():
    """
    Deletes all files in the specified temp folder.
    """
    folder = "Rplot/temp"
    files = glob.glob(os.path.join(folder, "*"))
    for file in files:
        try:
            os.remove(file)
        except Exception as e:
            st.warning(f"Could not delete file {file}: {e}")

# Register cleanup function to run on app exit
atexit.register(cleanup_temp_folder)

def main():

    ###############################################
    ##                 Sidebar                  ###
    ###############################################
    # with st.sidebar.form():           # use if want to add a Submit button before change everything
    
    st.sidebar.header('Plot Configuration')
    st.sidebar.markdown('## Please select a dataset:')
    dataset = st.sidebar.selectbox('Dataset', [ '---Please choose---','Late_seed', 'Leaf', 'Seedling','Early_seedling', 'Seminal_root', 'Panicle', 'Early_seed', 'Crown_root', 'Axillary_bud'])
    feature = st.sidebar.radio('Feature', ['gene', 'motif'], horizontal=True)
    
    ## Retrive the adata after data type selection
    if dataset != None and dataset != '---Please choose---':
        filename = f'../data/sobjs/{dataset}_{feature}_adata.h5ad'
        adata = sc.read_h5ad(filename)

        # read the genes
        featureIDs = adata.var.index.tolist()
        featureName = st.sidebar.selectbox('Enter gene name for expression plot', ['', *featureIDs])
    else:
        st.sidebar.write('Please select a data to explore the genes')

    ###############################################
    ##                Main page                 ###
    ###############################################    
    st.markdown('### Plot gene expression on single-cell UMAP')
    
    if dataset == None or dataset == '---Please choose---':
        st.error(':point_left: Please select the dataset and genes for plotting')
        return
    # select features for umap plot (feed geneID and cell annotation as feature to Scanpy)
    else:
        variables_to_plot = ['Final_annotation_TCP_up']        ## default cell type to plot
        if featureName:
            variables_to_plot.append(featureName)
            
        fig, axs = plt.subplots(len(variables_to_plot), 1, figsize=(5, 5 * len(variables_to_plot)))
        if len(variables_to_plot) == 1:
            axs = [axs]
        for ax, gene in zip(axs, variables_to_plot):
            sc.pl.umap(adata, color=gene, ax=ax, show=False)
            plt.subplots_adjust(wspace=1.2)
        st.pyplot(fig)

if __name__ == '__main__':
    main()
