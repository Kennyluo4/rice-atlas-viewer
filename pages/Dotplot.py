import streamlit as st
import matplotlib.pyplot as plt
import os, subprocess, glob, atexit

###############################################
##           Configuration                  ###
###############################################

## define function for generation figures using R
@st.cache_data
def generate_png_plot(configure_file, gene_lst_file, species, tissue):
    """
    Run the R script to generate the PNG plot.
    """
    # import hashlib

    # Create a unique hash based on input parameters
    # hash_id = hashlib.md5((species + tissue).encode()).hexdigest()
    # prefix = f"{species}_{tissue}_{hash_id}"
    prefix = f"{species}_{tissue}"
    out_fig = f"Rplot/temp/{prefix}_dotplot.tiff"
    
    cmd = ["Rscript", "Rplot/plot_dot_plot_scripts_all_species_121824.R", str(configure_file), str(gene_lst_file), str(out_fig)]
    # test = ' '.join(cmd)
    # st.write('...Running: ', test)
    
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    if result.returncode != 0:
        st.error(f"Error in R script: {result.stderr.decode('utf-8')}")
        return None
    # Return the path of the generated PNG file
    return out_fig

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
    
    
    st.sidebar.header('Plot Configuration')
    st.sidebar.markdown('## Please select a dataset to plot')
    
    ## set the plot setting using form submission: 
    with st.sidebar.form(key='Plot_config'):
        species = st.selectbox('Species', [ '---Please choose---','Osativa', 'test2'])
        tissue = st.selectbox('Tissue',['---Please choose---','bud', 'crownroot'])
        gene_list = st.text_area('Gene ID (e.g. LOC_Os10g42750)', 
                                         'LOC_Os01g06250\nLOC_Os04g55970')
        plot_h = st.slider('plot height(cm)', 1,30,10)
        plot_w = st.slider('plot width(cm)', 1,30,10)
        
        submitted = st.form_submit_button("Submit")
        
        # Clear cache when inputs change
        if submitted:
            st.cache_data.clear()  # Clear cache to avoid reusing old data
            st.write('Generating dotplot...')
            # st.write(f"Plot genes for {gene_list}")
    
    ## generate input files required by Rscript
    configs = f'ipt_target_species <- "{species}"\nipt_target_organ <- "{tissue}"\nplot_width <- {plot_w}\nplot_height <- {plot_h}\n'
    with open ('Rplot/temp/configure_file.txt', 'w') as handle:
        handle.write(configs)
        
    with open('Rplot/temp/gene_ID_list.txt', 'w') as handle:
        handle.write(gene_list)
    st.write('generated files')
    
    ###############################################
    ##                Main page                 ###
    ###############################################
    st.markdown('### Plot cell-type-specific gene expression using dotplot')
    if species != '---Please choose---' and tissue != '---Please choose---':
        fig = generate_png_plot('Rplot/temp/configure_file.txt','Rplot/temp/gene_ID_list.txt', species, tissue)
        st.image(fig)
    else:
        st.error(':point_left: Please select the dataset and genes for plotting')
    
                             
if __name__ == '__main__':
    main()

