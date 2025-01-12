import streamlit as st
# import matplotlib.pyplot as plt
import subprocess
import os, subprocess, glob, atexit
        
###############################################
##           Configuration                  ###
###############################################

## define function for generation figures using R
@st.cache_data
def generate_png_plot(species, feature, trajectory,candidateID):
    """
    Run the R script to generate the PNG plot.
    """
    import hashlib

    out_path = 'Rplot/temp'

    cmd = ["Rscript", "Rplot/plot_trajectory_genes_TFs_motifs_ACRs_only_heatmap_010925_HD.R", str(species), str(feature),str(trajectory), str(candidateID), str(out_path)]
    # test = ' '.join(cmd)
    # st.write('...Running: ', test)
    
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    if result.returncode != 0:
        st.error(f"Error in R script: {result.stderr.decode('utf-8')}")
        return None
    # Return the path of the generated tiff file
    
    out_fig = f"{out_path}/{trajectory}.trajectory.{feature}_heatmap.tiff"     ## same as the output name from Rscript
    
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
    
    ## set the plot settings: 
    species = st.sidebar.selectbox('Species', [ '---Please choose---','Osativa', 'test2'])
    trajectory = st.sidebar.selectbox('Trajectory',['---Please choose---','bud_ProcamPhloem', 'bud_ProcamXylem', 'Eseedling_SAMPhloem', 'Eseedling_SAMPhloemCC', 'Eseedling_SAMXylemPrecursor', 'crownroot_QCCortex', 'crownroot_QCMetaXylem', 'crownroot_QCProtoXylem','panicle_PBMSBMSMFM','panicle_PBMSMFM', 'semroot_QCCortex', 'semroot_QCMetaXylem', 'semroot_QCProtoXylem'])
    # feature = st.selectbox('Feature', [ '---Please choose---','Gene', 'TF', 'Mt', 'ACR'])
    feature = st.sidebar.radio('Feature', ['Gene', 'TF', 'Mt', 'ACR'], horizontal=True)
    ## Retrive the featurelist after feature type selection
    feature_file = f'Rplot/data/trajectory_features/opt_{feature}_feature_list.txt'
    featureIDs = [id.strip() for id in open(feature_file).readlines()]
    featureName = st.sidebar.selectbox('Enter feature name to highlight', ['', *featureIDs])
    
    # plot_h = st.slider('plot height(cm)', 1,30,10)
    # plot_w = st.slider('plot width(cm)', 1,30,10)

    ## Create a submit button, plot only when click submit
    with st.sidebar.form(key='my_form'):
        st.write('Click to generate the plot')
        submitted = st.form_submit_button("Submit")

    
    ###############################################
    ##                Main page                 ###
    ###############################################
    st.markdown('### Plot the trajectory')
    if submitted:
        st.cache_data.clear()  # Clear cache to avoid reusing old data
        if species != '---Please choose---' and feature != '---Please choose---' and trajectory != '---Please choose---':
            fig = generate_png_plot(species, feature, trajectory, featureName)
            st.image(fig)
        else:
            st.error(':point_left: Please select the dataset and genes for plotting')
        
if __name__ == '__main__':
    main()


