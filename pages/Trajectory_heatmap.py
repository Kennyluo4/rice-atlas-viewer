import streamlit as st
# import matplotlib.pyplot as plt
import subprocess
import os, subprocess, glob, atexit
        
###############################################
##           Configuration                  ###
###############################################

## define function for generation figures using R
@st.cache_data
def generate_png_plot(species, feature, trajectory,candidateID_file):
    """
    Run the R script to generate the PNG plot.
    """
    import hashlib


    out_path = 'Rplot/temp'

    cmd = ["Rscript", "Rplot/plot_trajectory_genes_TFs_motifs_ACRs_100824.R", str(species), str(feature),str(trajectory), str(candidateID_file), str(out_path)]
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
    
    ## set the plot setting using form submission: 
    with st.sidebar.form(key='Plot_config'):
        species = st.selectbox('Species', [ '---Please choose---','Osativa', 'test2'])
        feature = st.selectbox('Feature', [ '---Please choose---','Gene', 'TF', 'Mt', 'ACR'])
        trajectory = st.selectbox('Trajectory',['---Please choose---','bud_ProcamPhloem', 'bud_ProcamXylem', 'Eseedling_SAMPhloem', 'Eseedling_SAMPhloemCC', 'Eseedling_SAMXylemPrecursor', 'crownroot_QCCortex', 'crownroot_QCMetaXylem', 'crownroot_QCProtoXylem','panicle_PBMSBMSMFM','panicle_PBMSMFM', 'semroot_QCCortex', 'semroot_QCMetaXylem', 'semroot_QCProtoXylem'])
        feature_list = st.text_area('candidate feature ID (e.g. gene ID, motif ID, and ACRs)', 
                                         'NAC011\nNAC005\nMYB46\nWRKY65')
        # plot_h = st.slider('plot height(cm)', 1,30,10)
        # plot_w = st.slider('plot width(cm)', 1,30,10)
        
        submitted = st.form_submit_button("Submit")
        
        # Clear cache when inputs change
        if submitted:
            st.cache_data.clear()  # Clear cache to avoid reusing old data
            st.write('Generating heatmap...')
            # st.write(f"Plot genes for {gene_list}")
    
    ## generate input files required by Rscript
    if species != '---Please choose---' and feature != '---Please choose---' and trajectory != '---Please choose---':
        # candidate_feature_file = f'ipt_target_species <- "{species}"\nipt_target_organ <- "{tissue}"\nplot_width <- {plot_w}\nplot_height <- {plot_h}\nFDR <- {fdr}'   
        candidate_feature_file = f'Rplot/temp/input_3_users_provide_ID_fl_{feature}.txt'
        with open(candidate_feature_file, 'w') as handle:
            handle.write(feature_list)
    # st.write('Updated plot configuration')
    
    ###############################################
    ##                Main page                 ###
    ###############################################
    st.markdown('### Plot the trajectory')
    if species != '---Please choose---' and feature != '---Please choose---' and trajectory != '---Please choose---':
        fig = generate_png_plot(species, feature, trajectory, candidate_feature_file)
        st.image(fig)
    else:
        st.error(':point_left: Please select the dataset and genes for plotting')
    
                             
if __name__ == '__main__':
    main()


