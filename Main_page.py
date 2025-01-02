import streamlit as st
import os, atexit, glob
import streamlit as st


###############################################
##             Configuration                ###
###############################################
# Set the page config to wide mode
st.set_page_config(layout='centered',
                   page_title='Plant atlas-main page',
                   page_icon=':ear_of_rice:'
                    )

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
    ##                Main page                 ###
    ###############################################
    
    st.markdown("""
                <style>
                .big-font {
                    font-size:12px !important;
                }
                </style>
                """, unsafe_allow_html=True)
    
    st.markdown('# Welcome to the plant single-cell atlas of cis-regulatory elements ')
    st.image('image/rice_atlas.png', width=600)
    st.markdown('### Visualization tools: \n')
    
    st.page_link('pages/Dotplot.py', label='Gene dotplot', help='Plot the gene cell-type-specific expression level using dotplot')
    st.page_link('pages/UMAP.py', label='Feature UMAP plot', help='Plot the gene expression in the single-cell UMAP')
    st.page_link('pages/Motif_heatmap.py', label='Motif heatmap plot', help='Plot the heatmap of selected motifs')
    st.page_link('pages/Trajectory_heatmap.py', label='Trajectory plot', help='Plot the trajectory heatmap of selected motifs')
    
    
    
    ## html for footnote
    footnote = """
                <div class="footer">
                    <p style='font-size: small;'>For issues and questions, please contact 
                    <a href='mailto:schmitz@uga.edu'>Schmitz_lab</a>.</p>
                </div>
                """
    st.markdown(footnote, unsafe_allow_html=True)
    
    
if __name__ == '__main__':
    main()
