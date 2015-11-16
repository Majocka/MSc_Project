
## msc_project

The main aim was to extract evolutionary knowledge from protein structures via a structural alignment tool and in-house 2-stage hierarchical clustering. This was to obtain structural regions of protein that have remained preserved throughout the evolution. An example of such color-coded regions is shown below for protein 1f6m_E. 

<img src="https://github.com/Majocka/msc_project/blob/master/results_to_view/image_final_clusters.png" width="300" height="270" />

To reproduce the results for 1f6m_E after you clone this repository do the following:
* run *main_chimera.py* from within Chimera (https://www.cgl.ucsf.edu/chimera/). You will get two prompts:
      * to run *first_clustering.py* from normal Python terminal
      * to adjust the view of protein in Chimera, in which you want the images of first clusters to be stored
* run *second_clustering.py* again from within Chimera, where you will be prompted to adjust the view and save the final image with the color-coded regions shown above.


For a reason of this project being ongoing I am not contemplating on details here.
