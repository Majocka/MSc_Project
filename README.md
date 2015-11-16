
## msc_project

The main aim was to extract evolutionary knowledge from protein structures via a structural alignment tool and in-house 2-stage hierarchical clustering so we can obtain the structural regions of protein that has been preserved throughout the evolution. An example of such color-coded regions can be seen below for protein 1f6m_E. 

<img src="https://github.com/Majocka/msc_project/blob/master/results_to_view/image_final_clusters.png" width="300" height="270" />

To reproduce the results for 1f6m_E do the following:
* run *main_chimera.py* from within Chimera (https://www.cgl.ucsf.edu/chimera/). You will get two prompts:
      * to run *first_clustering.py* from normal Python terminal
      * to adjust the view of protein in Chimera, in which you want the images of the first clusters to be stored in
* run *second_clustering.py* again from within Chimera, where you will be prompted to adjust the view and save the final image wiht color-coded regions as shown above.
