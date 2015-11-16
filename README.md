
## msc_project

The aim was to extract evolutionary knowledge from protein structures via PDBeFold structural alignment tool and 2-stage hierarchical clustering. 

BACKGROUND:
The main goal of structural biology is to observe protein structures at atomic levels to understand their basic functions and ultimately the basic life processes at molecular and cellular level. Single technique to obtain the structures, however, does not exist. On one hand, we have cryo-EM technique that allows for determination of protein structures in various functional states but at too low resolution and on the other hand crystallography method that does not allow for functional information but acquires structures of sufficiently high-resolution. With flexible fitting it is possible to computationally bridge the gap between the two techniques by fitting the crystallography structures into EM-density maps obtained via cryo-EM. 

Prior to fitting identification of protein rigid bodies is essential to prevent the fitting process from getting trapped in many local minima and was the main aim of my MSc project. Rigid bodies are regions at domain or subdomain level, which during the fitting stay rigid as opposed to the flexible regions, which are repositioned during the search for optimal conformation of protein. The rigid bodies, here called 'final clusters' are predicted by the use of ‘evolutionary’ information, which is extracted from protein structures firstly by use of structural alignment tool PDBeFold. The tool’s output is an xml file with all necessary information about the target protein (PDB code 1f6m_E) and all its structural neighbors identified by the tool referred to as hits. The code here takes the xml file along with additional txt file that contains further information on the target structure and applies 2-stage hierarchical clustering to eventually produce the rigid bodies. In the first stage the clustering identifies local matches of every hit separately, whilst in the second stage it pools all hits together to produce information on rigid bodies. For example, two color coded regions of 1f6m_E are shown below.

![alt tag](https://github.com/Majocka/msc_project/blob/master/results_to_view/image_final_clusters.png)
