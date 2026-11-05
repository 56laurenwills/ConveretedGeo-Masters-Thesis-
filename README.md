# Converted Geometry Masters-Thesis-
Repository for scripts and files a part of Converted Geometry Master's Thesis

Links to other repositories: 
- Local ACTS build needed to create Augmented ACTS material maps [Local ACTS Build Pin 47.7.0](https://github.com/56laurenwills/acts/commit/e2d10875e8dcdc2250072490482d61f87d6bffd2)
- Person Athena fork for augmented material maps in Athena [Personal Athena Fork 25.0.73](https://gitlab.cern.ch/lwills/athena/-/commit/3cdf226d71a48dc8557e1c4bdd98c7e64fc14167)

<br>
This repository includes folders: <br>
<br>
<strong> Analysis: </strong> Files regarding the ACTS material map analysis 

- The least detailed script: 
``` python MM_json_summary_info.py -f ../path/to/material-map.json``` Prints a quick summary of the total number of surfaces, and their types. This is used to make sure the create material map has the expected surfaces included with the correct amount of null material (bins without any material assignment). The output looks like: 
    ```Total Number of Surfaces: 208, Type of Surface: N Cylindrical: 22, Number of Disc: 186
    Number of Entries with V/L/A: 185, with V/L: 1, with V/B: 22
    Number of Entries with mapMaterial = False : None vs. Entries with Null Material: 0
    Total number of bins: 522035
    ```
 
- The more detailed script: `python MM_json_information.py --material_map ../path/to/material-map.json -o acts-itk_mm_info.txt` Creates a txt files containing the average material values for every surface in the material map. This is particularly useful to confirm that the 'thickness' values of the material maps are physical and can be used in the conversion script. 
  - The output of the two JSON material maps (the internally released ITk group material map and the locally create material map) are under `acts-itk_mm_info.txt` and `local_mm_info.txt` respectively. 
- There is also a filtering script, used to remove any surfaces that you no longer want in your material map. This script will hopefully not have to be used in any future mapping.
- The JSON and ROOT jupyter notebooks are generic, open-ended analysis notebooks, used to study the produced material maps in a more detail.
     - The JSON jupyter notebook includes a summary of the comparison between a locally produced v.44.0 material map and the ITk-group material map, included to show the initial complications of material map production in standalone ACTS.

<br>

<strong> Carbon Toy Model: </strong> 
the macro and resulting 1000bin Carbon-toy model of the ITkPixel sub-detector.
   -  The macro is to be used with an augmented material track ROOT file (created using the linked local ACTS build) to convert the original material track recording to re-write every observed material to be only carbon (all other information from the Geantino scan is preserved)  
  <br>

<strong> Conversion Tool: </strong> A 'sanetized' version of the Conversion tool (all specific geometry information, such as the exact the end-cap z-value text file has been removed)
   - The input is an augmented ACTS material map created with the linked local ACTS build 
- the default it a GeoModel GMX output, which is to be called within ITKLayouts ITkPixel and ITkStrip.gmx. Modifications to `XML_production_assembly.py` must be made if you wish to use the output alone (within GeoModel). 


To create a geometry, run: 

``` python ACTS_to_G4.py --material /path/to/DRAFT-Run4-05-00_material-maps.root --regions Strip --output ConvertedStrip --geotag "Run4-05-00-00"```

- Then copy the created geometry to the Simplified ITkStrip pathway
- For best results, create Pixel and Strip separately (this will be fixed in the future)


<br>
<strong> Geometries: </strong> The preliminary converted ITkPixel and ITkStrip geometries with the original sensors included, converted into `.db` files. 
