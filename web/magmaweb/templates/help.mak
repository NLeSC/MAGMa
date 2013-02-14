<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>MAGMa - Help</title>
</head>
<body>
<h1>MAGMa Help page</h1>
<h2>Home page</h2>
Currently only contains the button to go to the start page.
<h2>Start page</h2>
The start page enables to run a new annotation. This requires three types of input:
<h3>Molecules</h3>
Molecule can be imported by:
<ul>
<li>Retrieval from a Database (PubChem, Kegg of HMDB (Human Metabolite Database))</li>
<li>Upload from smiles or SDF (Structure Data File format)</li>
<li>Drawing a molecule</li>
</ul>
<h3>MS Data</h3>
MS Data can be imported, either as file or via the text field, as:
<ul>
<li>mzXML: a public dataformat that can be generated with
<ul>
<li><a href="http://sourceforge.net/projects/sashimi/files/ReAdW%20%28Xcalibur%20converter%29/">READW</a></li>
<li><a href="http://proteowizard.sourceforge.net/downloads.shtml">MSConvert</a> module from ProteoWizard</li>
</ul>
</li>
<li>Tree: a custom format to enter a single spectral tree, e.g.<pre>
MS1_peak_mz:MS1_peak_intensity (
MS2_peak1_mz: MS2_peak1_intenstity,
MS2_peak2_mz: MS2_peak2_intenstity (
    MS3_peak1_mz: MS3_peak1_intenstity,
    MS3_peak2_mz: MS3_peak2_intenstity
),
…
MS2_peakX_mz: MS2_peakX_intenstity
)
</pre></li>
</ul>
During import the (mzXML) data can be filtered on:
<ul>
<li>MS level: enables to ignore deeper levels in a spectrum</li>
<li>Noise:  peaks below this intensity threshold are ignored</li>
</ul>

<h3>Annotate options</h3>
<ul>
<li>Ionisation mode</li>
<li>Substructure options
<ul>
<li>Bond breaks: maximum number of bonds to break when generating substructures</li>
<li>Additional water losses: maximum additional OH groups to remove from substructures</li>
</ul>
</li>
<li>Precision:<ul>
<li>Relative (ppm) and Absolute (Da) mass tolerance for matching (sub)structures with m/z values. (The maximum of these two is used)</li>
<li>Precursor m/z (Da): This is used to link precursor ion m/z values with the corresponding peaks in the parent spectrum</li>
</ul></li>
<li>Intensity thresholds:
<ul>
<li>MS1 (abs): can be used to only annotate the most intense peaks in the Chromatogram</li>
<li>MS2n>1 (% of base peak): threshold to filter fragment ion peaks</li>
</ul>
</li>
</ul>
Start the calculation with the Submit button on the lower right of the start page. When the calculation has finished it will automatically open the results page. If the browser is closed during the calculation, it will continue running. You will be able to access the results with the workspace page (see below).

<h2>Results page</h2>

The results page consists of 4 panels:
<h3>Molecules table (upper left)</h3>
<ul>
<li>Contains the candidate molecules, imported or retrieved from the structure database
</li><li>The Reference column provides links to original database entries
</li><li>You can sort on a column by clicking the header
</li><li>You can filter, hide and show columns via the drop-down menu that is available at the column header
</li><li>You can change order and width of the columns by drag and drop
</li><li>When you have imported LC-MSn data:
<ul><li>The Scans column indicates the number of spectral trees/MS1 scans containing matching (precursor) ions
</li><li>The Assigned column indicates if you have assigned this molecule to a (precursor) ion in your dataset (the “assign molecule” button is explained below)
</ul>
</li><li>By clicking a row of the table, you can select a candidate molecule, in case of LC-MSn data:
<ul><li>The Chromatogram panel will show an extracted ion chromatogram for the corresponding m/z
</li><li>The spectral trees matching the candidate molecule are labeled with green triangles in the Chromatogram panel
</li></ul>
</li><li>The buttons in the upper right corner can be used to export the metabolites and to clear the table filters</li>
</ul>
<h3>Chromatogram (bottom left)</h3>
<ul>
<li>Appears only when LC-MSn data has been imported
</li><li>A spectral tree can be selected by clicking one of the green triangles, as a result:
<ul><li>the Molecules table will show only candidate molecules that match the selected spectral tree
</li><li>two new columns appear in the Molecules table:
<ul><li>
Candidate (penalty) score: lower values indicate more likely candidates</li><li>
ΔMass indicates the deviation between the calculated and the observed m/z value</li></ul>
</ul></li><li>Zooming in/out can be performed using the mouse-wheel
</li><li>The buttons in the upper right corner can perform the following tasks
<ul>Select a spectral tree by scan number
</li><li>Clear the scan selection
</li><li>Reset the zoom
</li><li>Turn on/off zooming in X and Y direction
</ul>
</li></ul>
<i>When you have selected both a molecule and a spectral tree, the annotated spectra become available on the right.</i>
<h3>Substructures (upper right)</h3>

<ul><li>This table has a tree structure: the folder icons allow you to open up, or close, the next level of substructures. At the same time, the corresponding spectrum is shown, or hidden, at the bottom.
</li><li>Selecting a substructure highlights the corresponding peak in the spectrum below</li>
</ul>
<h3>MS spectral tree (bottom right)</h3>
<ul>
<li>Selecting a peak by clicking one of the green triangles highlights the corresponding substructure.
</li><li>You can zoom in using the mouse-wheel
</li>
</ul>

<h3>“Assign molecule to peak” button</h3>
With this button a selected candidate molecule can be assigned to an MS1 peak. The assignment can be undone by clicking the button again. The corresponding scan is highlighted in the Chromatochram panel. (To give some overview of all assignments), and the Assigned column in the Molecules table is updated.

<h2>Workspace</h2>
The workspace page gives a list of annotations that have performed. Old data can be removed by clicking the minus button in the right column. Double clicking the Description field allows typing a brief description of the purpose of the annotation.
</body>
</html>