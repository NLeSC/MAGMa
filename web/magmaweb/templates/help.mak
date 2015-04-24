<!DOCTYPE html>
<html>

<head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>MAGMa - Help</title>
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.4/css/bootstrap.min.css">
</head>

<body>
    <div class="container">
        <div class="row">
            <div class="span12">
                <h1>MAGMa Help page</h1>
                <h2 id="home">Home page</h2>
                The home page enables to run a new MAGMa job. This requires three types of input:
                <h3 id="msdata">MS Data</h3>
                MS Data can be imported, either as file or via the text field, as:
                <ul>
                    <li>mzXML: a public dataformat that can be generated with
                        <ul>
                            <li><a href="http://sourceforge.net/projects/sashimi/files/ReAdW%20%28Xcalibur%20converter%29/">READW</a></li>
                            <li><a href="http://proteowizard.sourceforge.net/downloads.shtml">MSConvert</a> module from ProteoWizard</li>
                        </ul>
                    </li>
                    <li>Mass Tree: a custom format to enter a single spectral tree based on m/z, e.g.
                        <pre>MS1_peak_mz:MS1_peak_intensity (
    MS2_peak1_mz: MS2_peak1_intenstity,
    MS2_peak2_mz: MS2_peak2_intenstity (
        MS3_peak1_mz: MS3_peak1_intenstity,
        MS3_peak2_mz: MS3_peak2_intenstity
    ),
    …
    MS2_peakX_mz: MS2_peakX_intenstity
)</pre>
                    </li>
                    <li>Formula Tree: a custom format to enter a single spectral tree based on molecular formulas, e.g.
                        <pre>MS1_peak_formula:MS1_peak_intensity (
    MS2_peak1_formula: MS2_peak1_intenstity,
    MS2_peak2_formula: MS2_peak2_intenstity (
        MS3_peak1_formula: MS3_peak1_intenstity,
        MS3_peak2_formula: MS3_peak2_intenstity
    ),
    …
    MS2_peakX_formula: MS2_peakX_intenstity
)</pre>
                    </li>
                </ul>
                During import the mzXML data can be filtered on:
                <ul>
                    <li>MS level: enables to ignore deeper levels in a spectrum</li>
                    <li>Noise: peaks below this intensity threshold are ignored</li>
                    <li>MS1 scan number: Only read spectral tree specified by MS1 scan number.
                        <br> Note: the public website requires specification of one scan if candidates are retrieved from a database.</li>
                </ul>
                When importing mzXML data, please:
                <br>
                <ul>
                    <li><b>Use centroid data!</b> (MAGMa considers all signals to be different ions.)
                    </li>
                    <li>Note that when multiple scans are provided for the same precursor (in the same precursor scan) MAGMa merges these scans into a single composite spectrum which includes only the most intense peak for a given m/z, using the given precision
                        to match peaks from the different scans.</li>
                </ul>
                <ul>
                </ul>
                <h3 id="molecules">Molecules</h3>
                Molecules can be imported by:
                <ul>
                    <li>Retrieval from a Database (PubChem, Kegg or Human Metabolite Database). Precursor ion m/z values of all MS2 scans are used to query the database on mono-isotopic mass, assuming [M-H]<sup>-</sup> or [M+H]<sup>+</sup> ions, depending
                        on the specified ionisation mode.</li>
                    <li>Upload from smiles or SDF (Structure Data File format).</li>
                    <li>Drawing a molecule.</li>
                </ul>
                Metabolizing molecules:
                <br>
                <br> Uploaded or drawn molecules can be submitted to <i>in silico</i> metabolite prediction to augment the set of candidate molecules. MAGMa metabolite prediction is based on reaction rules that are applied according to a specified editable
                "scenario". The buttons "Drugs" and "Polyphenols" provide predefined scenario's that can be used to generate human metabolites for these classes of compounds. Phase1 and phase2 human biotransformation rules are derived according to the
                method outlined in [<a href="http://dx.doi.org/10.1002/cmdc.200700312">Ridder
      and Wagener 2008</a>]. Phase1_selected and phase2_selected are subsets including the most common rules: more than 5% of the metabolites they predict for a large and diverse set of clinically studied compounds agreed with the clinical observations.
                The "glycosidase" and "gut" rules are described in [Ridder et al. submitted]. The first set of rules defined in the scenario is applied to all uploaded or drawn molecule(s). Subsequent stages are also applied to products of the preceding
                steps. When a mass filter is applied, the subsequent stage in the scenario are performed only on molecules below the specified mass limit. The value in the "Steps" column determines either the number of reaction steps to be applied with
                the specified set of rules, or the mass limit in case of a mass filter step. The glycosidase rule can also be used with the word "complete", which will iterate the glycosidase rule until no further metabolites are generated. In this case
                only the end products are used for further processing.
                <br>
                <ul>
                </ul>
                <h3 id="settings">Parameters settings</h3>
                <ul>
                </ul>
                <ul>
                    <li>Ionisation mode</li>
                    <li>Substructure options
                        <ul>
                            <li>Bond dissociations: maximum number of bonds to break when generating substructures</li>
                            <li>Additional small losses: maximum number of additional water (OH) and/or ammonia (NH2) losses from substructures</li>
                        </ul>
                    </li>
                    <li>Accuracy:
                        <ul>
                            <li>Relative (ppm) and Absolute (Da) mass tolerance for matching (sub)structures with m/z values. (The maximum of these two is used.)</li>
                            <li>Precursor m/z (Da): This is used to match precursor ion m/z values with the corresponding peaks in the parent spectrum</li>
                        </ul>
                    </li>
                    <li>Intensity thresholds:
                        <ul>
                            <li>MS<sup>1</sup> (abs): can be used to only annotate the most intense peaks in the Chromatogram</li>
                            <li>MS<sup>n&gt;1</sup> (% of base peak): threshold to filter fragment ion peaks</li>
                        </ul>
                    </li>
                </ul>
                For more information on the parameters used in the automatic annotation please refer to [<a href="http://dx.doi.org/10.1002/rcm.6364">Ridder et al.
      2012</a>].
                <br>
                <br>
                <br>
                <i><b>Start the calculation with the Submit button on the lower right of the
        start page. When the calculation has finished it will automatically open
        the <a href="#results">results page</a>. If the browser is closed
        during the calculation, it will continue running. You will be able to
        access the results from the <a href="#workspace">workspace</a>. When
        the MAGMa processing fails an error message will appear. You can check
        the log file via the Information panel (click "Information" on the
        top-left button group, then click the "Save" icon on the top right of
        the Information panel).</b></i>
                <br>
                <br>
                <h2 id="results">Results page</h2>
                <img width="100%" src="${request.static_url('magmaweb:static/img/metabolites.png')}" /> The results page consists of 4 panels:
                <h3 id="molpanel">Molecules table (upper left)</h3>
                <ul>
                    <li>Contains the candidate molecules, imported or retrieved from the structure database </li>
                    <li>The Reference column provides links to original database entries </li>
                    <li>You can sort on a column by clicking the header </li>
                    <li>You can filter, hide and show columns via the drop-down menu that is available at the column header. When a filter is applied to the molecules table, the green triangles in the Chromatogram panel are filtered accordingly, to show only
                        matching precursor peaks, i.e. the corresponding precursor MS1 scans.</li>
                    <li>You can change order and width of the columns by drag and drop </li>
                    <li>When you have imported LC-MSn data:
                        <ul>
                            <li>The Scans column indicates the number of MS2 scans with matching precursor m/z values in the data</li>
                            <li>The Assigned column indicates if you have assigned this molecule to a precursor ion in your dataset (the “assign molecule” button is explained below) </li>
                        </ul>
                    </li>
                    <li>By clicking a row of the table, you can select a candidate molecule. In case of LC-MSn data:
                        <ul>
                            <li>The Chromatogram panel will show an extracted ion chromatogram for the corresponding m/z </li>
                            <li>The precursor peaks (i.e. the corresponding precursor scans) matching the candidate molecule are labeled with green triangles in the Chromatogram panel </li>
                        </ul>
                    </li>
                    <li>When <i>in silico</i> metabolite prediction has been applied, the Reactions column provides information about the number of connected reactant and product molecules. Numbers in bold/italic indicate that these include metabolites which
                        match one or more precursor peaks in the data. When clicking this information the molecules are filtered to display the actual reactant or product structures. When using this reactant or product filter, it is recommended to turn
                        off other filters on the molecules table.</li>
                    <li>The buttons in the upper right corner can be used to export the displayed molecules and to clear the table filters.</li>
                </ul>
                <h3 id="chromatogram">Chromatogram (bottom left)</h3>
                <ul>
                    <li>Appears only when LC-MSn data has been imported </li>
                    <li>A spectral tree can be selected by clicking one of the green triangles, as a result:
                        <ul>
                            <li>the Molecules table will show only candidate molecules that match the selected spectral tree </li>
                            <li>two new columns appear in the Molecules table:
                                <ul>
                                    <li> Candidate (penalty) score: lower values indicate more likely candidates
                                    </li>
                                    <li> ΔMass indicates the deviation between the calculated and the observed m/z value</li>
                                </ul>
                            </li>
                        </ul>
                    </li>
                    <li>Zooming in/out can be performed using the mouse-wheel </li>
                    <li>The buttons in the upper right corner can perform the following tasks
                        <ul>
                            <li>Select a spectral tree by scan number</li>
                            <li>Clear the scan selection</li>
                            <li>Reset the zoom</li>
                            <li>Turn on/off zooming in X and Y direction</li>
                        </ul>
                    </li>
                </ul>
                <i><b>Only when you have selected both a molecule and a precursor, the
        annotated MSn spectra become available on the right.</b></i>
                <h3 id="substructures">Substructures (upper right)</h3>
                <ul>
                    <li>This table has a tree structure: the folder icons allow you to open up, or close, the next level of substructures. At the same time, the corresponding spectrum is shown, or hidden, at the bottom. </li>
                    <li>Selecting a substructure highlights the corresponding peak in the spectrum below.</li>
                </ul>
                <h3 id="scan">MS scans (bottom right)</h3>
                <ul>
                    <li>Selecting a peak by clicking one of the green triangles highlights the corresponding substructure. </li>
                    <li>You can zoom in using the mouse-wheel.</li>
                </ul>
                <h3>“Assign molecule to peak” button</h3>
                With this button a selected candidate molecule can be assigned to an MS1 peak. The assignment can be undone by clicking the button again. The corresponding scan is highlighted in the Chromatochram panel. (To provide an overview of all assignments), and
                the Assigned column in the Molecules table is updated.
                <h2>The top-left button group</h2>
                <p>The button group on the top-left of the screen provides links to different components of MAGMa.</p>
                <ul>
                    <li>Home: go to the MAGMa <a href="#home">home page</a> to start a new job
                    </li>
                    <li>Help: to to this help page</li>
                    <li>Workspace: go to the <a href="#workspace">workspace</a> to manage MAGMa annotated datasets</li>
                    <li>Information: opens an information panel with all parameters settings used for the current dataset.</li>
                    <ul>
                        <li>The "save" icon on the top-right of this panel will open the log file generated during the MAGMa processing.</li>
                    </ul>
                </ul>
                <h2 id="workspace">Workspace</h2>
                The workspace page gives an overview of (and links to) automatic annotations that have performed. Old data can be removed by clicking the minus button in the right column. Fields of the job list can be edited by clicking them:
                <ul>
                    <li>Description, allows typing a brief description for annotation purposes.
                    </li>
                    <li>MS filename, allows typing a file name of the MS data input file for annotation purposes.</li>
                    <li>Public? By default jobs are private, i.e. they can only be seen by the person who submitted it. The job can be toggled to Public. A public job can be seen by another authenticated user (or anyone in case of the public website) knowing
                        the job url.</li>
                </ul>
            </div>
        </div>
    </div>
    <footer class="footer">
        <div class="container">
            <a href="${request.route_url('home')}">Home</a></div>
    </footer>
</body>

</html>
