<div class="x-hidden" id="resultsinfo">
% if run!=None:
<p id="description">${run.description}</p>
<fieldset class="x-fieldset x-fieldset-default">
<legend class="x-fieldset-header x-fieldset-header-default">Generate metabolite options</legend>
Maximum number of reaction steps: ${run.n_reaction_steps}<br/>
Metabolism types: ${run.metabolism_types}<br/>
</fieldset>
<fieldset class="x-fieldset x-fieldset-default">
<legend class="x-fieldset-header x-fieldset-header-default">MS data options</legend>
MS Filename: ${run.ms_filename}<br/>
Maximum MS level: ${run.max_ms_level}<br/>
Precision for matching precursor mz with peak mz in parent scan: ${run.precursor_mz_precision}<br/>
Absolute intensity threshold for storing peaks in database: ${run.abs_peak_cutoff}<br/>
Fraction of basepeak intensity threshold threshold for storing peaks in database: ${run.rel_peak_cutoff}<br/>
</fieldset>
<fieldset class="x-fieldset x-fieldset-default">
<legend class="x-fieldset-header x-fieldset-header-default">Annotate options</legend>
Skip fragmentation:
% if run.skip_fragmentation:
Yes
% else:
No
% endif
<br/>
Absolute intensity minimum of lvl1 scan peaks which are matched with metabolites: ${run.ms_intensity_cutoff}<br/>
Ratio of basepeak intensity: ${run.msms_intensity_cutoff}<br/>
M/z offset which is allowed for matching a metabolite mass to m/z of a peak: ${run.mz_precision}<br/>
Maximum number of bonds broken in substructures generated from metabolites: ${run.max_broken_bonds}<br/>
Ionisation mode:
% if run.ionisation_mode == 1:
Positive
% else:
Negative
% endif
<br/>
Annotate all peaks, also peaks without fragmentation data:
% if run.use_all_peaks:
Yes
% else:
No
% endif
<br/>
</fieldset>
% else:
<p id="description">No information, nothing has been done</p>
% endif
</div>