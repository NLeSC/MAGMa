"""Module for job query"""
import json
import os
from cgi import FieldStorage
import colander


class JobQuery(object):
    """Perform actions on a :class:`Job` by parsing post params,
    staging files and building magma cli commands
    """

    class File(object):
        """ Colander schema type for file upload field"""
        def serialize(self, node, appstruct):
            return appstruct

        def deserialize(self, node, cstruct):
            if cstruct is colander.null:
                return colander.null
            if cstruct == '':
                return colander.null
            if not isinstance(cstruct, FieldStorage):
                raise colander.Invalid(
                    node,
                    '{} is not a cgi.FieldStorage'.format(cstruct)
                )
            return cstruct

    def __init__(self,
                 directory,
                 script='',
                 prestaged=None,
                 status_callback_url=None,
                 restricted=False,
                 ):
        """Contruct JobQuery

        Params:
        - directory, job directory
        - script, script to run in job directory
        - prestaged, list of files to prestage
        - status_callback_url, Url to PUT status of job to.

        """
        self.dir = directory
        self.script = script
        self.prestaged = prestaged or []
        self.status_callback_url = status_callback_url
        self.restricted = restricted

    def __eq__(self, other):
        return (self.dir == other.dir and
                self.script == other.script and
                self.prestaged == other.prestaged and
                self.status_callback_url == other.status_callback_url and
                self.restricted == other.restricted
                )

    def __repr__(self):
        """Return a printable representation."""
        s = "JobQuery({!r}, script={!r}, "
        s += "prestaged={!r}, status_callback_url={!r},"
        s += "restricted={!r})"
        return s.format(self.dir, self.script,
                        self.prestaged, self.status_callback_url,
                        self.restricted,
                        )

    def escape(self, string):
        """ Replaces single quote with its html escape sequence"""
        return str(string).replace("'", '&#39;')

    def _addAnnotateSchema(self, schema):
        schema.add(colander.SchemaNode(colander.Float(),
                                       missing=0.0,
                                       validator=colander.Range(0, 1),
                                       name='precursor_mz_precision'))
        schema.add(colander.SchemaNode(colander.Float(),
                                       missing=0.0,
                                       validator=colander.Range(0, 1000),
                                       name='mz_precision'))
        schema.add(colander.SchemaNode(colander.Float(),
                                       missing=0.0,
                                       validator=colander.Range(0, 1),
                                       name='mz_precision_abs'))
        schema.add(colander.SchemaNode(colander.Float(),
                                       missing=0.0,
                                       name='ms_intensity_cutoff'))
        schema.add(colander.SchemaNode(colander.Float(),
                                       missing=0.0,
                                       name='msms_intensity_cutoff'))
        schema.add(colander.SchemaNode(colander.Integer(),
                                       validator=colander.OneOf([-1, 1]),
                                       name='ionisation_mode'))
        schema.add(colander.SchemaNode(colander.Integer(),
                                       validator=colander.Range(0, 4),
                                       name='max_broken_bonds'))
        schema.add(colander.SchemaNode(colander.Integer(),
                                       validator=colander.Range(0, 4),
                                       name='max_water_losses'))
        schema.add(colander.SchemaNode(colander.String(),
                                       missing=colander.null,
                                       validator=colander.OneOf(['pubchem',
                                                                 'kegg',
                                                                 'hmdb',
                                                                 ]),
                                       name='structure_database'
                                       ))
        schema.add(colander.SchemaNode(colander.Integer(),
                                       missing=colander.null,
                                       validator=colander.Range(min=1),
                                       name='min_refscore'))
        schema.add(colander.SchemaNode(colander.Integer(),
                                       missing=colander.null,
                                       validator=colander.Range(min=1),
                                       name='max_mz'))

    def _addMetabolizeSchema(self, schema):
        scenario = colander.SchemaNode(colander.Mapping())
        transformation_types = colander.OneOf([
                                               'phase1',
                                               'phase2',
                                               'phase2_selected',
                                               'glycosidase',
                                               'mass_filter',
                                               'gut'])

        scenario.add(colander.SchemaNode(colander.String(),
                                         name='type',
                                         validator=transformation_types))
        # TODO add validator for steps, but what are valid options for steps?
        scenario.add(colander.SchemaNode(colander.String(),
                                         name='steps'))

        schema.add(colander.SchemaNode(colander.Sequence(),
                                       scenario,
                                       validator=colander.Length(1),
                                       name='scenario'))

    def _addAddStructuresSchema(self, has_ms_data, must_metabolize):
        def textarea_or_file(node, value):
            """Validator that either textarea or file upload is filled"""
            if not(not value['structures'] is colander.null or
                   not value['structures_file'] is colander.null):
                error = 'Either structures or structures_file must be set'
                exception = colander.Invalid(node)
                exception.add(colander.Invalid(structuresSchema, error))
                exception.add(colander.Invalid(structures_fileSchema, error))
                raise exception

        schema = colander.SchemaNode(colander.Mapping(),
                                     validator=textarea_or_file)
        formats = ['smiles', 'sdf']
        schema.add(colander.SchemaNode(colander.String(),
                                       validator=colander.OneOf(formats),
                                       name='structure_format'))
        filled = colander.Length(min=1)
        structuresSchema = colander.SchemaNode(colander.String(),
                                               validator=filled,
                                               missing=colander.null,
                                               name='structures')
        schema.add(structuresSchema)
        structures_fileSchema = colander.SchemaNode(self.File(),
                                                    missing=colander.null,
                                                    name='structures_file')
        schema.add(structures_fileSchema)
        schema.add(colander.SchemaNode(colander.Boolean(),
                                       default=False,
                                       missing=False,
                                       name='metabolize'))

        if has_ms_data:
            self._addAnnotateSchema(schema)
        if must_metabolize:
            self._addMetabolizeSchema(schema)

        return schema

    def add_structures(self, params, has_ms_data=False):
        """Configure job query to add_structures from params.

        ``params`` is a MultiDict from which the following keys are used:

        * structure_format, in which format is structures or structure_file
        * structures, string with structures
        * structures_file, file-like object with structures
        * metabolize, when key is set then
            :meth:`~magmaweb.job.JobQuery.metabolize` will be called.

        If ``has_ms_data`` is True then
            :meth:`~magmaweb.job.JobQuery.annotate` will be called.
        If both ``stuctures`` and ``structures_file`` is filled then
            ``structures`` is ignored.
        """
        must_metabolize = 'metabolize' in params
        schema = self._addAddStructuresSchema(has_ms_data, must_metabolize)

        if hasattr(params, 'mixed'):
            # unflatten multidict
            params = params.mixed()
        self._deserialize_scenario(params)

        params = schema.deserialize(params)

        metsfile = file(os.path.join(self.dir, 'structures.dat'), 'w')
        if not params['structures_file'] is colander.null:
            sf = params['structures_file'].file
            sf.seek(0)
            while 1:
                data = sf.read(2 << 16)
                if not data:
                    break
                metsfile.write(data)
            sf.close()
        else:
            metsfile.write(params['structures'])
        metsfile.close()

        script = "{{magma}} add_structures -t '{structure_format}'"
        script += " structures.dat {{db}}"
        sf = self.escape(params['structure_format'])
        self.script += script.format(structure_format=sf)
        self.prestaged.append('structures.dat')

        if must_metabolize:
            self.script += " |"
            self.metabolize(params, has_ms_data, True)
        elif (has_ms_data):
            # dont annotate when metabolize is also done
            # as metabolize will call annotate
            self.script += " |"
            self.annotate(params, True)
        else:
            self.script += "\n"

        return self

    def _getMsDataSchema(self):
        def textarea_or_file(node, value):
            """Validator that either textarea or file upload is filled"""
            if not(not value['ms_data'] is colander.null or
                   not value['ms_data_file'] is colander.null):
                error = 'Either ms_data or ms_data_file must be set'
                exception = colander.Invalid(node)
                exception.add(colander.Invalid(msdata_stringSchema, error))
                exception.add(colander.Invalid(msdata_fileSchema, error))
                raise exception

        schema = colander.SchemaNode(colander.Mapping(),
                                     validator=textarea_or_file)

        valid_formats = colander.OneOf(['mzxml',
                                        'mass_tree',
                                        'form_tree',
                                        ])
        schema.add(colander.SchemaNode(colander.String(),
                                       validator=valid_formats,
                                       name='ms_data_format'
                                       ))

        filled = colander.Length(min=1)
        msdata_stringSchema = colander.SchemaNode(colander.String(),
                                                  validator=filled,
                                                  missing=colander.null,
                                                  name='ms_data')
        schema.add(msdata_stringSchema)

        msdata_fileSchema = colander.SchemaNode(self.File(),
                                                missing=colander.null,
                                                name='ms_data_file')
        schema.add(msdata_fileSchema)

        schema.add(colander.SchemaNode(colander.Integer(),
                                       missing=0,
                                       validator=colander.Range(min=0),
                                       name='max_ms_level'))
        schema.add(colander.SchemaNode(colander.Float(),
                                       missing=0.0,
                                       validator=colander.Range(min=0),
                                       name='abs_peak_cutoff'))
        schema.add(colander.SchemaNode(colander.Integer(),
                                       missing=colander.null,
                                       validator=colander.Range(min=0),
                                       name='scan'))
        return schema

    def _writeMsFile(self, params):
        msfile = file(os.path.join(self.dir, 'ms_data.dat'), 'w')
        if not params['ms_data_file'] is colander.null:
            msf = params['ms_data_file'].file
            msf.seek(0)
            while 1:
                data = msf.read(2 << 16)
                if not data:
                    break
                msfile.write(data)
            msf.close()
        else:
            msfile.write(params['ms_data'])
        msfile.close()

    def _writeScenarioFile(self, params):
        with open(os.path.join(self.dir, 'scenario.csv'), 'w') as f:
            for transformation in params['scenario']:
                f.write(",".join([transformation['type'],
                                  transformation['steps']
                                  ])+"\n")

    def _deserialize_scenario(self, params):
        if 'scenario' in params:
            try:
                params['scenario'] = json.loads(params['scenario'])
            except TypeError:
                pass

    def _addIonisationToFormulaTree(self, params, schema, orig_params):
        if params['ms_data_format'] == 'form_tree':
            if 'ionisation_mode' in orig_params:
                if orig_params['ionisation_mode'] == "1":
                    params['ms_data_format'] = 'form_tree_pos'
                elif orig_params['ionisation_mode'] == "-1":
                    params['ms_data_format'] = 'form_tree_neg'
            else:
                sd = schema['ms_data_format']
                msg = 'Require ionisation_mode when ms_data_format=form_tree'
                raise colander.Invalid(sd, msg)

    def add_ms_data(self, params, has_metabolites=False):
        """Configure job query to add ms data from params.

        ``params`` is a dict from which the following keys are used:

        * ms_data_format
        * ms_data, string with MS data
        * ms_data_file, file-like object with MS data
        * max_ms_level
        * abs_peak_cutoff

        If ``has_metabolites`` is True then
            :meth:`~magmaweb.job.JobQuery.annotate` will be called.

        If both ``ms_data`` and ``ms_data_file`` is filled then
            ``ms_data`` is ignored.
        """
        schema = self._getMsDataSchema()

        if has_metabolites:
            self._addAnnotateSchema(schema)
        orig_params = params
        params = schema.deserialize(params)

        self._writeMsFile(params)

        self._addIonisationToFormulaTree(params, schema, orig_params)

        script__substitution = {
            'ms_data_format': self.escape(params['ms_data_format']),
            'max_ms_level': self.escape(params['max_ms_level']),
            'abs_peak_cutoff': self.escape(params['abs_peak_cutoff']),
            'call_back_url': self.status_callback_url
        }
        script = "{{magma}} read_ms_data --ms_data_format '{ms_data_format}' "
        script += "-l '{max_ms_level}' "
        script += "-a '{abs_peak_cutoff}' "

        is_mzxml = params['ms_data_format'] == 'mzxml'
        empty_scan = params['scan'] is colander.null
        if is_mzxml and self.restricted and empty_scan:
            sd = schema['scan']
            msg = 'Require MS1 scan number'
            raise colander.Invalid(sd, msg)
        if not empty_scan and is_mzxml:
            script += "--scan '{scan}' "
            script__substitution['scan'] = self.escape(params['scan'])

        script += " --call_back_url '{call_back_url}' "
        script += "ms_data.dat {{db}}\n"
        self.script += script.format(**script__substitution)

        self.prestaged.append('ms_data.dat')

        if (has_metabolites):
            self.annotate(params)

        return self

    def metabolize(self, params, has_ms_data=False, from_subset=False):
        """Configure job query to metabolize all structures.

        ``params`` is a MultiDict from which the following keys are used:

        * n_reaction_steps
        * metabolism_types, comma seperated string with metabolism types

        If ``has_ms_data`` is True then
            :meth:`~magmaweb.job.JobQuery.annotate` will be called.

        If ``from_subset`` is True then metids are read from stdin
        """
        schema = colander.SchemaNode(colander.Mapping())
        self._addMetabolizeSchema(schema)
        if has_ms_data:
            self._addAnnotateSchema(schema)
        if hasattr(params, 'mixed'):
            # unflatten multidict
            params = params.mixed()

        self._deserialize_scenario(params)
        params = schema.deserialize(params)

        self._writeScenarioFile(params)
        self.prestaged.append('scenario.csv')

        script = "{{magma}} metabolize"
        script += " --scenario scenario.csv"
        script += " --call_back_url '{call_back_url}' "
        script_substitutions = {
            'call_back_url': self.status_callback_url,
        }
        self.script += script.format(**script_substitutions)

        if from_subset:
            self.script += " -j -"

        if (has_ms_data):
            self.script += " {db} |"
            self.annotate(params, True)
        else:
            self.script += " {db}\n"

        return self

    def metabolize_one(self, params, has_ms_data=False):
        """Configure job query to metabolize one structure.

        ``params`` is a MultiDict from which the following keys are used:

        * metid, metabolite identifier to metabolize
        * scenario, transformation scenario

        If ``has_ms_data`` is True then
            :meth:`~magmaweb.job.JobQuery.annotate` will be called.
        """
        schema = colander.SchemaNode(colander.Mapping())
        schema.add(colander.SchemaNode(colander.Integer(),
                                       validator=colander.Range(min=0),
                                       name='metid'))
        self._addMetabolizeSchema(schema)
        if has_ms_data:
            self._addAnnotateSchema(schema)
        if hasattr(params, 'mixed'):
            # unflatten multidict
            params = params.mixed()
        self._deserialize_scenario(params)
        params = schema.deserialize(params)

        self._writeScenarioFile(params)
        self.prestaged.append('scenario.csv')

        script = "echo '{metid}' | {{magma}} metabolize -j - "
        script += "--scenario scenario.csv {{db}}"
        script += " --call_back_url '{call_back_url}' "
        script_substitutions = {
            'metid': self.escape(params['metid']),
            'call_back_url': self.status_callback_url,
        }
        self.script += script.format(**script_substitutions)

        if (has_ms_data):
            self.script += " |"
            self.annotate(params, True)
        else:
            self.script += "\n"

        return self

    def annotate(self,
                 params,
                 from_subset=False,
                 ):
        """Configure job query to annotate.

        ``params`` is a dict from which the following keys are used:

        * precursor_mz_precision
        * mz_precision
        * mz_precsion_abs
        * ms_intensity_cutoff
        * msms_intensity_cutoff
        * ionisation_mode
        * max_broken_bonds
        * max_water_losses
        * structure_database,
            only used when ``structure_database_location`` is given
        * min_refscore, only used when ``structure_database_location`` is given
        * max_mz, only used when ``structure_database_location`` is given

        If ``from_subset`` is True then metids are read from stdin

        Uses fast option by default.
        """
        schema = colander.SchemaNode(colander.Mapping())
        self._addAnnotateSchema(schema)
        params = schema.deserialize(params)

        script = "{{magma}} annotate"
        script += " -p '{mz_precision}' -q '{mz_precision_abs}'"
        script += " -c '{ms_intensity_cutoff}' -d '{msms_intensity_cutoff}'"
        script += " -i '{ionisation_mode}' -b '{max_broken_bonds}'"
        script += " --precursor_mz_precision '{precursor_mz_precision}'"
        script += " --max_water_losses '{max_water_losses}'"
        script += " --call_back_url '{call_back_url}' "
        pmzp = params['precursor_mz_precision']
        ms_ic = params['ms_intensity_cutoff']
        msms_ic = params['msms_intensity_cutoff']
        script_substitutions = {
            'precursor_mz_precision': self.escape(pmzp),
            'mz_precision': self.escape(params['mz_precision']),
            'mz_precision_abs': self.escape(params['mz_precision_abs']),
            'ms_intensity_cutoff': self.escape(ms_ic),
            'msms_intensity_cutoff': self.escape(msms_ic),
            'ionisation_mode': self.escape(params['ionisation_mode']),
            'max_broken_bonds': self.escape(params['max_broken_bonds']),
            'max_water_losses': self.escape(params['max_water_losses']),
            'call_back_url': self.status_callback_url,
        }

        if params['structure_database'] is not colander.null:
            script += "--structure_database '{structure_database}'"
            script += " --db_options "
            script += "',{max_mim},{max_64atoms},{min_refscore}' "
            sd = self.escape(params['structure_database'])
            script_substitutions['structure_database'] = sd
            db_options = {'max_mim': self.escape(params['max_mz']),
                          'min_refscore': self.escape(params['min_refscore']),
                          'max_64atoms': self.restricted,
                          }
            script_substitutions.update(db_options)

        script = script.format(**script_substitutions)

        if from_subset:
            script += '-j - '

        if self.restricted:
            script += '--time_limit 3 '

        script += "--fast {db}\n"
        self.script += script

        return self

    def allinone(self, params):
        """Configure job query to do all sub commands in one go.

        ``params`` is a MultiDict

        See
        :meth:`~magmaweb.job.JobQuery.add_ms_data`,
        :meth:`~magmaweb.job.JobQuery.add_structures`,
        :meth:`~magmaweb.job.JobQuery.metabolize`,
        :meth:`~magmaweb.job.JobQuery.annotate` for params.

        """
        # only metabolize when params['metabolize'] exists
        metabolize = False
        if 'metabolize' in params:
            metabolize = True
            del(params['metabolize'])

        allin = self.add_ms_data(params)
        try:
            allin = allin.add_structures(params)
        except colander.Invalid as e:
            # no structures given
            if 'structure_database' in params and params['structure_database']:
                # structures will be added by database lookup during annotate
                pass
            else:
                sd = colander.SchemaNode(colander.String(),
                                         name='structure_database')
                msg = 'Either structures or structures_file'
                msg += ' or structure_database must be set'
                e.add(colander.Invalid(sd, msg))
                raise e
        if metabolize:
            allin = allin.metabolize(params)
        return allin.annotate(params, False)

    @classmethod
    def defaults(cls, selection=None):
        """Returns dictionary with default params

        If selection=='example' then params for example are returned.
        """
        if selection == 'example':
            return cls._example()
        elif selection == 'example2':
            return cls._example2()

        return dict(scenario=[
                              {'type': 'phase1', 'steps': '2'},
                              {'type': 'phase2', 'steps': '1'}
                              ],
                    ms_data_format='mzxml',
                    ms_data_area='',
                    ionisation_mode=1,
                    ms_intensity_cutoff=0.0,
                    msms_intensity_cutoff=5,
                    mz_precision=5.0,
                    mz_precision_abs=0.001,
                    abs_peak_cutoff=5000,
                    max_ms_level=10,
                    precursor_mz_precision=0.005,
                    max_broken_bonds=3,
                    max_water_losses=1,
                    )

    @classmethod
    def _example(cls):
        """Returns dictionary with params for example MS data set"""
        example_tree = [
            '353.087494: 69989984 (',
            '    191.055756: 54674544 (',
            '        85.029587: 2596121,',
            '        93.034615: 1720164,',
            '        109.029442: 917026,',
            '        111.045067: 1104891 (',
            '            81.034691: 28070,',
            '            83.014069: 7618,',
            '            83.050339: 25471,',
            '            93.034599: 36300,',
            '            96.021790: 8453',
            '            ),',
            '        127.039917: 2890439 (',
            '            57.034718: 16911,',
            '            81.034706: 41459,',
            '            83.050301: 35131,',
            '            85.029533: 236887,',
            '            99.045074: 73742,',
            '            109.029404: 78094',
            '            ),',
            '        171.029587: 905226,',
            '        173.045212: 2285841 (',
            '            71.013992: 27805,',
            '            93.034569: 393710,',
            '            111.008629: 26219,',
            '            111.045029: 339595,',
            '            137.024292: 27668,',
            '            155.034653: 145773',
            '            ),',
            '        191.055725: 17000514',
            '        ),',
            '    353.087097: 4146696',
            '    )'
        ]
        return dict(ms_data="\n".join(example_tree),
                    ms_data_format='mass_tree',
                    ionisation_mode=-1,
                    )

    @classmethod
    def _example2(cls):
        """Returns dictionary with params for example MS data set"""
        example_tree = [
            'C16H17O9: 69989984 (',
            '    C7H11O6: 54674544 (',
            '        C4H5O2: 2596121,',
            '        C6H5O: 1720164,',
            '        C6H5O2: 917026,',
            '        C6H7O2: 1104891 (',
            '            C5H5O: 28070,',
            '            C4H3O2: 7618,',
            '            C5H7O: 25471,',
            '            C6H5O: 36300,',
            '            C5H4O2: 8453',
            '            ),',
            '        C6H7O3: 2890439 (',
            '            C3H5O: 16911,',
            '            C5H5O: 41459,',
            '            C5H7O: 35131,',
            '            C4H5O2: 236887,',
            '            C5H7O2: 73742,',
            '            C6H5O2: 78094',
            '            ),',
            '        C7H7O5: 905226,',
            '        C7H9O5: 2285841 (',
            '            C3H3O2: 27805,',
            '            C6H5O: 393710,',
            '            C5H3O3: 26219,',
            '            C6H7O2: 339595,',
            '            C7H5O3: 27668,',
            '            C7H7O4: 145773',
            '            ),',
            '        C7H11O6: 17000514',
            '        ),',
            '    C16H17O9: 4146696',
            '    )',
        ]
        return dict(ms_data="\n".join(example_tree),
                    ms_data_format='form_tree',
                    ionisation_mode=-1,
                    )
