.. quick_start:

Quick start
***********

Download and installation
=========================

Please :ref:`download_instructions` and install pyQms following these 
:ref:`installation_instructions` instructions. Please consider using a 
virtual environment (e.g. using the excellent `virtualenvwrapper`_) for using 
and developing pyQms.

.. _virtualenvwrapper:
    https://virtualenvwrapper.readthedocs.io/en/latest/



Matching a peak list
====================

Let's start with a most simple example: Mathing a single peptide on a predefined
peak list. Start a Python (3.4+) console and start quantifying in 4 steps:

First import pyQms: ::

    import pyqms

Second, initialize a isotopologue library (:py:class:`pyqms.IsotopologueLibrary`) 
using 'DDSPDLPK' as the example peptide (from BSA example file) and the charge 
state 2: ::
    
    lib = pyqms.IsotopologueLibrary(
        molecules  = [ 'DDSPDLPK' ],
        charges    = [ 2 ],
    )

Third, match the library on the provided peak list. You can find a peak list 
:ref:`here<example_peak_list>`, which will produce a match with this peptide.
Copy and paste the peak list into the Python console.

Fourth, use the :py:func:`pyqms.IsotopologueLibrary.match_all` function to quantify 
the peptide using the peak list: ::

    results = lib.match_all(
        mz_i_list = peak_list,
        file_name = 'test',
        spec_id   = 1165,
        spec_rt   = 29.10,
        results   = None
    )

Done! The peptide has been quantified in the given peak list. Please continue
with the next section to learn how to access and process the results.

.. note::

    The keyword arguments `file_name`, `spec_id` and `spec_rt` are hardcoded in
    this example case. In the advanced examples these information (as well as 
    the peak list) are parsed from the mzML file directly.


.. _access results:

Access and interpret the results
================================

The results object represents the :py:class:`pyqms.Results` class and is
organized as a dictionary: ::
    
    results.keys()

Will give the following output: ::

    dict_keys(
        [
            m_key(
                file_name='test',
                formula='C(37)H(59)N(9)O(16)',
                charge=2,
                label_percentiles=(('N', '0.000'),)
            )
        ]
    )

The keys of the :py:class:`pyqms.Results` class are :py:func:`namedtuple` with
the following field_names:

    * file_name
    * formula
    * charge
    * label_percentiles

`file_name` related to the original file name of the LC-MS/MS runs, `formula` is
the molecular formula of the input molecule/peptide, `charge` refers to the
charge state of the matched isotope envelope and `label_percentile` indicates
the labeling of the molecule. Default behaviour is to use the natural abundance
of the element isotopes (default this fieldname is set to 0% artificical
enrichment of nitrogen i.e. ('N','0.000') in a :py:obj:`tuple` of multiple
possible labeling percentiles i.e. (('N','0.000'),).

.. note::

    Every input molecule (e.g. peptide 'DDSPDLPK' ) will be converted to its
    molecular formula ('C(37)H(59)N(9)O(16)') in `Hill notation`_ by pyQms. To
    map between the peptide and formula, please use the integrated lookups, i.e. 
    `results.lookup['formula to molecule']` or results.lookup['molecule to formula']`.
    Please consider, that multiple molecules can have the same formula, therefor
    e.g. `results.lookup['formula to molecule']['C(37)H(59)N(9)O(16)']` is by 
    default a list.


.. _Hill Notation:
    https://en.wikipedia.org/wiki/Hill_system


For each of the keys one will get the following dict: ::

    {
        'data': [
            match(
                spec_id=1165,
                rt=29.1,
                score=0.9606609710868856,
                scaling_factor=40.75802642055527,
                peaks=(
                    (443.7112735313511, 2517650.0, 1.0, 443.7112648946701, 62091), (444.21248374593875, 1156173.75, 0.4459422196277157, 444.2127374486285, 27689), 
                    (444.71384916266277, 336326.96875, 0.12958327918547244, 444.7142840859656, 8046), 
                    (445.21533524843596, 58547.0703125, 0.02805309805863953, 445.21582563050043, 1742)
                )
            )
        ],
        'max_score': 0.9606609710868856, 
        'len_data': 1, 
        'max_score_index': 0
    }

The keys on the top level of this dictionary are:
    
    * data
    * max_score
    * len_data
    * max_score_index

While `len_data` will indicate how many spectra were matched for the formula in the
repective key, `max_score` and `max_score_index` provides the maximum score,
which was obtained during matching and the index of this match in the `data` 
list, respectively. The `data` list contains matches for all single spectra
as :py:func:`namedtuple`. The following fieldnames are contained in each `match`:

    * spec_id
    * rt
    * score
    * scaling_factor
    * peaks

Besides the given input information on the spectrum like the spectrum ID
(`spec_id`) and the retention time (`spec_rt`) the mScore of the match is
provided (`score`) as well as the determined amount/intensity of the molecule in
the spectrum (`scaling_factor`). Furthermore, detailed match information are
given in `peaks`. This tuple contains for each peak of the isotopologue the
following information in this order:

    * The measured (and matched) m/z value of the isotope peak in the spectrum
    * The measured intensity of the isotope peak in the spectrum
    * The relative intensity of the isotopologue peak to the monoisotopic peak
    * The calculated m/z value of the isotope peak of the input molecule
    * The calculated intensity of the isotope peak of the input molecule

These information can be processed to further analyze, besides the mScore, the 
quality of the `match`.

.. note::

    Please note, that measured m/z entry in `peaks` can be `None`, if this peak
    was not found in the input data.

We have now seen, how peptides/molecules can be quantified and how the results
can be accessed.

.. note::
    
    The :py:class:`pyqms.Results` class offers several functions to access,
    process and visualize the data. E.g. :py:func:`pyqms.Results.extract_results` 
    provides and iterator yielding `key`, `i`, `entry`.
    The `key` is the :py:func:`namedtuple` containing the molecules information,
    `i` is the position of `entry` in results[key]['data'] and `entry` is the 
    `match` :py:func:`namedtuple`.


Quantify peptides in a whole LC-MS run
======================================

This part will describe how to process a whole LC-MS/MS run and quantify 
multiple peptides in one batch. This example assumes you have started your
Python console in the pyqms base folder.

For this example we will use pymzML, which is used to parse mzML files and
retrieve the spectra and meta data used for quantification. pymzML will be
installed as a requirement (See: :ref:`installation_instructions`).

We start again by importing pyQms and initializing a isotopologue library 
( :py:class:`pyqms.IsotopologueLibrary` ): ::
    
    import pyqms
    lib = pyqms.IsotopologueLibrary(
        molecules        = [
            'HLVDEPQNLIK',
            'YICDNQDTISSK',
            'DLGEEHFK'
        ],
        charges          = [2, 3, 4, 5],
    )

We need to import pymzML and initialize the run. Note, that the path to the BSA1
mzML file ('data/BSA1.mzML') may have to be adjusted. This file can be 
downloaded using this example script `get_example_BSA_file` 
(See: :ref:`get bsa file`) and can then be found under the 'data' folder in the
pyqms base folder. ::
    
    import pymzml
    run = pymzml.run.Reader( 'data/BSA1.mzML' )

We now iterate over the spectra in the mzML file and quantify all peptides in
all MS1 spectra. Before we start the loop we set the results variable to `None`. 
Please note, that the `results` variable is iteratively passed to :py:func:`pyqms.IsotopologueLibrary.match_all`. This will lead to one `results`
object, which combines quantifications for all peptides in every spectra. See
also description above (see: `access results`_) or refer directly to the :py:class:`pyqms.Results`: class 
: ::
    
    results = None
    for spectrum in run:
        scan_time = spectrum['scan time']
        spec_id   = spectrum['id']
        if spectrum['ms level'] == 1:
            results = lib.match_all(
                mz_i_list = spectrum.centroidedPeaks,
                file_name = 'BSA1',
                spec_id   = spec_id,
                spec_rt   = scan_time,
                results   = results
            )
.. note ::
    
    pymzML centroids spectra if these are not already centroided, if
    `spectrum.centroidedPeaks` is accessed.

The `results` can now be accessed as described above (see: `access results`_).
Furthermore the :py:class:`pyqms.Results` class can be pickled: ::

    import pickle
    pickle.dump(
        results,
        open(
            'data/BSA1_pyQms_results.pkl',
            'wb'
        )
    )





For further examples and how to use the adaptor functions, please refer to the
next section.


Use the adaptors, Luke 
======================

The :ref:`adaptors` functions are useful for parsing a set of identified peptides 
(e.g. from Ursgal_ result files; `Ursgal documentation`_) including retention time 
information for determining the maximum intensity of every (identified) peptide 
in the LC-MS/MS measurement. Furthermore, adaptors can be added to e.g. read 
results of other analysis pipelines and tools.

The current adaptor to read Ursgal_ results can be used as follows for the
shipped identification result file of the database search engine OMSSA. Please
note, that if the adaptors are used one need to define fixed modifications like
Carbamidomathylation as presented. This modification and the molecules will then 
be correctly formatted as input for pyqms: ::
    
    import pyqms
    import pyqms.adaptors
    input_fixed_labels = {
        'C' : [
            {
                'element_composition' : {
                    'O'   : 1,
                    'H'   : 3,
                    '14N' : 1,
                    'C'   : 2
                },
                'evidence_mod_name': 'Carbamidomethyl'
            },
        ]
    }
    formatted_fixed_labels, evidence_lookup, molecules = pyqms.adaptors.parse_evidence(
        fixed_labels   = input_fixed_labels,
        evidence_files = [ 'data/BSA1_omssa_2_1_9_unified.csv' ],
    )

The returned objects can be used a direct input for the pyQms :py:class:`pyqms.IsotopologueLibrary`. The advantage of parsing evidence files is, that 
MS2 identification information is added to the results and can e.g. be used for 
defining RT windows for a correct quantification of every peptide: ::

    lib = pyqms.IsotopologueLibrary(
        molecules    = molecules,
        charges      = [1, 2, 3, 4, 5],
        fixed_labels = formatted_fixed_labels,
        evidences    = evidence_lookup
    )


.. _Ursgal:
    https://github.com/ursgal/ursgal

.. _Ursgal Documentation:
    http://ursgal.readthedocs.io/en/latest/


Further examples and more adavanced usage
=========================================

Please refer to the :ref:`example scripts` section for more usage examples and 
ready-to-go Python scripts for quantification, data analysis and visualization.




