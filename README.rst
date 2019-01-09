Introduction
############

*pyQms enables universal and accurate quantification of mass spectrometry data*


|doc-status| |build-status-travis| |build-status-appveyor| |pypi|

.. |doc-status| image:: https://readthedocs.org/projects/pyqms/badge/?version=latest
   :target: http://pyqms.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. |build-status-travis| image:: https://travis-ci.org/pyQms/pyqms.svg?branch=master
   :target: https://travis-ci.org/pyQms/pyqms
   :alt: Travis CI status

.. |build-status-appveyor| image:: https://ci.appveyor.com/api/projects/status/n4x2ug7h3ce4d49y?svg=true
   :target: https://ci.appveyor.com/project/fufezan-lab/pyqms
   :alt: AppVeyor CI status

.. |pypi| image:: https://img.shields.io/pypi/v/pyqms.svg
   :target: https://pypi.org/project/pyqms/


Summary
*******

pyQms is an extension to Python that offers amongst other things

    a) fast and accurate quantification of all high-res LC-MS data

    b) full labeling and modification flexibility

    c) full platform independence



Abstract
********

Quantitative mass spectrometry (MS) is a key technique in many research areas (Yates III et al. 2009), including proteomics, metabolomics, glycomics, and lipidomics. Because all of the corresponding molecules can be described by chemical formulas, universal quantification tools are highly desirable. Here we present pyQms, an open-source software for accurate quantification of all types of molecules measurable by MS. pyQms uses isotope pattern matching which offers accurate quality assessment of the quantification and the ability to directly incorporate mass spectrometer accuracy. pyQms is, due to its universal design, applicable to every research field, labeling strategy, and acquisition technique. This opens ultimate flexibility for researchers to design experiments employing innovative and hitherto unexplored labeling strategies. Importantly, pyQms performs very well to accurately quantify partially labeled proteomes in large-scale and high-throughput, the most challenging task for a quantification algorithm.

            -- Leufken, J., Niehues, A., Hippler, M., Sarin, L. P., Hippler, M., Leidel, S. A., and Fufezan, C. (2017) pyQms enables universal and accurate quantification of mass spectrometry data. MCP In Press

Link to manuscript.

http://www.mcponline.org/content/early/2017/07/20/mcp.M117.068007.abstract 


pyQms module
************
At its core, pyQms is a Python module that allows a isotope pattern library to
be initialized and any list of (mz, intensity) to be matched against the library,
yielding a mScore.

Documentation
*************

http://pyqms.readthedocs.io/en/latest/


Implementation
**************

pyQms requires Python3.4+ .


The module is freely available on pyqms.github.io or pypi,
published under MIT LGPL and requires no additional modules to be installed.
For fast spectra from mzML access we recommend pymzML (Bald et al. 2012).
For example scripts it is necessary to install pymzML as well or
change the code for alternated spectra access. For some scripts also the
openpyxl module is required.

.. _download_instructions:

Download
********

Get the latest version via github
    | https://github.com/pyQms/pyQms


Citation
********

Please cite us when using pyQms in your work.


The original publication can be found here:
    Leufken, J., Niehues, A., Hippler, M., Sarin, L. P., Hippler, M., Leidel, S. A., and Fufezan, C. (2017) pyQms enables universal and accurate quantification of mass spectrometry data. Mol. Cell. Proteomics 16, 1736–1745
    

Full article
============

http://www.mcponline.org/content/16/10/1736


Early access article version
============================


http://www.mcponline.org/content/early/2017/07/20/mcp.M117.068007.abstract 

DOI
===

10.1074/mcp.M117.068007


.. _installation_instructions:

Installation
************

Install requirements::

    user@localhost:~$ cd pyqms
    user@localhost:~/pyqms$ pip3.4 install -r requirements.txt


.. note:

    Pip is included in Python 3.4 and higher. However, it might not be
    included in in your system's PATH environment variable.
    If this is the case, you can either add the Python scripts directory to your
    PATH env variable or use the path to the pip.exe directly for the
    installation, e.g.: ~/Python34/Scripts/pip.exe install -r requirements.txt


Install pyQms::

    user@localhost:~/pyqms$ python3.4 setup.py install

.. note:

    Consider to use a Python virtual environment for easy installation and use. 
    Further, usage of python3.4+ is recommended.


pyQms can be also be installed via pip::
    
    pip install pyqms

.. note:
    
    For obtaining the latest version of pyQms please use the github repo.



(You might need administrator privileges to write in the Python site-package folder.
On Linux or OS X, use ```sudo python setup.py install``` or write into a user folder
by using this command ```python setup.py install --user```. On Windows, you have to
start the command line with administrator privileges.)

pyQms docs recompiling and extending
====================================

You will require sphinx and other packages to build the documentation from
scratch. We recommend to use a Python virtual environment for easy installation
and use.


Tests
*****

Run nosetests in root folder. You might need to install `nose`_ for Python3
first. Then just execute::

    user@localhost:~/pyqms$ nosetests3

to test the package.

.. _nose:
    https://nose.readthedocs.org/en/latest/




LICENSE
*******

This software is under MIT license, please refer to LICENSE for full license.



Publications and project using pyQms for quantification
*******************************************************
        
 | - Hohner, R., Barth, J., Magneschi, L., Jaeger, D., Niehues, A., Bald, T., Grossman, A., Fufezan, C., and Hippler, M. (2013) The Metabolic Status Drives Acclimation of Iron Deficiency Responses in Chlamydomonas reinhardtii as Revealed by Proteomics Based Hierarchical Clustering and Reverse Genetics. **Mol. Cell. Proteomics** 12, 2774–2790 `Pubmed <https://www.ncbi.nlm.nih.gov/pubmed/23820728>`_
 | - Barth, J., Bergner, S. V., Jaeger, D., Niehues, A., Schulze, S., Scholz, M., and Fufezan, C. (2014) The Interplay of Light and Oxygen in the Reactive Oxygen Stress Response of Chlamydomonas reinhardtii Dissected by Quantitative Mass Spectrometry. **Mol. Cell. Proteomics** 13, 969–989 `Pubmed <https://www.ncbi.nlm.nih.gov/pubmed/24482124>`_ 
 | - Kukuczka, B., Magneschi, L., Petroutsos, D., Steinbeck, J., Bald, T., Powikrowska, M., Fufezan, C., Finazzi, G., and Hippler, M. (2014) Proton Gradient Regulation5-Like1-Mediated Cyclic Electron Flow Is Crucial for Acclimation to Anoxia and Complementary to Nonphotochemical Quenching in Stress Adaptation. **Plant Physiol.** 165, 1604–1617 `Pubmed <https://www.ncbi.nlm.nih.gov/pubmed/24948831>`_ 
 | - Alings, F., Sarin, L. P., Fufezan, C., Drexler, H. C. A., and Leidel, S. A. (2015) An evolutionary approach uncovers a diverse response of tRNA 2-thiolation to elevated temperatures in yeast. **RNA** 21, 202–212 `Pubmed <https://www.ncbi.nlm.nih.gov/pubmed/25505025>`_ 
 | - Bergner, S. V., Scholz, M., Trompelt, K., Barth, J., Gäbelein, P., Steinbeck, J., Xue, H., Clowez, S., Fucile, G., Goldschmidt-Clermont, M., Fufezan, C., and Hippler, M. (2015) STATE TRANSITION7-Dependent Phosphorylation Is Modulated by Changing Environmental Conditions, and Its Absence Triggers Remodeling of Photosynthetic Protein Complexes. **Plant Physiol.** 168, 615–634 `Pubmed <https://www.ncbi.nlm.nih.gov/pubmed/25858915>`_ 
 | - Hochmal, A. K., Zinzius, K., Charoenwattanasatien, R., Gäbelein, P., Mutoh, R., Tanaka, H., Schulze, S., Liu, G., Scholz, M., Nordhues, A., Offenborn, J. N., Petroutsos, D., Finazzi, G., Fufezan, C., Huang, K., Kurisu, G., and Hippler, M. (2016) Calredoxin represents a novel type of calcium-dependent sensor-responder connected to redox regulation in the chloroplast. **Nat. Commun.** 7, 11847 `Pubmed <https://www.ncbi.nlm.nih.gov/pubmed/27297041>`_ 
 | - Pfannmüller, A., Leufken, J., Studt, L., Michielse, C. B., Sieber, C. M. K., Güldener, U., Hawat, S., Hippler, M., Fufezan, C., and Tudzynski, B. (2017) Comparative transcriptome and proteome analysis reveals a global impact of the nitrogen regulators AreA and AreB on secondary metabolism in Fusarium fujikuroi. PLoS One in press, 1–27 `Pubmed <https://www.ncbi.nlm.nih.gov/pubmed/28441411>`_ 

Contact information
*******************

Please refer to:

    | Dr. Christian Fufezan
    | Cellzome
    | Molecular Discovery Research
    | GlaxoSmithKline
    | 69117 Heidelberg
    | Germany
    | eMail: christian@fufezan.net
