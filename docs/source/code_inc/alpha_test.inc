.. code-block:: python

	#!/usr/bin/env python3.4
	# encoding: utf-8
	"""
	
	
	"""
	import sys
	import pyqms
	import pymzml
	from pymzml.plot_v2 import Factory
	import pickle
	import bisect
	
	
	def main():
	    EXMAPLES = [
	        ('ALELDSNNEK', 16531 , [566, 572]),
	        ('AVDDGVNTFK', 19453, [532, 538]),
	    ]
	    for peptide, scanID, mzRange in EXMAPLES:
	        match_and_plot( base_peptide= peptide, scanID = scanID, mzRange=mzRange)
	
	
	def match_and_plot( base_peptide = None, scanID = None, mzRange = None ):
	    MOLECULES = []
	    LABELS = [
	        "#Label:2H(4):1",
	        "#Label:13C(6):1",
	        # "#Label:2H(6):1",
	        ""
	    ]
	    for m in LABELS:
	        MOLECULES.append( base_peptide + m )
	
	    CHARGES = [ 2, 3, 4, 5 ]
	    lib = pyqms.IsotopologueLibrary(
	        molecules        = MOLECULES,
	        charges          = CHARGES,
	        params           = {
	            'REL_MZ_RANGE' : 20e-6,
	            'REL_I_RANGE' : 10,
	        }
	    )
	    r = pymzml.run.Reader('/Users/fu/Cloud/Sciebo/Collaborations/zmbh_2016_1/Q35D142-16-030-06.mzML')
	    # sweep( r, lib )
	    # exit(1)
	    spectrum = r[ scanID ]
	    condensed_results = lib.match_all(
	        mz_i_list = spectrum.centroidedPeaks,
	        file_name = "Q35D142-16-030-06",
	        spec_id   = spectrum['id'],
	        spec_rt   = spectrum['MS:1000016'],
	        results   = None
	    )
	    pf = Factory()
	    pf.newPlot( header='ZMBH // scanID {0}'.format( scanID ))
	    pf.add( spectrum.peaks, color=(225, 225, 225), name='raw data', mzRange=mzRange, plotType="Scatter" )
	    pf.add( spectrum.centroidedPeaks, color=(225, 225, 190), name='centroided raw data', mzRange=mzRange)
	    for m_key, i, match in condensed_results.extract_results():
	        molecule = lib.lookup['formula to molecule'][ m_key.formula ][0]
	        print( molecule, match.score)
	        # print( m_key, match )
	        matched_peaks = []
	        unmatched_peaks = []
	        for mmz, mi,ri, cmz, ci in match.peaks:
	            if mmz is None:
	                unmatched_peaks.append( ( cmz, match.scaling_factor * ci ) )
	            else:
	                matched_peaks.append( ( cmz, match.scaling_factor * ci ) )
	        if len(unmatched_peaks) != 0:
	            pf.add( unmatched_peaks , style='triangle.small', color=(200, 0, 0), name=molecule)
	        if len(matched_peaks) != 0:
	            pf.add( matched_peaks , style='triangle.small',  color=(0, 200, 0), name=molecule)
	    pf.save('akfufezan_match_example_{0}.html'.format( scanID ))
	    # pf.info()
	    return
	
	
	
	def sweep(run, lib ):
	    results = None
	    condensed_results = None
	    for offset in range(300):
	        spectrum = run[19376 + offset]
	        results = lib.match_all(
	            mz_i_list = spectrum.centroidedPeaks,
	            file_name = "Q35D142-16-030-06",
	            spec_id   = spectrum['id'],
	            spec_rt   = spectrum['MS:1000016'],
	            results   = results
	        )
	
	    pickle.dump(
	        results,
	        open('{0}_pyQms_results.pkl'.format('Q35D142-16-030-06'),'wb')
	    )
	    for m_key, i, match in results.extract_results():
	        molecule = lib.lookup['formula to molecule'][ m_key.formula ][0]
	        print( molecule, match.score, match.spec_id)
	    return
	
	
	def review_jpr():
	    '''
	    '''
	    TRIVIAL_NAMES = {
	        "GKGGKGLGKGGAKR#H3_Acetyl:1;H3_Acetyl:2;H3_Acetyl:3;H3_Acetyl:4" : 'HHHH',
	        "GKGGKGLGKGGAKR#H3_Acetyl:1;H3_Acetyl:2;H3_Acetyl:3;D3_Acetyl:4" : 'HHHD',
	        "GKGGKGLGKGGAKR#H3_Acetyl:1;H3_Acetyl:2;D3_Acetyl:3;D3_Acetyl:4" : 'HHDD',
	        "GKGGKGLGKGGAKR#H3_Acetyl:1;D3_Acetyl:2;D3_Acetyl:3;D3_Acetyl:4" : 'HDDD',
	        "GKGGKGLGKGGAKR#D3_Acetyl:1;D3_Acetyl:2;D3_Acetyl:3;D3_Acetyl:4" : 'DDDD'
	    }
	    CHARGES = [2,]
	    lib = pyqms.IsotopologueLibrary(
	        molecules        = list(TRIVIAL_NAMES.keys()),
	        charges          = CHARGES,
	        trivial_names = TRIVIAL_NAMES
	    )
	
	    pf = Factory()
	    pf.newPlot( header='Error')
	    for fn, formula in enumerate(sorted(lib.keys())):
	        for charge in CHARGES:
	            trivial_name = lib.lookup['formula to trivial name'][formula][0]
	            lib.print_overview(
	                formula,
	                charge=charge
	            )
	            for label_percentile in sorted(lib[ formula ]['env']):
	                number_of_isotopic_peaks = len(lib[ formula ]['env'][ label_percentile ]['abun'])
	                peaks = []
	                anno = []
	                for n in range(number_of_isotopic_peaks):
	                    # mass = lib[ formula ]['env'][ label_percentile ]['mass'][n]
	                    # abun = lib[ formula ]['env'][ label_percentile ]['abun'][n]
	                    mzc1 = lib[ formula ]['env'][ label_percentile ][charge]['mz'][n]
	                    relabun = lib[ formula ]['env'][ label_percentile ]['relabun'][n]
	                    # c_pos = lib[ formula ]['env'][ label_percentile ]['c_peak_pos'][ n ]
	                    peaks.append( (mzc1, relabun ))
	                    if mzc1 is not None:
	                        anno.append( (mzc1, mzc1, relabun, trivial_name ))
	
	                pf.add( peaks,color=(0, 50*fn, fn*20), name=trivial_name)
	                pf.add( anno, style='triangle.small', color=(0, 50*fn, fn*20), name='{0} annotation'.format( trivial_name ))
	    pf.save('')
	    pf.info()
	
	
	if __name__ == '__main__':
	    main()
