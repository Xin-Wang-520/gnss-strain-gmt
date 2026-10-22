"""Known block medians: coordinate alignment, component outliers and sigma columns."""
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from GNSS_analysis.gnss_processing import ROOT,blockmedian_input,locate_gmt,gmt_environment
from run2_blockmedian_eu import arguments
import numpy as np
import json


def main():
    # Different records supply the east and north medians. Both components
    # must nevertheless use the same median coordinates and block membership.
    data=np.array([[80.01,30.01,1,5,.5,.3],
                   [80.02,30.02,99,100,.2,.4],
                   [80.03,30.03,3,2,.8,.9],
                   [80.11,30.11,10,20,1,2]])
    with tempfile.TemporaryDirectory(prefix='test_blockmedian_') as temp:
        result,meta=blockmedian_input(locate_gmt(None),data,[80,80.2,30,30.2],'0.1',Path(temp),gmt_environment())
        result=result[np.argsort(result[:,0])]
        np.testing.assert_allclose(result,[[80.02,30.02,3,5,.5,.4],[80.11,30.11,10,20,1,2]],atol=1e-12)
        assert meta['raw_station_count']==4 and meta['occupied_block_count']==2
        assert meta['merged_block_count']==1 and meta['largest_block_count']==3
    with patch('sys.argv',['run2_blockmedian_eu.py']):args=arguments()
    assert args.spacing=='30m' and args.region==[67,113,22,45]
    print('PASS: aligned component medians, singleton, representative sigma, block counts and new defaults.')


if __name__=='__main__':main()
