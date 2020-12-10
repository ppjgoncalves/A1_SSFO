## Figures 5 and S2

1. download DRC data from [figshare](https://figshare.com/s/70c3d487f11680acc6d6):

- ``dataArray_complete_v2.mat``
- ``newNegPosIdcs.mat``
- ``signalPower_07_v2.mat``
- ``signalPower_79-80_v2.mat``
- ``signalPower_119-120-153-154_v2.mat``
- ``testMask.mat``

2. install Python module [**lnpy**](https://github.com/arnefmeyer/lnpy).

3. run ``estimate_STRFs.ipynb`` to estimate the spectro-temporal receptive fields (STRFs) for all units. As the code has not been parallelised, this operation takes around 30-40s per unit, so a total of 5 to 6h on all units. Estimated STRFs are saved in folder ``/strfs/`` as pickled files ``.pkl``, and take around 10GB of storage.

4. run ``figs_STRFs.ipynb`` to generate Figures 5 and S2.

### Requirements

Python 2.7.
