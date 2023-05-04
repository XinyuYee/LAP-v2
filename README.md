# Local Affine Preservation with Motion Consistency for Feature Matching
This is the Github repository for feature matching using local affine preservation matching with motion consistency.

# Requirement
MATLAB (test on 2021b)

[VLfeat](https://github.com/vlfeat/vlfeat) 

# Usage
- Demo dataset: ```.\127_r.JPG```,```.\127_l.JPG```,```.\127.mat```.
- Run ```demo_LAP_fast.m``` for mismatch removal using LAP (Local Affine Preservation) matching algorithm.
- Run ```evaluatePR.m```
- Run ```plot_matches.m```

# Citation

If you used the code, please cite:
```
@ARTICLE{ye2022local,
  author={Ye, Xinyu and Ma, Jiayi and Xiong, Huilin},
  journal={IEEE Transactions on Geoscience and Remote Sensing}, 
  title={Local Affine Preservation With Motion Consistency for Feature Matching of Remote Sensing Images}, 
  year={2022},
  volume={60},
  number={},
  pages={1-12},
  doi={10.1109/TGRS.2021.3128292}}
```
