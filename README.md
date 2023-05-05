# Local Affine Preservation with Motion Consistency for Feature Matching
This is the Github repository for feature matching using local affine preservation matching with motion consistency.

# Requirement
MATLAB (test on 2021a)

[VLfeat](https://github.com/vlfeat/vlfeat) 

# Usage
- Demo dataset: ```.\127_r.JPG```,```.\127_l.JPG```,```.\127.mat```.
- Run ```demo_LAP_fast.m``` for mismatch removal using LAP (Local Affine Preservation) matching algorithm.
- ```evaluatePR.m``` calculates the precision_rate, Recall_rate and F1_score of the matching result.
- ```plot_matches.m``` visualize the matching result on the image pair.

  ![image](https://user-images.githubusercontent.com/87520254/236364113-5a13b260-5b54-446c-b307-206ab9cee338.png)

- ```plot_4c.m``` shows the distribution of ture correspondences and false correspondences.
  
![untitled](https://user-images.githubusercontent.com/87520254/236365469-55d11f07-a1ee-406f-9ed2-16955f6af233.png)



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
