# Author
Author:Wenkang Xu

Last Modified: Jul, 2023

Contact person email: 22131113@zju.edu.cn
# SLA-VBI
The source code for the paper for "Successive Linear Approximation VBI for Joint Sparse Signal Recovery and Dynamic Grid Parameters Estimation"

If you use this code or any (modified) part of it in any publication, please cite the paper: W. Xu, A. Liu, B. Zhou, and M.-J. Zhao, "Successive linear approximation VBI for joint sparse signal recovery and dynamic grid parameters estimation," 2023, arxiv:2307.09149. [Online]. Avaiable: https://arxiv.org/abs/2307.09149

# Brief Description:
This is the code for the scenario: FDD massive MIMO channel estimation.

1.main.m: The main function

2.SLA_VBI.m: The propsoed SLA_VBI algorith.

3.SLA_VBI_simplify.m: A faster version of the proposed SLA_VBI, where only a few dymanic grid parameters are updated.

4.SLA_IFVBI: The propsoed IFSLA_VBI algorith.

5.IFSLA_VBI_simplify.m: A faster version of the proposed IFSLA_VBI, where only a few dymanic grid parameters are updated.

6.OGSBI: Baseline 1.

7.SBL: Baseline 2.

8.Turbo_CS: Baseline 3.

9.Turbo_VBI: Baseline 4.

You should run main.m to simulate different algorithms.
