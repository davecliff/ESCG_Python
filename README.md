# ESCG_Python
 Python code for Rock-paper-scissors-lizard-Spock (RPSLS) style of Evolutionary Spatial Cyclic Game (ESCG)

 Everything is in the one Python file. The comments at the top of the Python file explain what is going on. 
 You will need to install Graphviz, or comment-out (or delete) all the code that calls graphviz -- only needed if you want diagrams of the dominance networks.

 This started as an implementation and extension of the kind of ablated RPSLS models explored in the paper "Species coexistence in spatial cyclic game of five species" by Linwu Zhong, Liming Zhang, Haihong Li, Qionglin Dai & Junzhong Yang, _Chaos, Solitons and Fractals_, 156 (2022) 111806. Referred to in comments below as ZZLDY22.

The code incorporates the "Orginal Elementary Step" (OES) as used by ZZLDY22 and also the "Revised Elementary Step" (RES) introduced in this paper:
Cliff, D. (2024a) "Never Mind The No-Ops: Faster and Less Volatile Simulation Modelling of Co-Evolutionary Species Interactions via Spatial Cyclic Games". Accepted for presentation/publication in _Proceedings of the 36th European Modeling and Simulation Symposium (EMSS2024)_; Tenerife, Spain; 18th-20th September 2024.
Available at SSRN: https://ssrn.com/abstract=4883174

The code also implements parameterised circulant networks as reported in this paper:
Cliff, D. (2024b) "Tournament versus Circulant: On Simulating 7-Species Evolutionary Spatial Cyclic Games with Ablated Predator-Prey Networks as Models of Biodiversity". Accepted for presentation/publication in _Proceedings of the 36th European Modeling and Simulation Symposium (EMSS2024)_. Tenerife, Spain, 18-20 September 2024.
Available at SSRN: https://ssrn.com/abstract=4961889
