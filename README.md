# ESCG_Python
 Python code for Rock-paper-scissors-lizard-Spock (RPSLS) style of Evolutionary Spatial Cyclic Game (ESCG).

Everything is in the one Python file: ESCG_RPSLS.py. The comments at the top of the Python file explain what is going on. 
If you want to draw network diagrams, you will need to install Graphviz. If not, comment-out (or delete) all the code that calls graphviz. This is only needed if you want diagrams of the dominance networks; the ESCG simulation doesn't depend on this at all.

The code in this repo started as an implementation and extension of the kind of ablated RPSLS models explored in the paper "Species coexistence in spatial cyclic game of five species" by Linwu Zhong, Liming Zhang, Haihong Li, Qionglin Dai & Junzhong Yang, _Chaos, Solitons and Fractals_, 156 (2022) 111806. Referred to in comments below as ZZLDY22.

The code incorporates the "Orginal Elementary Step" (OES) as used by ZZLDY22 and also the "Revised Elementary Step" (RES) introduced in this paper:
Cliff, D. (2024a) "Never Mind The No-Ops: Faster and Less Volatile Simulation Modelling of Co-Evolutionary Species Interactions via Spatial Cyclic Games". Accepted for presentation/publication in _Proceedings of the 36th European Modeling and Simulation Symposium (EMSS2024)_; Tenerife, Spain; 18th-20th September 2024.
Available at SSRN: https://ssrn.com/abstract=4883174 .

The code also implements parameterised circulant networks as reported in this paper:
Cliff, D. (2024b) "Tournament versus Circulant: On Simulating 7-Species Evolutionary Spatial Cyclic Games with Ablated Predator-Prey Networks as Models of Biodiversity". Accepted for presentation/publication in _Proceedings of the 36th European Modeling and Simulation Symposium (EMSS2024)_. Tenerife, Spain, 18-20 September 2024.
Available at SSRN: https://ssrn.com/abstract=4961889 .

This code can be used to do the kind of studies of long-term coevolutionary dynamis explored in this paper:
Cliff, D. (2024c) "On Long-Term Species Coexistence in Five-Species Evolutionary Spatial Cyclic Games with Ablated and Non-Ablated Dominance Networks". This is currently an unpublished manuscript, submitted to a journal for peer-review, but it is available as a pre-print at SSRN: https://ssrn.com/abstract=4967488. Please beware that the Python code available here would take a _very_ long time to run those kind of experiments -- those experiments were generated instead from a re-implementation of the Python code available here, rewritten in the C programming language and compiled with full optimization. 
