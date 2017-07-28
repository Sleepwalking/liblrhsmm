liblrhsmm
===

*liblrhsmm* is a C implementation of inference and training algorithms for Hidden Semi-Markov Models (HSMM) with left-to-right topology. This variant of HSMM is particularly useful for time series alignment, for example, text-and-speech/audio/video alignment. Though such topology could **not** be directly used in speech recognition applications, *liblrhsmm* could possibly train speech recognizers and synthesizers.

The specific type of model supported by this library is Hidden Semi-Markov Model with left-to-right state transitions (and in fact, arbitrary state transitions are supported, though in a constrained manner), single normal distribution as state duration, and multi-stream diagonal-covariance Gaussian Mixture Model as output (emission) distribution. This restrictive topology allows **state tying configuration to be stored separately, outside of the model description file**; stream-level tying and tying of duration distributions are also supported. Moreover, *liblrhsmm* supports training/doing inference on the HSMM model as a left-to-right first-order HMM: the next-state transition probability is treated as the inverse of mean duration.

| Property | Description (as HSMM) | Description (as first-order HMM) |
| --- | --- | --- |
| Model Type | Hidden Semi-Markov Model | Hidden Markov Model |
| Topology | Left-to-right **with optional jumps** | Left-to-right, Reflexive |
| Duration Distribution | Single-mixture, 1-dim Normal Distribution | Geometric Distribution |
| Emission Distribution | Multi-stream, Multi-mixture, Multi-dim Normal Distribution | Multi-stream, Multi-mixture, Multi-dim Normal Distribution |
| Covariance Matrix | Diagonal | Diagonal |

Inference algorithms of *liblrhsmm* are slightly different from the ones found in many literatures. Forward and backward probabilities at time t, state j are defined as the log probability of the joint or posterior of having state j **ending at time t** and observing the rest of the output sequence and the results are 2D matrices (marginalizing the state durations); state occupancy probability is calculated from marginalization of durational occupancy probability. The training algorithm is modified from the standard Baum-Welch algorithm for Hidden Markov Models. Viterbi training is not included in this library, but its implementation is trivial.

`serial.h` and `serial.c` provide very basic serialization of model, data and segmentation. The serialized binary is stored in [*messagepack*](http://msgpack.org/index.html) format. You may enable/disable this feature when compiling *liblrhsmm* at your will.

License
---

GPLv3

How to use
---

Users are expected to be familiar with Hidden Markov Models. The reference section lists several publications which you may find helpful understanding the models.

`test/test-random-model.c` is a rather comprehensive example making use of most *liblrhsmm* features. It starts with creating a random model and some random data, then train the model on the generated data, increase the number of mixtures, and run the training for another few iterations.

`test/test-state-jumps.c` demonstrates adding skip, self-loop and backward transitions to a LR-HSMM.

Reference
---

* Hua, Kanru. ["Doing Inference in a LR-HSMM"](https://github.com/Sleepwalking/prometheus-spark/blob/master/writings/inference-lr-hsmm-hua-2015.pdf), 2015. Web.
* Fink, Gernot A. *Markov models for pattern recognition: from theory to applications*. Springer Science & Business Media, 2014.
* Young, Steve, et al. *The HTK Book (for HTK Version 3.4)*. Cambridge University Engineering Department, 2009.
* Yu, Shun-Zheng. "Hidden semi-Markov models". *Artificial Intelligence* 174 (2010): 215-243. Print.
* Zen, Heiga, et al. "Hidden Semi-Markov Model Based Speech Synthesis". *ICSLP* (2004). Conference.
