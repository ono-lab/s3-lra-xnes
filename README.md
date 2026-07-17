## S3-LRA-xNES: Step-size and Shape Separated Learning Rate Adaptation for Exponential Natural Evolution Strategy (GECCO2026 Best Paper Nominatation [[slide]](S3-LRA-xNES_GECCO2026.pdf))

### About

This repository contains the code for the paper
"[S3-LRA-xNES: Step-size and Shape Separated Learning Rate Adaptation for Exponential Natural Evolution Strategy](https://dl.acm.org/doi/10.1145/3795095.3805058)"
(Yes, OPEN ACCES !)
by Kosuke Ujihara, Masahiro Nomura, and Isao Ono, which has been accepted to [GECCO'26 (Best Paper Nominated at ENUM Track)](https://gecco-2026.sigevo.org/HomePage).

If you find this code useful in your research then please cite:

```bibtex
@inproceedings{ujihara2026s3,
author = {Ujihara, Kosuke and Nomura, Masahiro and Ono, Isao},
title = {S3-LRA-xNES: Step-size and Shape Separated Learning Rate Adaptation for Exponential Natural Evolution Strategy},
year = {2026},
isbn = {9798400724879},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
url = {https://doi.org/10.1145/3795095.3805058},
doi = {10.1145/3795095.3805058},
abstract = {The exponential natural evolution strategy (xNES) and the covariance matrix adaptation evolution strategy (CMA-ES) are effective for continuous black-box optimization, but their search performance strongly depends on hyperparameters such as population size and learning rates. Learning rate adaptation CMA-ES (LRA-CMA-ES) improves robustness by adapting the learning rates of the mean vector and covariance matrix; however, it converges slowly on ill-conditioned problems and often fails on multimodal problems. We argue that these issues partly arise because LRA-CMA-ES employs a single learning rate for the covariance matrix. This implicitly couples the updates of the distribution's step-size and shape components, ignoring the established principle in standard CMA-ES that these components require decoupled update dynamics. To address this limitation, we propose step-size and shape separated LRA-xNES (S3-LRA-xNES), which separately adapts learning rates for step-size and shape components. Experimental results demonstrate that S3-LRA-xNES achieves faster convergence on ill-conditioned problems and improved success rates on multimodal problems. Notably, S3-LRA-xNES outperforms LRA-CMA-ES by achieving a median 82.4\% reduction in function evaluations on ill-conditioned problems and consistently improving or maintaining 100\% success rates across 50 trials on multimodal problems, highlighting the effectiveness of fully separated learning rate adaptation in xNES.},
booktitle = {Proceedings of the Genetic and Evolutionary Computation Conference},
pages = {609–617},
numpages = {9},
keywords = {exponential natural evolution strategy, black-box optimization, learning rate adaptation},
location = {Centro Internacional de Convenciones CIC-ANDE, San Jose, Costa Rica},
series = {GECCO '26}
}
```

### Details (Dependencies, Example, and so on)

See the README files for our [Java](java/README.md) or [Python](python/README.md) implementations.
