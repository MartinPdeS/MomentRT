# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/MomentRT/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                            |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|-------------------------------- | -------: | -------: | -------: | -------: | ---------: | --------: |
| MomentRT/boundary.py            |       71 |        5 |       16 |        4 |     89.66% |105, 111, 124, 170, 173 |
| MomentRT/convergence.py         |      107 |        8 |       58 |       11 |     88.48% |80, 85, 94, 97, 105, 183-\>202, 204-\>207, 217-\>225, 230, 281, 286 |
| MomentRT/detector.py            |      116 |        9 |       32 |        9 |     87.84% |86, 108, 117, 134, 137, 195, 298, 405, 452 |
| MomentRT/fem.py                 |       92 |       11 |       16 |        1 |     87.04% |218-232, 284-294 |
| MomentRT/fields.py              |       84 |       10 |       26 |       10 |     81.82% |55, 58, 85, 188, 263, 268, 355, 362, 400, 434 |
| MomentRT/geometry.py            |       64 |        3 |        6 |        3 |     91.43% |30, 33, 280 |
| MomentRT/material.py            |       92 |        7 |       38 |        6 |     88.46% |134, 163, 174, 303, 348, 367-368 |
| MomentRT/materials/liquids.py   |       64 |       51 |       14 |        0 |     16.67% |50-76, 88-120, 124, 181-268, 272 |
| MomentRT/materials/skin.py      |       46 |       36 |       12 |        0 |     17.24% |57-109, 113, 169-217, 221 |
| MomentRT/medium.py              |      125 |        6 |       50 |        8 |     92.00% |165, 184, 268, 332, 335, 501, 514-\>513, 525-\>528 |
| MomentRT/mesh.py                |      446 |       89 |      132 |       17 |     79.58% |60, 105, 116, 170-269, 286-320, 415, 473, 491, 525, 547, 558, 643, 661, 681, 800, 823, 832, 860-861, 1042, 1060-\>1063, 1107, 1153, 1172-1184, 1204 |
| MomentRT/phase\_function.py     |       82 |        9 |       26 |        8 |     84.26% |13, 22, 53, 69, 97, 121, 177, 267, 321, 364-\>369 |
| MomentRT/plotting.py            |      209 |       14 |       82 |       15 |     90.03% |25, 156, 204, 403, 406, 411, 418, 421, 428, 532, 564-\>567, 632, 839, 884, 939-\>946, 983 |
| MomentRT/results.py             |      215 |       34 |       58 |        7 |     79.12% |219-\>222, 224-\>227, 246, 260-\>278, 281-\>286, 688-\>exit, 749-752, 818-939, 1005, 1075-1104 |
| MomentRT/simulation.py          |      246 |       60 |      130 |       30 |     71.81% |100, 105-113, 125, 193, 198, 201, 204, 207, 210, 215, 218, 220-\>235, 232-233, 238-\>253, 247, 250-251, 260, 283, 313-318, 320-\>348, 336-337, 348-\>523, 366, 390-417, 433, 443-446, 451-\>485, 466, 482-483, 508-\>485, 528-535, 649-654, 662-665, 781, 812 |
| MomentRT/source.py              |       99 |        9 |       18 |        7 |     84.62% |76, 79, 174, 185, 188, 524, 538, 568-569 |
| MomentRT/transport/diffusion.py |      156 |       16 |       50 |       10 |     85.44% |69, 74, 130, 214, 321, 324, 327, 373-392, 562, 569 |
| MomentRT/transport/p3.py        |      153 |        2 |       38 |        2 |     97.91% |  385, 390 |
| MomentRT/transport/sp3.py       |       81 |        2 |       14 |        2 |     95.79% |  216, 223 |
| MomentRT/validation.py          |       22 |        0 |        2 |        1 |     95.83% |   92-\>95 |
| **TOTAL**                       | **2703** |  **381** |  **836** |  **151** | **82.71%** |           |

5 files skipped due to complete coverage.


## Setup coverage badge

Below are examples of the badges you can use in your main branch `README` file.

### Direct image

[![Coverage badge](https://raw.githubusercontent.com/MartinPdeS/MomentRT/python-coverage-comment-action-data/badge.svg)](https://htmlpreview.github.io/?https://github.com/MartinPdeS/MomentRT/blob/python-coverage-comment-action-data/htmlcov/index.html)

This is the one to use if your repository is private or if you don't want to customize anything.

### [Shields.io](https://shields.io) Json Endpoint

[![Coverage badge](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/MartinPdeS/MomentRT/python-coverage-comment-action-data/endpoint.json)](https://htmlpreview.github.io/?https://github.com/MartinPdeS/MomentRT/blob/python-coverage-comment-action-data/htmlcov/index.html)

Using this one will allow you to [customize](https://shields.io/endpoint) the look of your badge.
It won't work with private repositories. It won't be refreshed more than once per five minutes.

### [Shields.io](https://shields.io) Dynamic Badge

[![Coverage badge](https://img.shields.io/badge/dynamic/json?color=brightgreen&label=coverage&query=%24.message&url=https%3A%2F%2Fraw.githubusercontent.com%2FMartinPdeS%2FMomentRT%2Fpython-coverage-comment-action-data%2Fendpoint.json)](https://htmlpreview.github.io/?https://github.com/MartinPdeS/MomentRT/blob/python-coverage-comment-action-data/htmlcov/index.html)

This one will always be the same color. It won't work for private repos. I'm not even sure why we included it.

## What is that?

This branch is part of the
[python-coverage-comment-action](https://github.com/marketplace/actions/python-coverage-comment)
GitHub Action. All the files in this branch are automatically generated and may be
overwritten at any moment.