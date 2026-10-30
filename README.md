# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/MomentRT/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                            |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|-------------------------------- | -------: | -------: | -------: | -------: | ---------: | --------: |
| MomentRT/boundary.py            |       71 |        5 |       16 |        4 |     89.66% |105, 111, 124, 170, 173 |
| MomentRT/detector.py            |      116 |        9 |       32 |        9 |     87.84% |86, 108, 117, 134, 137, 195, 298, 405, 452 |
| MomentRT/fem.py                 |       92 |       11 |       16 |        1 |     87.04% |218-232, 284-294 |
| MomentRT/fields.py              |      100 |       26 |       30 |       11 |     66.92% |56, 59, 85-98, 113-122, 141, 237-242, 244, 319, 324, 411, 418, 456, 490 |
| MomentRT/geometry.py            |       64 |        3 |        6 |        3 |     91.43% |30, 33, 280 |
| MomentRT/material.py            |       92 |        7 |       38 |        6 |     88.46% |134, 163, 174, 303, 348, 367-368 |
| MomentRT/materials/liquids.py   |       64 |       51 |       14 |        0 |     16.67% |50-76, 88-120, 124, 181-268, 272 |
| MomentRT/materials/skin.py      |       46 |       36 |       12 |        0 |     17.24% |57-109, 113, 169-217, 221 |
| MomentRT/medium.py              |      125 |        6 |       50 |        8 |     92.00% |165, 184, 268, 332, 335, 501, 514-\>513, 525-\>528 |
| MomentRT/mesh.py                |      453 |      110 |      132 |       15 |     75.90% |60, 105, 116, 170-269, 286-320, 372, 438, 496, 514, 548, 570, 581, 639, 666, 684, 704, 823, 846, 855, 883-884, 1065, 1083-\>1086, 1125-1178, 1195-1207, 1236-1239, 1251 |
| MomentRT/phase\_function.py     |       82 |        9 |       26 |        8 |     84.26% |13, 22, 53, 69, 97, 121, 177, 267, 321, 364-\>369 |
| MomentRT/plotting.py            |      236 |       15 |       86 |       17 |     90.06% |80, 215, 263, 466, 469, 474, 481, 484, 491, 523, 597, 629-\>632, 697, 908, 953, 1010-\>1017, 1054 |
| MomentRT/results.py             |      207 |       31 |       54 |        7 |     80.08% |219-\>222, 224-\>227, 246, 260-\>278, 281-\>286, 680-\>exit, 789-910, 976, 1046-1075 |
| MomentRT/simulation.py          |      310 |       66 |      164 |       36 |     75.11% |100, 105-113, 125, 192, 197, 200, 203, 206, 209, 214, 217, 219-\>234, 231-232, 237-\>252, 246, 249-250, 259, 282, 312-317, 319-\>347, 335-336, 347-\>522, 365, 389-416, 432, 442-445, 450-\>484, 465, 481-482, 507-\>484, 527-534, 647-652, 660-663, 780, 785, 788, 848, 899, 904, 948, 979 |
| MomentRT/source.py              |       99 |        9 |       18 |        7 |     84.62% |76, 79, 174, 185, 188, 524, 538, 568-569 |
| MomentRT/transport/diffusion.py |      154 |       16 |       50 |       10 |     85.29% |63, 68, 124, 201, 308, 311, 314, 360-379, 549, 556 |
| MomentRT/transport/p3.py        |      152 |        2 |       38 |        2 |     97.89% |  379, 384 |
| MomentRT/transport/sp3.py       |       80 |        2 |       14 |        2 |     95.74% |  209, 216 |
| MomentRT/validation.py          |       20 |        0 |        2 |        1 |     95.45% |   82-\>85 |
| **TOTAL**                       | **2697** |  **414** |  **816** |  **147** | **81.53%** |           |

6 files skipped due to complete coverage.


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