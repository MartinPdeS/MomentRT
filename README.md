# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/MomentRT/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                            |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|-------------------------------- | -------: | -------: | -------: | -------: | ---------: | --------: |
| MomentRT/boundary.py            |       71 |        5 |       16 |        4 |     89.66% |105, 111, 124, 170, 173 |
| MomentRT/convergence.py         |      193 |       18 |      104 |       21 |     86.87% |31, 33, 46, 52, 54, 113, 118, 139, 162, 170, 259, 264, 273, 276, 284, 368-\>387, 389-\>392, 402-\>410, 415, 466, 471 |
| MomentRT/detector.py            |      154 |       13 |       50 |       13 |     87.25% |44, 46, 136, 216, 238, 247, 264, 267, 272, 341, 458, 565, 612 |
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
| MomentRT/results.py             |      371 |       43 |      104 |       20 |     83.37% |221-\>224, 226-\>229, 248, 262-\>280, 283-\>288, 699-\>exit, 760-763, 829-950, 1049, 1089-\>1070, 1119-\>exit, 1120-\>1119, 1248-\>1253, 1308-\>1312, 1315, 1377, 1380, 1398, 1407, 1425, 1499, 1530, 1596, 1666-1695 |
| MomentRT/simulation.py          |      288 |       60 |      154 |       29 |     76.24% |100, 105-113, 125, 202, 207, 210, 213, 216, 219, 224, 227, 229-\>244, 241-242, 247-\>262, 256, 259-260, 277, 300, 330-335, 337-\>365, 353-354, 365-\>540, 383, 407-434, 450, 460-463, 468-\>502, 483, 499-500, 545-552, 711-716, 724-727, 916, 947 |
| MomentRT/source.py              |      165 |       16 |       50 |       14 |     85.12% |45, 56, 104, 108, 163, 180, 187, 319, 322, 420, 429, 432, 777, 796, 826-827 |
| MomentRT/transport/diffusion.py |      156 |       16 |       50 |       10 |     85.44% |69, 74, 130, 214, 321, 324, 327, 373-392, 562, 569 |
| MomentRT/transport/p3.py        |      153 |        2 |       38 |        2 |     97.91% |  385, 390 |
| MomentRT/transport/sp3.py       |       81 |        2 |       14 |        2 |     95.79% |  216, 223 |
| MomentRT/validation.py          |       22 |        0 |        2 |        1 |     95.83% |   92-\>95 |
| **TOTAL**                       | **3096** |  **411** | **1002** |  **184** | **83.53%** |           |

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