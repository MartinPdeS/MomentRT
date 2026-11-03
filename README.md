# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/MomentRT/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                               |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|----------------------------------- | -------: | -------: | -------: | -------: | ---------: | --------: |
| MomentRT/boundary.py               |       71 |        5 |       16 |        4 |     89.66% |105, 111, 124, 170, 173 |
| MomentRT/convergence.py            |      266 |       24 |      128 |       25 |     87.56% |106, 152, 154, 167, 173, 175, 234, 239, 260, 283, 291, 386, 391, 400, 403, 411, 503-\>522, 522-\>525, 530, 584-\>591, 631-634, 656, 661, 710-711, 745 |
| MomentRT/detector.py               |      154 |       13 |       50 |       13 |     87.25% |44, 46, 136, 216, 238, 247, 264, 267, 272, 341, 458, 565, 612 |
| MomentRT/fem.py                    |       92 |       11 |       16 |        1 |     87.04% |218-232, 284-294 |
| MomentRT/fields.py                 |       84 |       10 |       26 |       10 |     81.82% |55, 58, 85, 188, 263, 268, 355, 362, 400, 434 |
| MomentRT/geometry.py               |       64 |        3 |        6 |        3 |     91.43% |30, 33, 280 |
| MomentRT/material.py               |       92 |        7 |       38 |        6 |     88.46% |134, 163, 174, 303, 348, 367-368 |
| MomentRT/materials/liquids.py      |       64 |       51 |       14 |        0 |     16.67% |50-76, 88-120, 124, 181-268, 272 |
| MomentRT/materials/skin.py         |       46 |       36 |       12 |        0 |     17.24% |57-109, 113, 169-217, 221 |
| MomentRT/medium.py                 |      125 |        6 |       50 |        8 |     92.00% |165, 184, 268, 332, 335, 501, 514-\>513, 525-\>528 |
| MomentRT/mesh.py                   |      443 |       53 |      132 |       18 |     86.26% |60, 170-267, 301, 410, 468, 486, 520, 638, 656, 676, 795, 818, 827, 855-856, 1037, 1055-\>1058, 1102, 1148, 1195 |
| MomentRT/phase\_function.py        |       82 |        9 |       26 |        8 |     84.26% |13, 22, 53, 69, 97, 121, 177, 267, 321, 364-\>369 |
| MomentRT/plotting.py               |      209 |       14 |       82 |       15 |     90.03% |25, 156, 204, 403, 406, 411, 418, 421, 428, 532, 564-\>567, 632, 839, 884, 939-\>946, 983 |
| MomentRT/results/convergence.py    |      289 |       17 |      100 |       23 |     88.69% |318-\>exit, 382, 446, 450, 460-465, 467-\>470, 472, 482, 483-\>489, 492, 514, 518-519, 620-\>638, 764, 766, 788, 889, 929-\>910, 959-\>exit, 960-\>959, 1088-\>1093, 1148-\>1152, 1155 |
| MomentRT/results/simulation.py     |       67 |        1 |       22 |        5 |     93.26% |218-\>221, 223-\>226, 245, 259-\>277, 280-\>285 |
| MomentRT/results/spectral.py       |       90 |       19 |       30 |        7 |     73.33% |70, 73, 91, 100, 118, 192, 223, 290, 360-389 |
| MomentRT/simulation.py             |      288 |       60 |      154 |       29 |     76.24% |100, 105-113, 125, 202, 207, 210, 213, 216, 219, 224, 227, 229-\>244, 241-242, 247-\>262, 256, 259-260, 277, 300, 330-335, 337-\>365, 353-354, 365-\>540, 383, 407-434, 450, 460-463, 468-\>502, 483, 499-500, 545-552, 711-716, 724-727, 916, 947 |
| MomentRT/source.py                 |      165 |       16 |       50 |       14 |     85.12% |45, 56, 104, 108, 163, 180, 187, 319, 322, 420, 429, 432, 777, 796, 826-827 |
| MomentRT/transport/infinite\_p1.py |       22 |        2 |        8 |        2 |     86.67% |    64, 69 |
| MomentRT/transport/p1.py           |      137 |       14 |       42 |        8 |     85.47% |32, 116, 223, 226, 229, 275-294, 464, 471 |
| MomentRT/transport/p3.py           |      153 |        2 |       38 |        2 |     97.91% |  385, 390 |
| MomentRT/transport/sp3.py          |       81 |        2 |       14 |        2 |     95.79% |  216, 223 |
| MomentRT/validation.py             |       22 |        0 |        2 |        1 |     95.83% |   92-\>95 |
| **TOTAL**                          | **3274** |  **375** | **1074** |  **204** | **85.07%** |           |

8 files skipped due to complete coverage.


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