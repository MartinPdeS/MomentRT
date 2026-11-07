# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/MomentRT/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                                           |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|----------------------------------------------- | -------: | -------: | -------: | -------: | ---------: | --------: |
| MomentRT/\_face\_geometry.py                   |       32 |        2 |        8 |        2 |     90.00% |    34, 37 |
| MomentRT/boundary.py                           |       43 |        3 |        8 |        2 |     90.20% |104, 110, 123 |
| MomentRT/convergence.py                        |      267 |       24 |      128 |       25 |     87.59% |107, 153, 155, 168, 174, 176, 240, 245, 266, 289, 297, 392, 397, 406, 409, 417, 514-\>534, 534-\>537, 542, 596-\>603, 643-646, 668, 673, 722-723, 757 |
| MomentRT/detector.py                           |      155 |       13 |       50 |       13 |     87.32% |49, 51, 146, 226, 258, 267, 284, 287, 292, 361, 478, 585, 632 |
| MomentRT/fem/kernels.py                        |       70 |       11 |       14 |        1 |     83.33% |217-231, 283-293 |
| MomentRT/fields.py                             |       85 |       10 |       26 |       10 |     81.98% |56, 59, 86, 189, 264, 271, 376, 383, 421, 455 |
| MomentRT/geometry.py                           |       65 |        3 |        6 |        3 |     91.55% |31, 34, 310 |
| MomentRT/material.py                           |       93 |        7 |       38 |        6 |     88.55% |141, 170, 181, 319, 366, 385-386 |
| MomentRT/materials/liquids.py                  |       65 |       43 |       14 |        1 |     29.11% |68, 90-123, 127, 184-272, 276 |
| MomentRT/materials/skin.py                     |       47 |        8 |       12 |        6 |     76.27% |64, 74, 82, 115, 177, 187, 195, 224 |
| MomentRT/medium.py                             |      126 |        6 |       50 |        7 |     92.61% |177, 196, 280, 346, 349, 515, 529-\>528 |
| MomentRT/mesh.py                               |      444 |       53 |      132 |       18 |     86.28% |61, 182-279, 313, 423, 483, 501, 535, 653, 671, 691, 812, 835, 844, 872-873, 1054, 1072-\>1075, 1119, 1165, 1214 |
| MomentRT/models/spectrum.py                    |       31 |        1 |        8 |        1 |     94.87% |        29 |
| MomentRT/phase\_function/base.py               |       19 |        4 |        4 |        2 |     73.91% |13, 22, 53, 69 |
| MomentRT/phase\_function/henyey\_greenstein.py |       17 |        2 |        4 |        2 |     80.95% |    35, 59 |
| MomentRT/phase\_function/isotropic.py          |       13 |        1 |        2 |        1 |     86.67% |        57 |
| MomentRT/phase\_function/tabulated.py          |       39 |        2 |       14 |        3 |     90.57% |41, 95, 138-\>143 |
| MomentRT/plotting.py                           |      209 |       14 |       82 |       15 |     90.03% |25, 156, 204, 403, 406, 411, 418, 421, 428, 532, 564-\>567, 632, 839, 884, 939-\>946, 983 |
| MomentRT/results/convergence.py                |      289 |       17 |      100 |       23 |     88.69% |318-\>exit, 382, 446, 450, 460-465, 467-\>470, 472, 482, 483-\>489, 492, 514, 518-519, 620-\>638, 764, 766, 788, 889, 929-\>910, 959-\>exit, 960-\>959, 1088-\>1093, 1148-\>1152, 1155 |
| MomentRT/results/simulation.py                 |       67 |        1 |       22 |        5 |     93.26% |218-\>221, 223-\>226, 245, 259-\>277, 280-\>285 |
| MomentRT/results/spectral.py                   |       91 |       19 |       30 |        7 |     73.55% |75, 78, 97, 106, 124, 199, 230, 301, 371-400 |
| MomentRT/results/transient.py                  |       29 |        3 |        8 |        3 |     83.78% |75, 81, 120 |
| MomentRT/simulation.py                         |      311 |       65 |      176 |       33 |     76.59% |105, 110-118, 130, 214, 219, 222, 225, 228, 231, 236, 239, 241-\>257, 254-255, 260-\>279, 273, 276-277, 294, 317, 347-352, 370-371, 382-\>553, 400, 424-451, 467, 477-480, 485-\>515, 496, 512-513, 558-565, 595, 597, 602, 654, 765-770, 778-781, 862, 975, 1006 |
| MomentRT/source.py                             |      166 |       16 |       50 |       14 |     85.19% |50, 62, 111, 115, 170, 187, 194, 336, 339, 437, 446, 449, 803, 822, 852-853 |
| MomentRT/transport/infinite\_p1.py             |       22 |        2 |        8 |        2 |     86.67% |    64, 69 |
| MomentRT/transport/p1.py                       |      138 |       14 |       42 |        8 |     85.56% |33, 117, 224, 227, 230, 277-296, 463, 470 |
| MomentRT/transport/p3.py                       |      155 |        2 |       38 |        2 |     97.93% |  388, 393 |
| MomentRT/transport/sp3.py                      |       83 |        2 |       14 |        2 |     95.88% |  219, 226 |
| MomentRT/transport/transient\_p1.py            |      109 |        7 |       28 |        7 |     89.78% |84, 89, 139, 141, 143, 205, 214 |
| MomentRT/validation.py                         |       22 |        0 |        2 |        1 |     95.83% |   92-\>95 |
| **TOTAL**                                      | **3571** |  **355** | **1138** |  **225** | **86.49%** |           |

15 files skipped due to complete coverage.


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