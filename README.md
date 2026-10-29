# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/MomentRT/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                    |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|------------------------ | -------: | -------: | -------: | -------: | ---------: | --------: |
| MomentRT/boundary.py    |       66 |        5 |       16 |        4 |     89.02% |108, 114, 127, 173, 176 |
| MomentRT/convergence.py |       69 |        5 |       32 |        5 |     90.10% |55, 60, 63, 166, 171 |
| MomentRT/detector.py    |      107 |        7 |       28 |        7 |     89.63% |79, 101, 110, 172, 223, 328, 368 |
| MomentRT/diffusion.py   |      140 |        8 |       40 |        8 |     91.11% |66, 124, 193, 300, 303, 306, 512, 519 |
| MomentRT/fem.py         |       84 |        4 |       12 |        0 |     95.83% |   270-280 |
| MomentRT/fields.py      |       78 |       19 |       20 |        7 |     69.39% |57, 60, 85-98, 112-117, 248, 253, 340, 347, 385 |
| MomentRT/geometry.py    |       50 |        2 |        4 |        2 |     92.59% |    30, 33 |
| MomentRT/material.py    |       64 |        3 |       20 |        3 |     92.86% |114, 218, 263 |
| MomentRT/medium.py      |      115 |        5 |       42 |        6 |     92.99% |167, 186, 270, 334, 337, 489-\>488 |
| MomentRT/mesh.py        |       24 |        2 |        2 |        1 |     88.46% |   55, 107 |
| MomentRT/plotting.py    |      215 |       10 |       84 |       12 |     92.64% |79, 214, 245, 351, 354, 359, 366, 369, 376, 474-\>477, 783, 840-\>847 |
| MomentRT/results.py     |      188 |       31 |       48 |        8 |     77.54% |174, 183-\>168, 188-\>206, 201-\>206, 209-\>214, 214-\>219, 608-\>exit, 717-838, 904, 958-987 |
| MomentRT/simulation.py  |       69 |        2 |       20 |        2 |     95.51% |   96, 106 |
| MomentRT/source.py      |       65 |        1 |        4 |        1 |     97.10% |       142 |
| MomentRT/sp3.py         |       78 |        3 |       16 |        3 |     93.62% |144, 214, 221 |
| MomentRT/units.py       |       16 |        1 |       10 |        1 |     92.31% |        60 |
| MomentRT/validation.py  |      158 |       33 |       84 |       28 |     72.31% |87-\>90, 123, 128, 131, 134, 137, 140, 143, 145-\>160, 157-158, 163-\>178, 172, 175-176, 185, 187-\>200, 196, 205, 223-231, 233-\>247, 237-\>242, 239-240, 247-\>375, 267, 283-\>310, 297, 305-308, 310-\>375, 313-\>342, 335-\>342, 339-340, 363-\>342, 380-387 |
| **TOTAL**               | **1682** |  **141** |  **490** |   **98** | **87.89%** |           |

2 files skipped due to complete coverage.


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