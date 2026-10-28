# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/MomentRT/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                    |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|------------------------ | -------: | -------: | -------: | -------: | ---------: | --------: |
| MomentRT/convergence.py |       65 |        5 |       32 |        5 |     89.69% |16, 19, 22, 125, 130 |
| MomentRT/detector.py    |       92 |        6 |       26 |        6 |     89.83% |25, 47, 56, 105, 186, 205 |
| MomentRT/diffusion.py   |      137 |        8 |       40 |        8 |     90.96% |26, 68, 106, 129, 132, 135, 317, 324 |
| MomentRT/fem.py         |       81 |        4 |       12 |        0 |     95.70% |   134-144 |
| MomentRT/fields.py      |       74 |       19 |       20 |        7 |     68.09% |13, 16, 27-40, 43-48, 100, 105, 160, 167, 183 |
| MomentRT/geometry.py    |       47 |        2 |        4 |        2 |     92.16% |    12, 15 |
| MomentRT/material.py    |       44 |        2 |       12 |        2 |     92.86% |    28, 60 |
| MomentRT/medium.py      |      111 |        5 |       42 |        6 |     92.81% |60, 68, 123, 170, 173, 258-\>257 |
| MomentRT/mesh.py        |       22 |        2 |        2 |        1 |     87.50% |    14, 44 |
| MomentRT/plotting.py    |      208 |       10 |       84 |       12 |     92.47% |38, 135, 150, 200, 203, 208, 215, 218, 225, 317-\>320, 537, 592-\>599 |
| MomentRT/results.py     |      150 |       30 |       32 |        2 |     74.73% |277-\>exit, 336-457, 481, 507-536 |
| MomentRT/simulation.py  |       62 |        2 |       20 |        2 |     95.12% |    15, 25 |
| MomentRT/source.py      |       55 |        1 |        4 |        1 |     96.61% |        60 |
| MomentRT/units.py       |       15 |        1 |       10 |        1 |     92.00% |        19 |
| **TOTAL**               | **1269** |   **97** |  **350** |   **55** | **89.50%** |           |

3 files skipped due to complete coverage.


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