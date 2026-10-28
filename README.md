# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/MomentRT/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                   |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|----------------------- | -------: | -------: | -------: | -------: | ---------: | --------: |
| MomentRT/detector.py   |       78 |        6 |       26 |        6 |     88.46% |21, 39, 48, 98, 127, 148 |
| MomentRT/diffusion.py  |      121 |        7 |       38 |        7 |     91.19% |25, 84, 108, 111, 116, 269, 276 |
| MomentRT/fields.py     |       79 |       20 |       20 |        7 |     68.69% |11, 14, 25-38, 41-48, 96, 101, 156, 163, 181 |
| MomentRT/geometry.py   |       48 |        2 |        4 |        2 |     92.31% |    11, 14 |
| MomentRT/material.py   |       44 |        2 |       12 |        2 |     92.86% |    28, 58 |
| MomentRT/medium.py     |      124 |        5 |       42 |        6 |     93.37% |55, 66, 116, 166, 169, 260-\>259 |
| MomentRT/mesh.py       |       23 |        4 |        2 |        1 |     80.00% | 14, 39-43 |
| MomentRT/plotting.py   |      170 |       14 |       64 |       12 |     88.89% |18-25, 122, 129, 132, 137, 144, 147, 154, 170, 174-\>179, 246-\>249, 444, 499-\>506 |
| MomentRT/results.py    |       54 |       13 |       12 |        1 |     69.70% |90, 117-148 |
| MomentRT/simulation.py |       87 |       25 |       30 |        5 |     72.65% |15, 25, 122-131, 162-181, 209-236 |
| MomentRT/source.py     |       26 |        1 |        2 |        1 |     92.86% |        37 |
| MomentRT/units.py      |       15 |        1 |       10 |        1 |     92.00% |        19 |
| **TOTAL**              | **1000** |  **100** |  **282** |   **51** | **87.29%** |           |

4 files skipped due to complete coverage.


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