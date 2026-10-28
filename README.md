# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/MomentRT/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                   |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|----------------------- | -------: | -------: | -------: | -------: | ---------: | --------: |
| MomentRT/diffusion.py  |      110 |        7 |       34 |        7 |     90.28% |25, 82, 104, 107, 110, 263, 270 |
| MomentRT/fields.py     |       55 |       19 |       16 |        6 |     59.15% |11, 14, 25-38, 41-48, 78, 85, 103, 137 |
| MomentRT/geometry.py   |       48 |        2 |        4 |        2 |     92.31% |    11, 14 |
| MomentRT/material.py   |       11 |        1 |        2 |        1 |     84.62% |        28 |
| MomentRT/medium.py     |       96 |        4 |       30 |        5 |     92.86% |52, 70, 120, 123, 204-\>203 |
| MomentRT/mesh.py       |       23 |        4 |        2 |        1 |     80.00% | 14, 39-43 |
| MomentRT/plotting.py   |      154 |       16 |       54 |       12 |     86.54% |18-25, 122, 129, 132, 137, 153, 157-\>162, 217, 308-310, 357, 372, 427-\>434, 471 |
| MomentRT/simulation.py |       44 |       14 |       12 |        3 |     66.07% |38-47, 78-97, 125 |
| MomentRT/source.py     |       26 |        1 |        2 |        1 |     92.86% |        37 |
| MomentRT/units.py      |       15 |        3 |       10 |        3 |     76.00% |16, 19, 25 |
| **TOTAL**              |  **710** |   **71** |  **186** |   **41** | **86.83%** |           |

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