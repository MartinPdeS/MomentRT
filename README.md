# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/MomentRT/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                        |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|---------------------------- | -------: | -------: | -------: | -------: | ---------: | --------: |
| MomentRT/boundary.py        |       71 |        5 |       16 |        4 |     89.66% |108, 114, 127, 173, 176 |
| MomentRT/convergence.py     |       69 |        5 |       32 |        5 |     90.10% |56, 61, 64, 167, 172 |
| MomentRT/detector.py        |      110 |        7 |       28 |        7 |     89.86% |79, 101, 110, 172, 265, 372, 419 |
| MomentRT/diffusion.py       |      143 |        8 |       40 |        8 |     91.26% |66, 124, 201, 308, 311, 314, 531, 538 |
| MomentRT/fem.py             |       84 |        4 |       12 |        0 |     95.83% |   270-280 |
| MomentRT/fields.py          |       78 |       20 |       20 |        8 |     67.35% |57, 60, 85-98, 112-117, 248, 253, 340, 347, 385, 419 |
| MomentRT/geometry.py        |       50 |        2 |        4 |        2 |     92.59% |    30, 33 |
| MomentRT/material.py        |       91 |        7 |       38 |        6 |     88.37% |115, 144, 155, 284, 329, 345-346 |
| MomentRT/medium.py          |      116 |        5 |       42 |        6 |     93.04% |168, 187, 271, 335, 338, 497-\>496 |
| MomentRT/mesh.py            |       24 |        2 |        2 |        1 |     88.46% |   55, 107 |
| MomentRT/p3.py              |      142 |        2 |       34 |        2 |     97.73% |  350, 355 |
| MomentRT/phase\_function.py |       82 |        9 |       26 |        8 |     84.26% |16, 25, 56, 72, 100, 122, 171, 256, 310, 350-\>355 |
| MomentRT/plotting.py        |      225 |       13 |       84 |       15 |     90.94% |79, 214, 245, 461, 464, 469, 476, 479, 486, 584-\>587, 651, 862, 907, 964-\>971, 1008 |
| MomentRT/results.py         |      207 |       31 |       54 |        7 |     80.08% |222-\>225, 227-\>230, 249, 263-\>281, 284-\>289, 683-\>exit, 792-913, 979, 1049-1078 |
| MomentRT/simulation.py      |       69 |        3 |       20 |        3 |     93.26% |96, 106, 366 |
| MomentRT/source.py          |       72 |        1 |        4 |        1 |     97.37% |       142 |
| MomentRT/sp3.py             |       80 |        2 |       14 |        2 |     95.74% |  212, 219 |
| MomentRT/units.py           |       16 |        1 |       10 |        1 |     92.31% |        60 |
| MomentRT/validation.py      |      159 |       31 |       84 |       23 |     75.31% |87-\>90, 125, 130, 133, 136, 139, 142, 145, 147-\>162, 159-160, 165-\>180, 174, 177-178, 187, 189-\>203, 208, 233-238, 240-\>261, 251-252, 261-\>389, 281, 311, 319-322, 324-\>389, 327-\>356, 353-354, 377-\>356, 394-401 |
| **TOTAL**                   | **1992** |  **158** |  **572** |  **109** | **88.57%** |           |

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