# Repository Coverage

[Full report](https://htmlpreview.github.io/?https://github.com/MartinPdeS/MomentRT/blob/python-coverage-comment-action-data/htmlcov/index.html)

| Name                            |    Stmts |     Miss |   Branch |   BrPart |      Cover |   Missing |
|-------------------------------- | -------: | -------: | -------: | -------: | ---------: | --------: |
| MomentRT/boundary.py            |       71 |        5 |       16 |        4 |     89.66% |105, 111, 124, 170, 173 |
| MomentRT/detector.py            |      116 |        9 |       32 |        9 |     87.84% |86, 108, 117, 134, 137, 195, 298, 405, 452 |
| MomentRT/fem.py                 |       92 |       11 |       16 |        1 |     87.04% |218-232, 284-294 |
| MomentRT/fields.py              |       78 |       20 |       20 |        8 |     67.35% |54, 57, 82-95, 109-114, 247, 252, 339, 346, 384, 418 |
| MomentRT/geometry.py            |       64 |       13 |        6 |        2 |     75.71% |30, 33, 276-294 |
| MomentRT/material.py            |       92 |        7 |       38 |        6 |     88.46% |134, 163, 174, 303, 348, 367-368 |
| MomentRT/materials/liquids.py   |       64 |       51 |       14 |        0 |     16.67% |50-76, 88-120, 124, 181-268, 272 |
| MomentRT/materials/skin.py      |       46 |       36 |       12 |        0 |     17.24% |57-109, 113, 169-217, 221 |
| MomentRT/medium.py              |      116 |        5 |       42 |        6 |     93.04% |165, 184, 268, 332, 335, 494-\>493 |
| MomentRT/mesh.py                |       24 |        2 |        2 |        1 |     88.46% |   53, 105 |
| MomentRT/phase\_function.py     |       82 |        9 |       26 |        8 |     84.26% |13, 22, 53, 69, 97, 121, 177, 267, 321, 364-\>369 |
| MomentRT/plotting.py            |      225 |       13 |       82 |       15 |     90.88% |76, 211, 259, 460, 463, 468, 475, 478, 485, 583-\>586, 650, 861, 906, 963-\>970, 1007 |
| MomentRT/results.py             |      207 |       31 |       54 |        7 |     80.08% |219-\>222, 224-\>227, 246, 260-\>278, 281-\>286, 680-\>exit, 789-910, 976, 1046-1075 |
| MomentRT/simulation.py          |      301 |       64 |      158 |       36 |     74.73% |100, 105-113, 125, 192, 197, 200, 203, 206, 209, 212, 214-\>229, 226-227, 232-\>247, 241, 244-245, 254, 256-\>272, 277, 307-312, 314-\>342, 330-331, 342-\>513, 363, 380-407, 423, 433-436, 438-\>513, 441-\>475, 456, 472-473, 498-\>475, 518-525, 638-643, 651-654, 771, 776, 779, 839, 890, 895, 957 |
| MomentRT/source.py              |       99 |       16 |       18 |        5 |     76.92% |76, 79, 174, 185, 188, 519-553, 568-569 |
| MomentRT/transport/diffusion.py |      151 |       15 |       46 |        9 |     85.79% |63, 121, 198, 305, 308, 311, 363-382, 552, 559 |
| MomentRT/transport/p3.py        |      152 |        2 |       38 |        2 |     97.89% |  379, 384 |
| MomentRT/transport/sp3.py       |       80 |        2 |       14 |        2 |     95.74% |  209, 216 |
| MomentRT/units.py               |       16 |        1 |       10 |        1 |     92.31% |        59 |
| MomentRT/validation.py          |       20 |        0 |        2 |        1 |     95.45% |   82-\>85 |
| **TOTAL**                       | **2213** |  **312** |  **654** |  **123** | **82.18%** |           |

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