"""Embed each new Plotly figure after the gallery cell that constructs it."""

import html
import os
from pathlib import Path


def plotly_scraper(block, block_vars, gallery_conf):
    from plotly.graph_objects import Figure

    rendered = block_vars.setdefault("momentrt_figures", [])

    output = Path(gallery_conf["src_dir"]) / "_static" / "gallery"

    output.mkdir(parents=True, exist_ok=True)

    script = Path(block_vars["target_file"])

    markup = []

    for figure in block_vars["example_globals"].values():
        if not isinstance(figure, Figure) or any(
            figure is previous for previous in rendered
        ):
            continue

        rendered.append(figure)

        filename = output / f"{script.stem}-{len(rendered)}.html"

        figure.write_html(
            file=str(filename),
            include_plotlyjs="directory",
            config={"responsive": True, "scrollZoom": True},
        )

        relative = Path(os.path.relpath(filename, script.parent)).as_posix()

        title = html.escape(figure.layout.title.text or "Interactive view", quote=True)

        markup.append(
            "\n.. raw:: html\n\n"
            f'   <iframe src="{relative}" title="{title}" width="100%" '
            'height="600" loading="lazy" style="border:0;"></iframe>\n'
        )

    return "".join(markup)
