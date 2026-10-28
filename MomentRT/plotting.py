"""Optional interactive Plotly figures for geometry and physical fields."""

from collections.abc import Mapping
from html import escape

import numpy as np

import MomentRT

try:
    import plotly.graph_objects as go
    from plotly.colors import qualitative
    from plotly.subplots import make_subplots
except ImportError:
    go = None
    qualitative = None
    make_subplots = None

try:
    import pyvista as pv
except ImportError:
    pv = None


def require_plotly():
    if go is None:
        raise ImportError('Install 3D support with: pip install "MomentRT[viz]"')

    return go


def require_pyvista():
    if pv is None:
        raise ImportError(
            'Install VTK export support with: pip install "MomentRT[vtk]"'
        )

    return pv


def make_figure(*, length_unit):
    go = require_plotly()

    figure = go.Figure()

    figure.update_layout(
        template="plotly_white",
        scene={
            "xaxis_title": f"x [{length_unit}]",
            "yaxis_title": f"y [{length_unit}]",
            "zaxis_title": f"z [{length_unit}]",
            "aspectmode": "data",
        },
        margin={"l": 0, "r": 0, "t": 40, "b": 0},
    )

    return figure


def plot_medium(*, medium, mesh, length_unit, show):
    go = require_plotly()

    points = mesh.points.to(length_unit).magnitude

    lower = points[0, 0, 0]

    upper = points[-1, -1, -1]

    vertices = np.array(
        [
            [x, y, z]
            for x in (lower[0], upper[0])
            for y in (lower[1], upper[1])
            for z in (lower[2], upper[2])
        ]
    )

    triangles = np.array(
        [
            [0, 1, 3],
            [0, 3, 2],
            [4, 6, 7],
            [4, 7, 5],
            [0, 4, 5],
            [0, 5, 1],
            [2, 3, 7],
            [2, 7, 6],
            [0, 2, 6],
            [0, 6, 4],
            [1, 5, 7],
            [1, 7, 3],
        ]
    )

    figure = make_figure(length_unit=length_unit)

    figure.add_trace(
        go.Mesh3d(
            x=vertices[:, 0],
            y=vertices[:, 1],
            z=vertices[:, 2],
            i=triangles[:, 0],
            j=triangles[:, 1],
            k=triangles[:, 2],
            color="lightsteelblue",
            opacity=0.15,
            name="Observation window",
            showlegend=True,
            hovertemplate="Observation window (infinite medium)<extra></extra>",
        )
    )

    edge_points = []

    for a in range(8):
        for b in range(a + 1, 8):
            if np.count_nonzero(vertices[a] != vertices[b]) == 1:
                edge_points.extend([vertices[a], vertices[b], [None, None, None]])

    edges = np.asarray(edge_points, dtype=object)

    figure.add_trace(
        go.Scatter3d(
            x=edges[:, 0],
            y=edges[:, 1],
            z=edges[:, 2],
            mode="lines",
            line={"color": "slategray", "width": 3},
            showlegend=False,
            hoverinfo="skip",
        )
    )

    if show:
        figure.show()

    return figure


def _add_overlays(*, figure, overlays, context):
    """Add each object's spatial traces using a shared plotting context."""

    if overlays is None:
        return

    if isinstance(overlays, Mapping):
        entries = tuple(overlays.items())

        if any(not isinstance(name, str) or not name.strip() for name, _ in entries):
            raise ValueError("Overlay labels must be nonempty strings.")
    else:
        if isinstance(overlays, (str, bytes)):
            raise TypeError("overlays must be an iterable of drawable objects.")

        try:
            entries = tuple((None, item) for item in overlays)
        except TypeError as error:
            raise TypeError(
                "overlays must be an iterable of drawable objects."
            ) from error

    for _, item in entries:
        if not callable(getattr(item, "_add_to_figure", None)):
            raise TypeError(
                "Each overlay must implement _add_to_figure(figure=..., context=...)."
            )

    group_offset = len(figure.data)

    for number, (name, item) in enumerate(entries):
        item._add_to_figure(
            figure=figure,
            context={
                **context,
                "name": None if name is None else escape(name),
                "legendgroup": f"momentrt-overlay-{group_offset}-{number}",
            },
        )


def _circle_points(*, position, radius, face):
    """Build a display-coordinate circle in an axis-aligned face plane."""
    angle = np.linspace(0, 2 * np.pi, 65)

    ring = np.repeat(position[None, :], len(angle), axis=0).astype(float)

    axes = [axis for axis in range(3) if axis != "xyz".index(face[0])]

    ring[:, axes[0]] += radius * np.cos(angle)

    ring[:, axes[1]] += radius * np.sin(angle)

    return ring


def plot_field(
    *, field, view, scale, length_unit, opacity, isosurface_count, show, overlays=None
):
    if view not in {"slices", "volume", "isosurfaces"}:
        raise ValueError("view must be slices, volume, or isosurfaces.")

    if scale not in {"linear", "log"}:
        raise ValueError("scale must be linear or log.")

    if opacity is not None and (
        not np.isscalar(opacity) or not np.isfinite(opacity) or not 0 <= opacity <= 1
    ):
        raise ValueError("opacity must be between zero and one.")

    if isosurface_count is not None and (
        isinstance(isosurface_count, (bool, np.bool_))
        or not isinstance(isosurface_count, (int, np.integer))
        or isosurface_count < 1
    ):
        raise ValueError("isosurface_count must be a positive integer.")

    if view == "slices" and isosurface_count is not None:
        raise ValueError(
            "isosurface_count applies only to volume and isosurfaces views."
        )

    raw = np.asarray(field.values.magnitude)

    if scale == "log" and np.any(raw <= 0):
        raise ValueError(
            "Log display requires strictly positive values; no clipping is applied."
        )

    go = require_plotly()

    points = field.mesh.points.to(length_unit).magnitude

    # Only display coordinates are transformed; physical values remain intact.
    displayed = np.log10(raw) if scale == "log" else raw

    minimum = float(displayed.min())

    maximum = float(displayed.max())

    if minimum == maximum:
        maximum = minimum + max(abs(minimum), 1) * 1e-6

    colorbar = {"title": {"text": f"{field.name} [{field.values.units:~}]"}}

    if scale == "log":
        ticks = np.linspace(minimum, maximum, 5)

        colorbar.update(tickvals=ticks, ticktext=[f"{10**tick:.3g}" for tick in ticks])

    figure = make_figure(length_unit=length_unit)

    figure.update_layout(title=f"{field.name} — {view} ({scale} scale)")

    common = {
        "colorscale": "Inferno",
        "colorbar": colorbar,
        "customdata": raw.ravel(),
        "hovertemplate": f"{field.name}: %{{customdata:.4g}} {field.values.units:~}<extra></extra>",
    }

    if view == "slices":
        for axis in range(3):
            selection = [slice(None)] * 3

            # Use an actual central grid plane; do not imply interpolation.
            selection[axis] = field.mesh.shape[axis] // 2

            selection = tuple(selection)

            plane = points[selection]

            options = {**common, "customdata": raw[selection], "showscale": axis == 0}

            figure.add_trace(
                go.Surface(
                    x=plane[..., 0],
                    y=plane[..., 1],
                    z=plane[..., 2],
                    surfacecolor=displayed[selection],
                    cmin=minimum,
                    cmax=maximum,
                    name=f"{'xyz'[axis]} slice",
                    **options,
                )
            )
    else:
        trace = go.Volume if view == "volume" else go.Isosurface

        default_opacity = 0.12 if view == "volume" else 0.5

        display_opacity = default_opacity if opacity is None else opacity

        display_surface_count = (
            (8 if view == "volume" else 6)
            if isosurface_count is None
            else isosurface_count
        )

        figure.add_trace(
            trace(
                x=points[..., 0].ravel(),
                y=points[..., 1].ravel(),
                z=points[..., 2].ravel(),
                value=displayed.ravel(),
                isomin=minimum,
                isomax=maximum,
                opacity=display_opacity,
                surface_count=display_surface_count,
                caps={"x_show": False, "y_show": False, "z_show": False},
                **common,
            )
        )

        if view == "volume":
            opacity_values = [step / 20 for step in range(21)]

            if not any(np.isclose(display_opacity, value) for value in opacity_values):
                opacity_values.append(display_opacity)

            opacity_values.sort()

            active_index = next(
                index
                for index, value in enumerate(opacity_values)
                if np.isclose(value, display_opacity)
            )

            figure.update_layout(
                sliders=[
                    {
                        "active": active_index,
                        "currentvalue": {
                            "font": {"size": 12},
                            "prefix": "Opacity: ",
                        },
                        "font": {"size": 9},
                        "len": 0.55,
                        "minorticklen": 0,
                        "pad": {"t": 12},
                        "ticklen": 2,
                        "x": 0.225,
                        "steps": [
                            {
                                "label": (
                                    f"{value:.2f}"
                                    if np.isclose(value, display_opacity)
                                    or np.isclose(value * 4, round(value * 4))
                                    else ""
                                ),
                                "method": "restyle",
                                "args": [{"opacity": value}, [0]],
                            }
                            for value in opacity_values
                        ],
                    }
                ]
            )

    field_trace_count = len(figure.data)

    _add_overlays(
        figure=figure,
        overlays=overlays,
        context={
            "length_unit": length_unit,
            "sample_origin": field.metadata.get("sample_origin"),
            "sample_size": field.metadata.get("sample_size"),
        },
    )

    if len(figure.data) > field_trace_count:
        figure.update_layout(
            legend={
                "orientation": "h",
                "x": 0.02,
                "xanchor": "left",
                "y": 0.98,
                "yanchor": "top",
                "groupclick": "togglegroup",
            }
        )

    if show:
        figure.show()

    return figure


def _box_surface(*, origin, size, name, color, opacity, hover):
    go = require_plotly()

    vertices = np.array(
        [
            [x, y, z]
            for x in (origin[0], origin[0] + size[0])
            for y in (origin[1], origin[1] + size[1])
            for z in (origin[2], origin[2] + size[2])
        ]
    )

    triangles = np.array(
        [
            [0, 1, 3],
            [0, 3, 2],
            [4, 6, 7],
            [4, 7, 5],
            [0, 4, 5],
            [0, 5, 1],
            [2, 3, 7],
            [2, 7, 6],
            [0, 2, 6],
            [0, 6, 4],
            [1, 5, 7],
            [1, 7, 3],
        ]
    )

    return go.Mesh3d(
        x=vertices[:, 0],
        y=vertices[:, 1],
        z=vertices[:, 2],
        i=triangles[:, 0],
        j=triangles[:, 1],
        k=triangles[:, 2],
        name=name,
        color=color,
        opacity=opacity,
        showlegend=True,
        hovertemplate=hover,
    )


def plot_sample(*, sample, length_unit, opacity, show):
    """Render construction surfaces with stable colors and material hover labels."""

    if not np.isscalar(opacity) or not np.isfinite(opacity) or not 0 <= opacity <= 1:
        raise ValueError("opacity must be between zero and one.")

    go = require_plotly()

    figure = make_figure(length_unit=length_unit)

    origin = sample.origin.to(length_unit).magnitude

    size = sample.size.to(length_unit).magnitude

    colors = ["#E5AD79", "#CB7582", "#EDCC71", "#65B3AD", "#7C93C3"]

    for region_number, region in enumerate((None, *sample.structures)):
        material = sample.background if region is None else region.material

        name = "Background / sample extent" if region is None else region.name

        # User-supplied names are text, not HTML instructions.
        hover = (
            f"{escape(name)}<br>Absorption: {material.absorption:~P}"
            f"<br>Reduced scattering: {material.reduced_scattering:~P}"
            f"<br>Refractive index: {material.refractive_index:g}"
            f"<br>Priority: {region_number} (later replaces earlier)<extra></extra>"
        )

        color = (
            "#A7B5C8" if region is None else colors[(region_number - 1) % len(colors)]
        )

        if region is None or isinstance(
            region, (MomentRT.geometry.Box, MomentRT.geometry.Layer)
        ):
            region_origin = origin.copy()

            region_size = size.copy()

            if isinstance(region, MomentRT.geometry.Box):
                region_origin = region.origin.to(length_unit).magnitude

                region_size = region.size.to(length_unit).magnitude
            elif isinstance(region, MomentRT.geometry.Layer):
                region_origin[2] = region.depth.to(length_unit).magnitude

                region_size[2] = region.thickness.to(length_unit).magnitude

            figure.add_trace(
                _box_surface(
                    origin=region_origin,
                    size=region_size,
                    name=name,
                    color=color,
                    opacity=min(opacity, 0.08) if region is None else opacity,
                    hover=hover,
                )
            )
        else:
            azimuth, polar = np.meshgrid(
                np.linspace(0, 2 * np.pi, 49), np.linspace(0, np.pi, 25)
            )

            center = region.center.to(length_unit).magnitude

            radius = region.radius.to(length_unit).magnitude

            figure.add_trace(
                go.Surface(
                    x=center[0] + radius * np.sin(polar) * np.cos(azimuth),
                    y=center[1] + radius * np.sin(polar) * np.sin(azimuth),
                    z=center[2] + radius * np.cos(polar),
                    colorscale=[[0, color], [1, color]],
                    showscale=False,
                    name=name,
                    opacity=opacity,
                    showlegend=True,
                    hovertemplate=hover,
                )
            )

    figure.update_layout(
        title="Sample geometry",
        scene_camera={"eye": {"x": 1.8, "y": 1.8, "z": 1.8}},
        scene_zaxis={"autorange": "reversed"},
        legend={"title": {"text": "Regions (later replaces earlier)"}},
    )

    if show:
        figure.show()

    return figure


def plot_section(*, section, quantity, length_unit, value_unit, scale, show):
    """Plot resolved material or property samples with physical hover values."""

    if quantity not in ("material", *section.fields):
        raise ValueError(
            "quantity must be material, absorption, reduced_scattering, or refractive_index."
        )

    if scale not in ("linear", "log"):
        raise ValueError("scale must be linear or log.")

    go = require_plotly()

    names = np.asarray([escape(name) for name in section.names])[section.labels]

    if quantity == "material":
        if scale != "linear" or value_unit is not None:
            raise ValueError(
                "Material labels require a linear scale and no value_unit."
            )

        colors = ["#A7B5C8", "#E5AD79", "#CB7582", "#EDCC71", "#65B3AD", "#7C93C3"]

        count = len(section.names)

        colorscale = []

        for region in range(count):
            color = colors[region % len(colors)]

            colorscale.extend([[region / count, color], [(region + 1) / count, color]])

        values = section.labels

        colorbar = {
            "title": {"text": "Region"},
            "tickvals": list(range(count)),
            "ticktext": [escape(name) for name in section.names],
        }

        limits = {"zmin": -0.5, "zmax": count - 0.5}

        hover = "%{customdata}<extra></extra>"

        customdata = names
    else:
        physical = section.fields[quantity]

        if value_unit is not None:
            physical = physical.to(value_unit)

        raw = physical.magnitude

        if scale == "log" and np.any(raw <= 0):
            raise ValueError(
                "Log display requires strictly positive values; no clipping is applied."
            )

        values = np.log10(raw) if scale == "log" else raw

        colorscale = "Inferno"

        colorbar = {"title": {"text": f"{quantity} [{physical.units:~}]"}}

        if scale == "log":
            ticks = np.linspace(values.min(), values.max(), 5)

            colorbar.update(
                tickvals=ticks, ticktext=[f"{10**tick:.3g}" for tick in ticks]
            )

        limits = {}

        customdata = np.stack([names, raw.astype(object)], axis=-1)

        hover = f"%{{customdata[0]}}<br>{quantity}: %{{customdata[1]:.4g}} {physical.units:~}<extra></extra>"

    figure = go.Figure(
        go.Heatmap(
            x=section.horizontal.to(length_unit).magnitude,
            y=section.vertical.to(length_unit).magnitude,
            z=values,
            colorscale=colorscale,
            colorbar=colorbar,
            customdata=customdata,
            hovertemplate=hover,
            zsmooth=False,
            **limits,
        )
    )

    extent = section.extent.to(length_unit).magnitude

    vertical_range = extent[1, ::-1] if section.axes[1] == "z" else extent[1]

    figure.update_layout(
        template="plotly_white",
        title=f"{quantity} — {section.normal} = {section.position.to(length_unit):~P}",
        xaxis={"title": f"{section.axes[0]} [{length_unit}]", "range": extent[0]},
        yaxis={
            "title": f"{section.axes[1]} [{length_unit}]",
            "range": vertical_range,
            "scaleanchor": "x",
            "scaleratio": 1,
        },
    )

    if show:
        figure.show()

    return figure
