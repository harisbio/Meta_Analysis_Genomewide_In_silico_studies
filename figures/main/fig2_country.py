import sys as _sys
from pathlib import Path as _P
_sys.path.insert(0, str(_P(__file__).resolve().parents[2]))
from config import DATA_FILE, GLOBAL_DATA_FILE, OUTPUT_DIR  # noqa: E402
import math
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

DATA = str(DATA_FILE)
OUT = Path(OUTPUT_DIR / "Figure_2_Country_Distribution")
d = pd.read_csv(DATA, dtype=str, keep_default_na=False)
N = len(d)
counts = d["Country of Origin"].value_counts()
order = sorted(counts.index, key=lambda c: (-counts[c], c))

INK, MUTED, GRID = "#16324f", "#5f7186", "#e6ebf0"
RAMP = [[0, "#a9d2f0"], [0.3, "#7fb8e6"], [0.6, "#3a8fd0"], [0.85, "#0a5fa3"], [1, "#06305c"]]
BAR = "#0a5fa3"
FONT = "Arial"

fig = make_subplots(rows=2, cols=1, specs=[[{"type": "choropleth"}], [{"type": "xy"}]], row_heights=[0.57, 0.43], vertical_spacing=0.05)
z = [math.log10(counts[c]) for c in order]
ticks = [1, 2, 5, 10, 25, 50, 100, int(counts.max())]
fig.add_trace(go.Choropleth(locations=order, locationmode="country names", z=z, zmin=0, zmax=math.log10(counts.max()), colorscale=RAMP,
                            marker_line_color="#26323f", marker_line_width=0.8, customdata=[[counts[c], counts[c] / N * 100] for c in order],
                            hovertemplate="<b>%{location}</b><br>%{customdata[0]} papers (%{customdata[1]:.1f}%)<extra></extra>",
                            colorbar=dict(
                                          tickvals=[math.log10(t) for t in ticks], ticktext=[str(t) for t in ticks], orientation="v",
                                          x=0.012, xanchor="left", y=0.76, yanchor="middle", len=0.36, thickness=22,
                                          tickfont=dict(size=18, color=INK), outlinewidth=0)),
              row=1, col=1)
fig.update_geos(projection_type="natural earth", showcountries=True, countrycolor="#9aa4af", countrywidth=0.5, showland=True, landcolor="#e3e6ea",
                showocean=False, showlakes=False, showcoastlines=False, showframe=False, lataxis_range=[-56, 84], domain_x=[0.055, 1], projection_scale=1.04, center=dict(lon=12, lat=10), bgcolor="rgba(0,0,0,0)")

ys = order[::-1]
fig.add_trace(go.Bar(y=ys, x=[counts[c] for c in ys], orientation="h", marker_color=BAR, marker_line_width=0, width=0.68, showlegend=False,
                     text=[f"{counts[c]}  ({counts[c] / N * 100:.1f}%)" for c in ys], textposition="outside", cliponaxis=False,
                     textfont=dict(size=19, color=INK), hovertemplate="%{y}: %{x} papers<extra></extra>"), row=2, col=1)
fig.update_xaxes(row=2, col=1, range=[0, counts.max() * 1.16], showgrid=True, gridcolor=GRID, zeroline=False, tickfont=dict(size=18, color=INK),
                 title=dict(text="Number of papers", font=dict(size=21, color=INK)))
fig.update_yaxes(row=2, col=1, tickfont=dict(size=19, color=INK), showgrid=False, ticks="", automargin=True, ticksuffix="   ")

fig.add_annotation(xref="paper", yref="paper", x=0.004, y=0.76, xanchor="center", yanchor="middle", textangle=-90, showarrow=False, text="Number of papers (log scale)", font=dict(size=20, color=INK, family=FONT))
for letter, y in (("A", 1.0), ("B", 0.435)):
    fig.add_annotation(xref="paper", yref="paper", x=0, y=y, xanchor="left", yanchor="bottom", showarrow=False, text=f"<b>{letter}</b>", font=dict(size=32, color=INK, family=FONT))
fig.update_layout(width=1800, height=1560, paper_bgcolor="white", plot_bgcolor="white", font=dict(family=FONT, color=INK), margin=dict(l=20, r=30, t=55, b=30), showlegend=False)
OUT.mkdir(parents=True, exist_ok=True)


def save_png(path, tries=8):
    import os, time
    tmp = path.with_name(path.stem + "_tmp.png")
    fig.write_image(tmp, scale=4)
    for i in range(tries):
        try:
            os.replace(tmp, path)
            return
        except OSError:
            time.sleep(5)
    raise OSError(f"{path.name} is locked; new version left in {tmp.name}")


save_png(OUT / "Figure_2_Country_Distribution.png")
fig.write_image(OUT / "Figure_2_Country_Distribution.pdf")
fig.write_image(OUT / "Figure_2_Country_Distribution.svg")
fig.write_html(OUT / "Figure_2_Country_Distribution_interactive.html", include_plotlyjs=True)
pd.DataFrame({"Papers": counts.reindex(order), "Percent": (counts.reindex(order) / N * 100).round(1)}).to_csv(OUT / "Figure_2_data.csv", encoding="utf-8-sig")
print(N, len(order), counts.head(3).to_dict())
