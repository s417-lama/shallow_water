import glob
import re
import pandas as pd
import plotly.express as px

filenames = sorted(glob.glob("out3d/h_*.txt"))
file2t = lambda x: int(re.match(r'out3d/h_([0-9]+).txt', x).groups()[0])
df = pd.concat([pd.read_csv(filename, header=None, names=("x", "y", "h")).assign(t=file2t(filename))
                for filename in filenames], ignore_index=True)

fig = px.scatter_3d(df,
                    x="x",
                    y="y",
                    z="h",
                    color="h",
                    animation_frame="t",
                    range_z=[0, 3.5],
                    range_color=[2.5, 3.5],
                    color_continuous_scale=px.colors.sequential.ice)
fig.layout.updatemenus[0].buttons[0].args[1]["frame"]["duration"] = 100
fig.show()
