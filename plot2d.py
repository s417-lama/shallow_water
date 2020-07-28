import glob
import re
import pandas as pd
import plotly.express as px

filenames = sorted(glob.glob("out2d/h_*.txt"))
file2t = lambda x: int(re.match(r'out2d/h_([0-9]+).txt', x).groups()[0])
df = pd.concat([pd.read_csv(filename, header=None, names=("x", "h")).assign(t=file2t(filename))
                for filename in filenames], ignore_index=True)

fig = px.scatter(df, x="x", y="h", animation_frame="t")
fig.layout.updatemenus[0].buttons[0].args[1]["frame"]["duration"] = 100
fig.show()
