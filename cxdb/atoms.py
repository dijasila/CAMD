import plotly.graph_objects as go

fig = go.Figure(
        data=[
            go.Isosurface(
                x=X.flatten(),
                y=Y.flatten(),
                z=Z.flatten(),
                value=d.flatten(),
                isomin=-0.1,
                isomax=0.4,
                caps=dict(x_show=False, y_show=False)),
            go.Scatter3d(
                x=x,
                y=y,
                z=z,
            mode='markers')])
